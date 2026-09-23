"""Recompute physical/algebraic checks from saved arrays, not reported maxima.

Checks the canonical 36-hour tutorial. This is a data consistency audit, not
experimental validation or an additional spatial convergence level.
"""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
from EK1D_tutorial import P, grid, electrical, snapshot_prefix


def check(folder):
    folder = Path(folder)
    d = json.loads((folder/'diagnostics.json').read_text())
    p = P(**d['parameters'])
    dx, x, diffusivity = grid(p)
    rows = []
    with np.load(folder/'states.npz') as a:
        np.testing.assert_allclose(a['x'], x, rtol=0, atol=1e-14)
        assert a['times_s'][0] == 0 and a['times_s'][-1] == p.tend == 129600.
        assert np.all(np.diff(a['times_s']) > 0)
        assert all(np.isfinite(a[k]).all() for k in a.files if a[k].dtype.kind in 'biuf')
        phases = a['phase_inventory_series']
        assert phases.min() >= 0
        np.testing.assert_allclose(phases.sum(axis=1), p.theta*p.L*p.Cu_tot0, rtol=1e-8)
        np.testing.assert_allclose(phases.sum(axis=1), a['inventory_series'][:, 0], rtol=1e-12)
        np.testing.assert_allclose(a['signed_charge_C_m2']/96485., a['faradaic_moles_series'], rtol=1e-10, atol=1e-12)
        assert a['currents'].min() > 0
        times = a['snapshot_times_s'] if 'snapshot_times_s' in a else (0.,21600.,86400.,129600.)
        for t in times:
            key = snapshot_prefix(t)
            c, cs, cp = (a[key+'_'+s] for s in ('c', 'Cu_s', 'Cu_ppt'))
            assert c.min() >= 0 and cs.min() >= 0 and cp.min() >= 0
            h, oh, cu = c[:3]/1000.
            water = float(np.max(abs(h*oh/p.Kw-1)))
            ratio = (10**(-p.pH50)/h)**p.nedge
            total = c[2]+cs+cp
            sorption = float(np.max(abs(cs-ratio*c[2])/np.maximum(total,1e-15)))
            saturation = cu*oh**2/p.Ksp
            mineral = float(np.max(np.where(cp > 1e-12, abs(saturation-1), np.maximum(saturation-1,0))))
            charge = float(np.max(abs(p.z@c)/np.sum(abs(p.z[:,None]*c),axis=0)))
            assert max(water,sorption,mineral) < 1e-9
            assert charge < 1e-5
            np.testing.assert_allclose(a[key+'_pH'], -np.log10(h), rtol=0, atol=1e-12)
            # Python face closure is recomputed from C++-generated aqueous fields.
            el = electrical(c,dx,diffusivity,p)
            np.testing.assert_allclose(a[key+'_j'],el['j'],rtol=1e-10,atol=1e-10)
            np.testing.assert_allclose(a[key+'_E'],el['E'],rtol=1e-9,atol=1e-8)
            np.testing.assert_allclose(a[key+'_phi'],el['phi'],rtol=1e-9,atol=1e-9)
            assert el['current_residual'] < 1e-9 and el['voltage_residual'] < 1e-9
            idx = np.flatnonzero(a['times_s'] == t)
            assert len(idx) == 1
            pools = p.theta*dx*np.array([c[2].sum(),cs.sum(),cp.sum()])
            np.testing.assert_allclose(pools,phases[idx[0]],rtol=1e-12,atol=1e-13)
            np.testing.assert_allclose(a['currents'][idx[0]],el['j'],rtol=1e-10,atol=1e-10)
            rows.append(dict(time_s=float(t),water_residual=water,sorption_residual=sorption,
                             mineral_residual=mineral,normalized_charge=charge))
    return dict(folder=str(folder),passed=True,snapshots=rows,
                states_sha256=hashlib.sha256((folder/'states.npz').read_bytes()).hexdigest())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folders',nargs='+')
    parser.add_argument('--report',required=True)
    args=parser.parse_args()
    rows=[]
    for folder in args.folders:
        try:
            rows.append(check(folder))
        except Exception as exc:
            rows.append(dict(folder=folder,passed=False,error=f'{type(exc).__name__}: {exc}'))
    report=dict(passed=all(row['passed'] for row in rows),runs=rows,
                scope='Saved-state verification; not experimental validation or pointwise convergence.',
                verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(args.report).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
