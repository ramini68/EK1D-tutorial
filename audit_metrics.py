"""Grid-independent inventory/distribution audit of saved EK1D outputs.

No interpolation or averaging is applied to the retained concentration data.
CDFs are the exact integrals of the finite-volume piecewise-constant states.
Their differences are piecewise linear on the union of both grids, allowing
exact supremum and absolute-integral calculations. Unequal-mass CDF distance
normalised by initial Cu has units of length; it is NOT a Wasserstein metric.
The separately labelled unit-mass shape distance is the 1-D Wasserstein-1
distance between the normalised precipitate distributions.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from EK1D_tutorial import snapshot_prefix


def absolute_linear_integral(x, y):
    """Exact integral of absolute value of a piecewise-linear function."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.ndim != 1 or x.shape != y.shape or x.size < 2 or not np.isfinite(x).all() or not np.isfinite(y).all() or np.any(np.diff(x) <= 0):
        raise ValueError("Require matching vectors on a strictly increasing grid.")
    left, right, width = y[:-1], y[1:], np.diff(x)
    area = .5 * width * (np.abs(left) + np.abs(right))
    crossing = left * right < 0
    area[crossing] = width[crossing] * (
        left[crossing]**2 + right[crossing]**2
    ) / (2 * (np.abs(left[crossing]) + np.abs(right[crossing])))
    return float(np.sum(area))


def cdf_on(edges, density, positions):
    """Integrate piecewise-constant areal-inventory density onto positions."""
    edges, density, positions = map(lambda v:np.asarray(v,dtype=float),(edges,density,positions))
    if (edges.ndim != 1 or edges.size < 2 or density.shape != (edges.size-1,)
            or not np.isfinite(edges).all() or not np.isfinite(density).all()
            or not np.isfinite(positions).all() or np.any(np.diff(edges)<=0)
            or np.any(density<0) or np.any(positions<edges[0]) or np.any(positions>edges[-1])):
        raise ValueError('Require a finite increasing grid, nonnegative cell densities and in-domain positions.')
    cumulative = np.r_[0., np.cumsum(density * np.diff(edges))]
    return np.interp(positions, edges, cumulative)


def distribution_comparison(edges_a, density_a, edges_b, density_b, initial_cu):
    if not np.isfinite(initial_cu) or initial_cu <= 0 or not np.allclose(
        [edges_a[0], edges_a[-1]], [edges_b[0], edges_b[-1]], rtol=0, atol=1e-14
    ):
        raise ValueError("Positive initial Cu and matching physical domains required.")
    points = np.unique(np.r_[edges_a, edges_b])
    a, b = cdf_on(edges_a, density_a, points), cdf_on(edges_b, density_b, points)
    delta = a - b
    result = {
        "inventory_a_mol_m2": float(a[-1]),
        "inventory_b_mol_m2": float(b[-1]),
        "inventory_difference_relative_initial_Cu": float(abs(delta[-1]) / initial_cu),
        "CDF_sup_difference_relative_initial_Cu": float(np.max(abs(delta)) / initial_cu),
        "CDF_L1_over_initial_Cu_mm": 1000 * absolute_linear_integral(points, delta) / initial_cu,
        "normalised_shape_Wasserstein1_mm": None,
    }
    if a[-1] > 0 and b[-1] > 0:
        result["normalised_shape_Wasserstein1_mm"] = 1000 * absolute_linear_integral(
            points, a / a[-1] - b / b[-1]
        )
    return result


def read_run(folder):
    folder = Path(folder)
    diagnostics = json.loads((folder / "diagnostics.json").read_text())
    with np.load(folder / "states.npz") as archive:
        arrays = {key: archive[key] for key in archive.files}
    p = diagnostics["parameters"]
    if arrays['x'].shape != (p['N'],):
        raise ValueError('Saved centre count does not match N.')
    times=np.asarray(arrays['times_s'])
    if (times.ndim!=1 or len(times)<2 or not np.isfinite(times).all()
            or times[0]!=0 or times[-1]!=p['tend'] or np.any(np.diff(times)<=0)):
        raise ValueError('Saved times must start at zero, increase strictly and reach the declared final time.')
    if arrays['currents'].shape!=times.shape or not np.isfinite(arrays['currents']).all():
        raise ValueError('Saved instantaneous currents must be finite and match saved times.')
    if not np.isclose(arrays['currents'][-1],diagnostics['final_current_A_m2'],rtol=1e-12,atol=1e-14):
        raise ValueError('Final saved current and diagnostics disagree.')
    edges = np.linspace(0., p["L"], p["N"] + 1)
    if not np.allclose(arrays["x"], .5 * (edges[:-1] + edges[1:]), rtol=0, atol=1e-14):
        raise ValueError("Saved centres do not match the specified uniform FV grid.")
    final = snapshot_prefix(p['tend']) + '_'
    if (arrays[final+'c'].shape!=(5,p['N']) or arrays[final+'Cu_s'].shape!=(p['N'],)
            or arrays[final+'Cu_ppt'].shape!=(p['N'],) or arrays[final+'pH'].shape!=(p['N'],)):
        raise ValueError('Final state dimensions do not match the declared five-species grid.')
    fields = np.array([arrays[final + "c"][2], arrays[final + "Cu_s"], arrays[final + "Cu_ppt"]])
    if np.any(fields < 0) or not np.isfinite(fields).all():
        raise ValueError("Nonfinite or negative final copper state.")
    if 'phase_inventory_series' in arrays:
        series=arrays['phase_inventory_series']
        masses=p['theta']*np.sum(fields*np.diff(edges),axis=1)
        if series.shape!=(len(times),3) or not np.allclose(series[-1],masses,rtol=1e-12,atol=1e-14):
            raise ValueError('Saved phase-inventory history does not reconcile with the final spatial states.')
    if 'signed_charge_C_m2' in arrays:
        charge=arrays['signed_charge_C_m2']
        if charge.shape!=times.shape or not np.isclose(charge[-1],diagnostics['integrated_signed_current_C_m2'],rtol=1e-12,atol=1e-12):
            raise ValueError('Saved integrated charge and diagnostics disagree.')
    return dict(folder=str(folder), diagnostics=diagnostics, arrays=arrays,
                edges=edges, fields=fields, pH=arrays[final + "pH"])


def threshold_crossings(x, values, threshold=7.):
    hits = [float(x[i]) for i in np.flatnonzero(values == threshold)]
    for i in np.flatnonzero((values[:-1] - threshold) * (values[1:] - threshold) < 0):
        hits.append(float(x[i] + (threshold - values[i]) * (x[i+1] - x[i]) / (values[i+1] - values[i])))
    return sorted(set(hits))


def describe_run(folder):
    run = read_run(folder)
    d, arrays, fields = run["diagnostics"], run["arrays"], run["fields"]
    p = d["parameters"]
    masses = np.sum(fields * np.diff(run["edges"]), axis=1) * p["theta"]
    ppt = fields[2]
    return dict(
        folder=str(folder), N=p["N"], spatial_order=p.get("spatial_order", 1),
        time_order=p.get("time_order", 1), cfl=p["cfl"], final_time_h=p["tend"] / 3600,
        phase_inventory_mol_m2=dict(zip(("aqueous", "sorbed", "precipitated"), masses.tolist())),
        precipitate_raw_peak_mol_m3=float(ppt.max()),
        precipitate_peak_position_mm=float(arrays["x"][np.argmax(ppt)] * 1000) if ppt.max() > 0 else None,
        precipitate_centroid_mm=float(np.dot(ppt, arrays["x"]) / ppt.sum() * 1000) if ppt.sum() > 0 else None,
        pH7_front_positions_mm=[v * 1000 for v in threshold_crossings(arrays["x"], run["pH"])],
        final_current_A_m2=d["final_current_A_m2"],
        integrated_absolute_current_C_m2=d["faradaic_moles_each_electrode_per_m2"] * 96485.,
        max_relative_inventory_errors=d["max_relative_inventory_errors"],
        source_sha256=d.get('source_sha256',{}),
        recorded_species=d.get('species'),recorded_physical_constants=d.get('physical_constants'),
    )


def compare_runs(folder_a, folder_b, bin_width_m=.002, bin_offsets_m=None):
    if not np.isfinite(bin_width_m) or bin_width_m<=0:
        raise ValueError('Reporting-bin width must be finite and positive.')
    a, b = read_run(folder_a), read_run(folder_b)
    pa, pb = a["diagnostics"]["parameters"], b["diagnostics"]["parameters"]
    if pa["tend"] != pb["tend"]:
        raise ValueError("Final simulation times must match.")
    initial = a["diagnostics"]["initial_inventories_mol_m2"][0]
    if not np.isclose(initial, b["diagnostics"]["initial_inventories_mol_m2"][0], rtol=1e-12, atol=0):
        raise ValueError("Initial copper inventories differ.")
    profiles = {}
    for name, field_a, field_b in zip(("aqueous", "sorbed", "precipitated"), a["fields"], b["fields"]):
        profiles[name] = distribution_comparison(
            a["edges"], field_a * pa["theta"], b["edges"], field_b * pb["theta"], initial
        )
    da, db = describe_run(folder_a), describe_run(folder_b)
    centroid_delta = None
    if da["precipitate_centroid_mm"] is not None and db["precipitate_centroid_mm"] is not None:
        centroid_delta = abs(da["precipitate_centroid_mm"] - db["precipitate_centroid_mm"])
    front_delta = None
    if len(da["pH7_front_positions_mm"]) == len(db["pH7_front_positions_mm"]) == 1:
        front_delta = abs(da["pH7_front_positions_mm"][0] - db["pH7_front_positions_mm"][0])
    bins = []
    if bin_offsets_m is None:
        bin_offsets_m = np.arange(8) * bin_width_m / 8
    for offset in bin_offsets_m:
        if not 0 <= offset < bin_width_m:
            raise ValueError("Bin offsets must lie in [0, bin_width).")
        edges = np.unique(np.r_[0., np.arange(offset, pa["L"], bin_width_m), pa["L"]])
        ma = np.diff(cdf_on(a["edges"], a["fields"][2] * pa["theta"], edges))
        mb = np.diff(cdf_on(b["edges"], b["fields"][2] * pb["theta"], edges))
        bins.append(dict(offset_mm=float(offset * 1000),
                         profile_L1_relative_initial_Cu=float(np.sum(abs(ma - mb)) / initial)))
    return dict(
        run_a=str(folder_a), run_b=str(folder_b), phase_distributions=profiles,
        precipitate_centroid_difference_mm=centroid_delta,
        pH7_front_difference_mm=front_delta,
        precipitate_raw_peak_a_mol_m3=da["precipitate_raw_peak_mol_m3"],
        precipitate_raw_peak_b_mol_m3=db["precipitate_raw_peak_mol_m3"],
        final_current_relative_difference=abs(da["final_current_A_m2"] / db["final_current_A_m2"] - 1)
            if db["final_current_A_m2"] != 0 else None,
        integrated_absolute_current_relative_difference=abs(da["integrated_absolute_current_C_m2"] / db["integrated_absolute_current_C_m2"] - 1)
            if db["integrated_absolute_current_C_m2"] != 0 else None,
        reporting_bin_width_mm=bin_width_m * 1000,
        shifted_precipitate_bin_comparisons=bins,
        note="Grid differences are observed discrepancies, not estimates of error relative to an exact solution."
    )


def self_test():
    assert absolute_linear_integral([0, 1], [-1, 1]) == .5
    assert absolute_linear_integral([0, 1], [1, 2]) == 1.5
    identical = distribution_comparison([0, .5, 1], [2, 0], [0, .25, .5, 1], [2, 2, 0], 1)
    assert identical["CDF_L1_over_initial_Cu_mm"] == 0
    shifted = distribution_comparison([0, .5, 1], [2, 0], [0, .5, 1], [0, 2], 1)
    assert shifted["normalised_shape_Wasserstein1_mm"] == 500
    unequal = distribution_comparison([0, 1], [1], [0, 1], [2], 1)
    assert unequal["normalised_shape_Wasserstein1_mm"] == 0
    assert unequal["CDF_L1_over_initial_Cu_mm"] == 500
    print("audit_metrics: analytic integration and distribution tests passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("folders", nargs="*")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--out")
    args = parser.parse_args()
    if args.self_test:
        self_test()
    if args.folders:
        result = dict(runs=[describe_run(f) for f in args.folders],
                      comparisons=[compare_runs(a, b) for a, b in zip(args.folders[:-1], args.folders[1:])])
        rendered = json.dumps(result, indent=2)
        if args.out:
            Path(args.out).write_text(rendered + "\n")
        print(rendered)
