"""Evidence-bound assessment of complete EK1D runs, with explicit failed checks."""
from pathlib import Path
import argparse
import json
import numpy as np
from audit_metrics import compare_runs, describe_run, read_run


# The retained CLI matrix used the first two audited driver/wrapper pairs.
# The additional packaged driver changes documentation and generic plotting
# only; SOURCE_UPDATE.json records the reversible patch and unchanged solver
# AST. Historical source hashes and missing species/constants stay unmodified.
AUDITED_SOURCE_PAIRS={
    ('8b87254fb69370e5cf39f2761be89fb6ea92759bdff3714cdb9e102563e09c59',
     '83aee1890774107bbfeab10f7bfd7982e30f29af39fb19cb4a7a6492fc4164ad'),
    ('1ac126ca67736fa233fcac9622a84d2cf6dbb965199228701f3c363e73c4720d',
     '5487be0b51d976cbd03ca5a0f2b4fa55ecc31005eaee5ecf8a41d7f60605c70c'),
    ('873cc7eb62b623e4d1a734784cdf502b2dea18ae6fb833d952f42e5829dbbd0b',
     '83aee1890774107bbfeab10f7bfd7982e30f29af39fb19cb4a7a6492fc4164ad'),
}
AUDITED_EQUATION_HASHES={
    'chemistry.py':'7308ac360f1ac02fec9495253d25fa770ee8d7c4498bfa686521ac8f01fa8125',
    'native_kernel.cpp':'a0bf0ee427b6ac38dd76759632dc44bb1a38fd58009e90f2c907daafd81f8923',
}
CANONICAL_FIXED_PARAMETERS=dict(L=.1,theta=.45,tau=1.6,V=20.,keo=1.2e-9,
    Kw=1e-14,Ksp=2.2e-20,pH50=5.5,nedge=2.,pH0=5.,Cu_tot0=5.,Na0=10.,
    tend=129600.,spatial_order=2,time_order=2,algebraic_chemistry=True)
CANONICAL_SPECIES=dict(names=['H','OH','Cu','Na','NO3'],charges=[1.,-1.,2.,1.,-1.],
    D0_m2_s=[9.31e-9,5.27e-9,.714e-9,1.33e-9,1.90e-9])
CANONICAL_CONSTANTS=dict(F_C_mol=96485.,R_J_mol_K=8.314,T_K=298.15)


def validate_pair(folder_a,folder_b,kind):
    """Reject vacuous or confounded refinement comparisons before scoring."""
    if kind not in ('spatial','temporal'):
        raise ValueError('Comparison kind must be spatial or temporal.')
    if Path(folder_a).resolve()==Path(folder_b).resolve():
        raise ValueError('Refinement requires two distinct result directories.')
    a,b=read_run(folder_a),read_run(folder_b)
    da,db=a['diagnostics'],b['diagnostics'];pa,pb=da['parameters'],db['parameters']
    allowed={'N'} if kind=='spatial' else {'cfl','dt'}
    if set(pa)!=set(pb):
        raise ValueError('Parameter records have different fields; resolve provenance before comparing.')
    differences={k:[pa[k],pb[k]] for k in pa if pa[k]!=pb[k]}
    confounded=set(differences)-allowed
    if confounded:
        raise ValueError('Refinement changes non-target parameters: '+', '.join(sorted(confounded)))
    if kind=='spatial' and not pb['N']>pa['N']:
        raise ValueError('Spatial refinement requires a strictly larger fine-grid N.')
    if kind=='temporal' and not (pb['dt']<pa['dt'] and pb['cfl']<pa['cfl']):
        raise ValueError('This temporal test requires both maximum dt and CFL to decrease at fixed N.')
    if da.get('backend')!=db.get('backend'):
        raise ValueError('Use one backend for a refinement pair; backend equivalence is a separate check.')
    if not np.array_equal(a['arrays']['times_s'],b['arrays']['times_s']):
        raise ValueError('Refinement comparisons require identical output times.')
    hashes=[d.get('source_sha256',{}) for d in (da,db)]
    for name in ('chemistry.py','native_kernel.cpp'):
        if not hashes[0].get(name) or hashes[0].get(name)!=hashes[1].get(name):
            raise ValueError(f'Refinement requires matching {name} equation hashes.')
    missing={str(folder):[k for k in ('species','physical_constants') if d.get(k) is None]
             for folder,d in ((folder_a,da),(folder_b,db))}
    missing={k:v for k,v in missing.items() if v}
    source_pairs=[(h.get('EK1D_tutorial.py'),h.get('native_backend.py')) for h in hashes]
    compatibility_note='Matching recorded species/constants and source identities.'
    if missing or source_pairs[0]!=source_pairs[1]:
        for d,p,h,pair in zip((da,db),(pa,pb),hashes,source_pairs):
            if pair not in AUDITED_SOURCE_PAIRS or any(h.get(k)!=v for k,v in AUDITED_EQUATION_HASHES.items()):
                raise ValueError('Unknown missing-metadata or mixed-version source; an explicit source audit is required.')
            if d.get('backend')!='native' or any(p.get(k)!=v for k,v in CANONICAL_FIXED_PARAMETERS.items()):
                raise ValueError('Legacy compatibility evidence applies only to the documented canonical native CLI matrix.')
            for key,expected in (('species',CANONICAL_SPECIES),('physical_constants',CANONICAL_CONSTANTS)):
                if d.get(key) is not None and d[key]!=expected:
                    raise ValueError(f'Recorded {key} disagree with the audited canonical CLI configuration.')
        compatibility_note=('Explicitly audited driver/wrapper hash pairs; identical native-kernel and chemistry '
            'hashes plus the documented canonical-CLI launches support compatibility. Missing legacy '
            'species/constants metadata remain unavailable, not reconstructed stored measurements.')
    else:
        if source_pairs[0][0] is None or source_pairs[0][1] is None:
            raise ValueError('Driver and wrapper source hashes are required.')
        for key in ('species','physical_constants'):
            if da[key]!=db[key]:
                raise ValueError(f'Refinement changes recorded {key}.')
    return dict(kind=kind,changed_controls=differences,missing_recorded_metadata=missing,
                source_compatibility=compatibility_note,
                run_a_source_sha256=hashes[0],run_b_source_sha256=hashes[1])


def temporal_series_comparison(folder_a,folder_b):
    aa,bb=read_run(folder_a),read_run(folder_b)
    a,b=aa['arrays'],bb['arrays']
    if not np.array_equal(a['times_s'],b['times_s']):
        raise ValueError('Compare identical output times; no interpolation of transient spikes.')
    t=a['times_s'];initial=bb['diagnostics']['initial_inventories_mol_m2'][0]
    ja,jb=a['currents'],b['currents']
    qa,qb=a['signed_charge_C_m2'],b['signed_charge_C_m2']
    ia,ib=a['phase_inventory_series'],b['phase_inventory_series']
    if (initial<=0 or qb[-1]<=0 or jb[-1]<=0
            or np.min(a['accepted_stage_current_envelopes'])<=0
            or np.min(b['accepted_stage_current_envelopes'])<=0
            or np.min(ja)<=0 or np.min(jb)<=0):
        raise ValueError('These relative-current screening metrics require the positive-current tutorial case.')
    k=int(np.argmax(abs(ja-jb)))
    ma,mb=np.diff(qa)/np.diff(t),np.diff(qb)/np.diff(t)
    return dict(
        final_current_relative=float(abs(ja[-1]/jb[-1]-1)),
        integrated_current_relative=float(abs(qa[-1]/qb[-1]-1)),
        max_final_phase_inventory_difference_over_initial_Cu=float(np.max(abs(ia[-1]-ib[-1]))/initial),
        max_phase_history_difference_over_initial_Cu=float(np.max(abs(ia-ib))/initial),
        max_charge_history_difference_over_final_charge=float(np.max(abs(qa-qb))/abs(qb[-1])),
        sampled_current_history_L1_relative=float(np.trapezoid(abs(ja-jb),t)/np.trapezoid(abs(jb),t)),
        interval_mean_current_history_L1_relative=float(np.sum(abs(ma-mb)*np.diff(t))/abs(qb[-1])),
        maximum_sampled_current_difference_A_m2=float(abs(ja[k]-jb[k])),
        maximum_current_difference_time_h=float(t[k]/3600),
        paired_currents_at_max_difference_A_m2=[float(ja[k]),float(jb[k])],
        accepted_stage_current_extrema_a_A_m2=[float(a['accepted_stage_current_envelopes'][:,0].min()),float(a['accepted_stage_current_envelopes'][:,1].max())],
        accepted_stage_current_extrema_b_A_m2=[float(b['accepted_stage_current_envelopes'][:,0].min()),float(b['accepted_stage_current_envelopes'][:,1].max())])


def hard_checks(folder,limits):
    run=read_run(folder);d=run['diagnostics'];arrays=run['arrays']
    observed={
        'relative_inventory_error':max(d['max_relative_inventory_errors'].values()),
        'normalised_charge_residual':d['max_normalized_charge'],
        'relative_equilibrium_residual':max(d['max_chemistry_residuals'].values()),
        'current_residual_A_m2':d['max_current_residual_A_m2'],
        'voltage_residual_V':d['max_voltage_residual_V']}
    conditions={k:bool(v<=limits[k]) for k,v in observed.items()}
    conditions['finite_saved_states']=all(bool(np.isfinite(v).all()) for v in arrays.values())
    conditions['positive_transport']=d['minimum_aqueous_mol_m3']>=0
    conditions['complete_36h']=arrays['times_s'][-1]==d['parameters']['tend']==129600.
    conditions['nonempty_integration']=d['accepted_steps']>0
    conditions['phase_inventory_reconciliation']=bool(np.max(abs(arrays['phase_inventory_series'].sum(axis=1)-arrays['inventory_series'][:,0]))<1e-12)
    conditions['positive_current_this_case']=bool(arrays['accepted_stage_current_envelopes'][:,0].min()>0 and arrays['currents'].min()>0)
    conditions['charge_faradaic_identity']=bool(abs(d['integrated_signed_current_C_m2']/96485.-d['faradaic_moles_each_electrode_per_m2'])<1e-10)
    return dict(folder=str(folder),passed=all(conditions.values()),observed=observed,conditions=conditions)


def compare_case(a,b,kind,criteria):
    design=validate_pair(a,b,kind)
    distribution=compare_runs(a,b)
    values=temporal_series_comparison(a,b)
    values['pH7_front_difference_mm']=distribution['pH7_front_difference_mm']
    checks={k:dict(value=values[k],limit=v,passed=values[k] is not None and values[k]<=v)
            for k,v in criteria[kind].items()}
    return dict(run_a=str(a),run_b=str(b),kind=kind,passed=all(v['passed'] for v in checks.values()),
                comparison_design=design,checks=checks,series=values,distributions=distribution)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--coarse',required=True)
    parser.add_argument('--fine',required=True)
    parser.add_argument('--time-refined',required=True)
    parser.add_argument('--extra-coarse')
    parser.add_argument('--extra-time-refined',help='Optional temporal refinement of --coarse, also included in hard checks and the overall gate.')
    parser.add_argument('--out',default='publication_assessment.json')
    args=parser.parse_args(argv)
    criteria=json.loads(Path(__file__).with_name('verification_criteria.json').read_text())
    folders=list(dict.fromkeys([args.coarse,args.fine,args.time_refined]
        +([args.extra_coarse] if args.extra_coarse else [])
        +([args.extra_time_refined] if args.extra_time_refined else [])))
    result=dict(criteria=criteria,runs=[describe_run(f) for f in folders],
        hard_checks=[hard_checks(f,criteria['hard']) for f in folders],
        spatial=compare_case(args.coarse,args.fine,'spatial',criteria),
        temporal=compare_case(args.fine,args.time_refined,'temporal',criteria))
    if args.extra_coarse:
        result['earlier_spatial']=compare_case(args.extra_coarse,args.coarse,'spatial',criteria)
    if args.extra_time_refined:
        result['earlier_temporal']=compare_case(args.coarse,args.extra_time_refined,'temporal',criteria)
    comparison_keys=[k for k in ('spatial','temporal','earlier_spatial','earlier_temporal') if k in result]
    result['gated_comparisons']=comparison_keys
    result['all_declared_checks_pass']=all(v['passed'] for v in result['hard_checks']) and all(result[k]['passed'] for k in comparison_keys)
    Path(args.out).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:result[k] for k in ['all_declared_checks_pass','hard_checks']},indent=2))
    for kind in comparison_keys:
        print(kind,json.dumps(result[kind]['checks'],indent=2))
    return 0 if result['all_declared_checks_pass'] else 1


if __name__=='__main__':
    raise SystemExit(main())
