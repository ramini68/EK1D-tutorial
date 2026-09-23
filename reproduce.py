"""Reproduce the declared tutorial refinement matrix without overwriting runs.

Run from an unpacked release: python reproduce.py --output rerun --backend native
This can be computationally expensive. Existing result directories are refused;
use a different --output path for an independent rerun. No tuning/calibration.
"""
from pathlib import Path
import argparse
import json
import subprocess
import sys
from EK1D_tutorial import P,run,save
from assess_publication import hard_checks
from verify_saved_states import check as check_saved_states


CASES=[('N100_cfl02',100,.2,7.5),('N200_cfl02',200,.2,7.5),
       ('N200_cfl01',200,.1,3.75),('N400_cfl02',400,.2,7.5),
       ('N400_cfl01',400,.1,3.75)]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',default='rerun')
    parser.add_argument('--backend',choices=('native','python'),default='native')
    parser.add_argument('--case',choices=[c[0] for c in CASES])
    args=parser.parse_args()
    root=Path(args.output)
    criteria=json.loads(Path(__file__).with_name('verification_criteria.json').read_text())
    chosen=[c for c in CASES if args.case is None or c[0]==args.case]
    for name,_,_,_ in chosen:
        if (root/name).exists():
            raise FileExistsError(f'Refusing to overwrite {root/name}; choose a new --output.')
    for name,n,cfl,dt in chosen:
        print(f'Starting {name}: 36 simulated hours',flush=True)
        result=run(P(N=n,cfl=cfl,dt=dt,spatial_order=2,time_order=2,output_dt=60.),backend=args.backend)
        save(result,root/name)
        checks=hard_checks(root/name,criteria['hard'])
        (root/name/'run_checks.json').write_text(json.dumps(checks,indent=2,allow_nan=False)+'\n')
        if not checks['passed']:
            raise RuntimeError(f'{name} failed numerical checks; inspect run_checks.json.')
        state_checks=check_saved_states(root/name)
        (root/name/'saved_state_checks.json').write_text(json.dumps(state_checks,indent=2,allow_nan=False)+'\n')
        print(json.dumps(result['diagnostics'],indent=2),flush=True)
    if args.case is None:
        subprocess.run([sys.executable,str(Path(__file__).with_name('assess_publication.py')),
            '--coarse',str(root/'N200_cfl02'),'--fine',str(root/'N400_cfl02'),
            '--time-refined',str(root/'N400_cfl01'),'--extra-coarse',str(root/'N100_cfl02'),
            '--extra-time-refined',str(root/'N200_cfl01'),
            '--out',str(root/'publication_assessment.json')],check=True)
        from plot_article_figures import render_article_figures
        render_article_figures(root,root/'article_figures',prefix='')
        print(f'Article figures: {root / "article_figures"}. Read the assessment scope before use.',flush=True)


if __name__=='__main__':
    main()
