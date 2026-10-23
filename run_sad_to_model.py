"""Unified SAD -> density modification -> model building, with audited restart."""
import argparse
import json
from pathlib import Path
import sys
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_thread_key]='1'
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['sad','native','sequence','heavy','reference','fpp','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--runtime');p.add_argument('--labels',default='')
    for name,default in [('copies',1),('phase-cycles',5),('dm-cycles',5),('rounds',5),('cycles',10),('free-value',0),('seed',20261023)]:
        p.add_argument('--'+name,type=int,default=default)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--check',action='store_true',help='Read-only preflight and plan; create no output')
    args=vars(p.parse_args(argv));resume=args.pop('resume');check=args.pop('check')
    from phasebuild.sad_pipeline import validate,stage_plan,run_pipeline
    cfg=validate(args)
    if check:
        if Path(cfg['output']).exists() and not resume:p.error('Output exists; choose a new directory')
        print(json.dumps(dict(config=cfg,stages=stage_plan(cfg),output_created=False),indent=2));return 0
    report=run_pipeline(cfg,resume=resume)
    print(json.dumps({'status':report['status'],'primary_model':report['primary_model']},indent=2));return 0


if __name__=='__main__':raise SystemExit(main())
