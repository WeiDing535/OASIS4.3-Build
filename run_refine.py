"""Refine/evaluate a PhaseBuild model using native F/SIGF and existing free flags."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.refinement import refmac_refine


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mtz');p.add_argument('--model',required=True);p.add_argument('--output',required=True)
    p.add_argument('--cycles',type=int,default=10)
    p.add_argument('--runtime','--ccp4',dest='ccp4',help='Bundled runtime or existing CCP4 directory')
    p.add_argument('--amplitude',default='F');p.add_argument('--sigma',default='SIGF')
    p.add_argument('--free-column',default='FreeR_flag');p.add_argument('--free-value',type=int,default=0)
    args=p.parse_args()
    report=refmac_refine(args.mtz,args.model,args.output,args.cycles,args.ccp4,
                         args.amplitude,args.sigma,args.free_column,args.free_value)
    print(json.dumps(report['statistics'],indent=2))


if __name__=='__main__':
    main()
