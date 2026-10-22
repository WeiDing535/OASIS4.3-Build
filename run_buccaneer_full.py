import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.benchmark import buccaneer_full


def main():
    p=argparse.ArgumentParser(description='Full sequence and sidechain Buccaneer control')
    p.add_argument('mtz');p.add_argument('--sequence',required=True);p.add_argument('--output',required=True)
    p.add_argument('--cycles',type=int,default=3);p.add_argument('--ccp4',default='/Applications/ccp4-9')
    p.add_argument('--correlation-mode',action='store_true')
    args=p.parse_args()
    report=buccaneer_full(args.mtz,args.sequence,args.output,args.cycles,args.ccp4,correlation=args.correlation_mode)
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
