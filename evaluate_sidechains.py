import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.sidechain_evaluate import evaluate_sidechains


def main():
    p=argparse.ArgumentParser()
    p.add_argument('model');p.add_argument('--reference',required=True);p.add_argument('--output',required=True)
    args=p.parse_args()
    report=evaluate_sidechains(args.model,args.reference)
    Path(args.output).write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['residues','backbone_summary']},indent=2))


if __name__=='__main__':
    main()
