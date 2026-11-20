"""Execute the CCP4 buccaneer_pipeline building recipe as a pipeline stage."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild import ccp4_pipeline


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--map',required=True)
    p.add_argument('--native',required=True)
    p.add_argument('--sequence',required=True)
    p.add_argument('--runtime',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--rounds',type=int,default=5)
    p.add_argument('--free-value',type=int,default=0)
    a=p.parse_args()
    report=ccp4_pipeline.run(a.native,a.map,a.sequence,a.runtime,a.output,
        cycles=a.rounds,free_value=a.free_value)
    print(json.dumps({'status':report['status'],'primary_metric':report['primary_metric']}))


if __name__=='__main__':main()
