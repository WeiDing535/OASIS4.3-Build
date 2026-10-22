"""PhaseBuild-seeded Buccaneer/REFMAC composite method (external dependency)."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.hybrid_recycle import recycle

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--map',required=True);p.add_argument('--native',required=True);p.add_argument('--model',required=True);p.add_argument('--sequence',required=True);p.add_argument('--output',required=True);p.add_argument('--rounds',type=int,default=5);p.add_argument('--cycles',type=int,default=10);p.add_argument('--build-semet',action='store_true');p.add_argument('--runtime','--ccp4',dest='runtime',help='Bundled engine runtime');a=p.parse_args();r=recycle(a.map,a.native,a.model,a.sequence,a.output,a.rounds,a.cycles,ccp4=a.runtime,build_semet=a.build_semet,progress=lambda s:print(s,flush=True));print(json.dumps(r['primary_metric'],indent=2))
if __name__=='__main__':main()
