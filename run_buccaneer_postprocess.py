"""Explicit hybrid postprocessing; never relabelled as native PhaseBuild."""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.buccaneer_stages import postprocess,STAGES

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('mtz');p.add_argument('--model',required=True);p.add_argument('--sequence',required=True);p.add_argument('--output',required=True);p.add_argument('--stages',nargs='+',choices=STAGES,default=['sequence','rebuild']);p.add_argument('--cycles',type=int,default=1);p.add_argument('--build-semet',action='store_true');a=p.parse_args();print(json.dumps(postprocess(a.mtz,a.model,a.sequence,a.output,a.stages,a.cycles,build_semet=a.build_semet),indent=2))
if __name__=='__main__':main()
