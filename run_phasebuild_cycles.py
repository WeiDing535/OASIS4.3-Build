import argparse
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.cycling import run_cycles


def main():
    p=argparse.ArgumentParser(description='Fixed-budget refinement/rebuilding experiment')
    p.add_argument('--native',required=True);p.add_argument('--sequence',required=True)
    p.add_argument('--start-refinement',required=True);p.add_argument('--sequence-report',required=True)
    p.add_argument('--original-map',required=True);p.add_argument('--profile',required=True)
    p.add_argument('--residue-profile',required=True);p.add_argument('--rotamers',required=True)
    p.add_argument('--mode',choices=['refine-only','fixed-map','updated-map'],default='updated-map')
    p.add_argument('--rebuild-stage',choices=['local','extension'],default='local')
    p.add_argument('--additional-rounds',type=int,default=3);p.add_argument('--initial-cycles',type=int,default=20)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    run_cycles(a.native,a.sequence,a.start_refinement,a.sequence_report,a.original_map,a.profile,
               a.residue_profile,a.rotamers,a.output,a.mode,a.additional_rounds,a.initial_cycles,
               progress=lambda message:print(message,flush=True),rebuild_stage=a.rebuild_stage)


if __name__=='__main__':
    main()
