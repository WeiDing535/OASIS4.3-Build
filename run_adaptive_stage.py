"""Run an existing construction stage with adaptive, batched sequencing.

Only the sequence permutation implementation changes; all stage arguments and
geometry construction routines are shared with the existing command-line tools.
"""
import argparse,runpy,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'src'))
def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('stage',choices=['improve','complete','rebuild'])
 args,rest=parser.parse_known_args()
 import phasebuild.sequencing as sequencing
 from phasebuild.fast_sequencing import assign_sequence
 sequencing.assign_sequence=assign_sequence
 script=ROOT/('run_'+args.stage+'.py')
 sys.argv=[str(script)]+rest
 runpy.run_path(str(script),run_name='__main__')
if __name__=='__main__':main()
