"""Composite PhaseBuild/Buccaneer recycling with a fixed experimental HL prior.

External Buccaneer performs the reconstruction. REFMAC retains the original
experimental phase likelihood; missing high-resolution priors remain missing.
"""
import argparse,json,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild import hybrid_recycle
from phasebuild.fixed_phase_prior import make_fixed_prior_refiner

def main(argv=None):
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--map',required=True);p.add_argument('--native',required=True);p.add_argument('--model',required=True);p.add_argument('--sequence',required=True);p.add_argument('--phase-prior',required=True);p.add_argument('--phase-prior-sha256');p.add_argument('--output',required=True);p.add_argument('--rounds',type=int,default=5);p.add_argument('--cycles',type=int,default=10);p.add_argument('--free-value',type=int,default=0);p.add_argument('--build-semet',action='store_true');p.add_argument('--runtime','--ccp4',dest='runtime',help='Bundled engine runtime; auto-detected when omitted');p.add_argument('--existing-first-build',help='Execution recovery: reuse a verified completed first3cycle construction');a=p.parse_args(argv)
 hybrid_recycle.refmac_refine=make_fixed_prior_refiner(a.phase_prior,expected_sha256=a.phase_prior_sha256,free_value=a.free_value)
 report=hybrid_recycle.recycle(a.map,a.native,a.model,a.sequence,a.output,a.rounds,a.cycles,ccp4=a.runtime,build_semet=a.build_semet,progress=lambda s:print(s,flush=True),existing_first_build=a.existing_first_build);print(json.dumps(report['primary_metric'],indent=2))
if __name__=='__main__':main()
