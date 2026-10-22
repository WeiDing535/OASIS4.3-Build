"""Contextual PhaseBuild seed search; outputs only map-derived candidates."""
import argparse,json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
from phasebuild.data import Density,fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.seed_search import find_context_seeds

def main():
 p=argparse.ArgumentParser();p.add_argument('mtz');p.add_argument('--profile',required=True);p.add_argument('--output',required=True);p.add_argument('--count',type=int,default=512);a=p.parse_args();out=Path(a.output);out.mkdir(parents=True,exist_ok=False);start=time.time();d=Density(a.mtz);d.profile=DensityProfile.load(a.profile);sources=[fingerprint(p) for p in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
 seeds,stats=find_context_seeds(d,refine_count=a.count,progress=lambda m:print(m,flush=True))
 np.savez_compressed(out/'candidates.npz',seeds=np.array([s[1] for s in seeds]),seed_scores=np.array([s[0] for s in seeds]),seed_psi=np.array([s[2] for s in seeds]),grown=np.empty((0,4,3)),grown_lengths=np.empty(0,dtype=int),metadata=json.dumps({'input':d.source,'profile_training':d.profile.metadata}))
 report=dict(engine='phasebuild-contextual-tripeptide-search',parameters=vars(a),statistics=stats,seconds=time.time()-start,sources=sources,source_changed_during_run=any(fingerprint(s['path'])!=s for s in sources));(out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
