"""Independent short growth from cached seeds; defer inter-fragment pruning."""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import Density,fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.weak_trace import grow_weak


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mtz')
    parser.add_argument('--profile',required=True)
    parser.add_argument('--candidate-archive',required=True)
    parser.add_argument('--output',required=True)
    parser.add_argument('--maximum-length',type=int,default=12)
    parser.add_argument('--beam-width',type=int,default=2)
    parser.add_argument('--seed-count',type=int,default=512)
    args=parser.parse_args()
    if args.maximum_length<3 or args.beam_width<1 or args.seed_count<1:
        parser.error('Length must be >=3; beam width and seed count must be positive')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    started=time.time();density=Density(args.mtz)
    density.profile=DensityProfile.load(args.profile)
    sources=[fingerprint(p) for p in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
    with np.load(args.candidate_archive,allow_pickle=False) as archive:
        metadata=json.loads(str(archive['metadata']))
        if metadata['input']['sha256']!=density.source['sha256']:
            raise ValueError('Map and cached candidate fingerprints differ')
        seeds=list(zip(archive['seed_scores'],archive['seeds'],archive['seed_psi']))[:args.seed_count]
    seeds=sorted(seeds,key=lambda item:float(density.score(item[1]).mean()),reverse=True)
    grown=[];rejected=0
    for i,(score,xyz,psi) in enumerate(seeds):
        values=density.score(xyz)
        if values.mean()<2.5 or values.min()<1.5:
            rejected+=1;continue
        chain=grow_weak(density,xyz,args.maximum_length,args.beam_width,used=())
        chain=grow_weak(density,chain,args.maximum_length,args.beam_width,backward=True,used=())
        grown.append(chain)
        print('Independent seed %d/%d: %d residues'%(i+1,len(seeds),len(chain)),flush=True)
    np.savez_compressed(out/'candidates.npz',seeds=np.asarray([s[1] for s in seeds]),
                        seed_scores=np.asarray([s[0] for s in seeds]),seed_psi=np.asarray([s[2] for s in seeds]),
                        grown=np.concatenate(grown) if grown else np.empty((0,4,3)),
                        grown_lengths=np.asarray(list(map(len,grown))),
                        metadata=json.dumps({'input':density.source,'profile_training':density.profile.metadata}))
    report={'engine':'phasebuild-independent-short-growth','parameters':vars(args),
            'input':density.source,'candidate_archive':fingerprint(args.candidate_archive),
            'profile':fingerprint(args.profile),'input_seeds':len(seeds),'grown_fragments':len(grown),
            'grown_residues':sum(map(len,grown)),'grown_lengths':list(map(len,grown)),
            'density_rejected_seeds':rejected,'sources':sources,
            'source_changed_during_run':any(fingerprint(s['path'])!=s for s in sources),
            'seconds':time.time()-started,
            'scope':'No inter-fragment used-coordinate exclusion during growth; final global assembly/pruning required. '
                    'No reference, Buccaneer model, or sequence labels are accepted.'}
    (out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
