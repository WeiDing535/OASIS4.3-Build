"""Experimental candidate selection with frozen unrelated-model discrimination."""
import argparse
import json
import sys
import time
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.data import Density,fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.backbone_discriminant import BackboneDiscriminant
from phasebuild.discriminant_assembly import assemble_discriminant_fragments
from phasebuild.trace import write_backbone
from phasebuild.validation import geometry_report
from phasebuild.weak_trace import supported_intervals


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mtz')
    parser.add_argument('--profile',required=True)
    parser.add_argument('--discriminant',required=True)
    parser.add_argument('--candidate-archive',required=True,action='append')
    parser.add_argument('--output',required=True)
    parser.add_argument('--capacity',required=True,type=int)
    parser.add_argument('--maximum-length',default=40,type=int)
    parser.add_argument('--beam-width',default=3,type=int)
    args=parser.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    started=time.time();density=Density(args.mtz)
    density.profile=DensityProfile.load(args.profile)
    classifier=BackboneDiscriminant.load(args.discriminant)
    root=Path(__file__).resolve().parent
    sources=[fingerprint(root/'src/phasebuild'/name) for name in
             ['discriminant_assembly.py','backbone_discriminant.py','fragment_join.py','data.py',
              'profile.py','geometry.py','selection.py','weak_trace.py','trace.py','validation.py']]
    fragments=[];archives=[];seed_count=0;grown_count=0
    for path in args.candidate_archive:
        with np.load(path,allow_pickle=False) as z:
            metadata=json.loads(str(z['metadata']))
            if metadata['input']['sha256']!=density.source['sha256']:
                raise ValueError('Candidate archive and map fingerprints differ')
            for xyz in z['seeds']:
                scores=density.score(xyz)
                if scores.mean()>=2.5 and scores.min()>=1.5:
                    fragments.append(xyz);seed_count+=1
            offset=0
            for length in z['grown_lengths']:
                chain=z['grown'][offset:offset+int(length)];offset+=int(length)
                for a,b in supported_intervals(density.score(chain)):
                    fragments.append(chain[a:b]);grown_count+=1
            archives.append(fingerprint(path))
    chains,report=assemble_discriminant_fragments(density,fragments,classifier,args.capacity,
                                                 maximum_length=args.maximum_length,
                                                 beam_width=args.beam_width,
                                                 progress=lambda s:print(s,flush=True))
    write_backbone(out/'backbone.cif',chains,density)
    report.update(engine='phasebuild-classifier-before-conflict-assembly',parameters=vars(args),
                  input=density.source,candidate_archives=archives,
                  profile=fingerprint(args.profile),discriminant=fingerprint(args.discriminant),
                  discriminant_training=classifier.metadata,
                  features='Backbone-frame density patches, normalized contrast and its square, '
                           'one-quarter raw patch, local mean/std, N/CA/C/O density, minimum and mean atom density',
                  seed_candidates=seed_count,grown_pieces=grown_count,residues=sum(map(len,chains)),
                  sources=sources,source_changed_during_run=any(fingerprint(s['path'])!=s for s in sources),
                  seconds=time.time()-started,geometry=geometry_report(chains,density))
    (out/'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['geometry','sources','discriminant_training']}),flush=True)


if __name__=='__main__':
    main()
