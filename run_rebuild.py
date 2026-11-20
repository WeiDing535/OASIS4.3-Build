"""A bounded post-refinement pruning, loop and sidechain rebuilding stage."""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import gemmi
from phasebuild import __version__
from phasebuild.data import Density,fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.residue_profile import ResidueProfile
from phasebuild.rotamers import load_library
from phasebuild.model import read_backbone
from phasebuild.loops import sequence_bindings,build_short_loops
from phasebuild.sequencing import read_sequences,assign_sequence,annotate_structure
from phasebuild.rebuilding import prune_weak_fragments,transfer_properties
from phasebuild.sidechains import build_sidechains
from phasebuild.trace import write_backbone
from phasebuild.validation import geometry_report
from phasebuild.atom_validation import atom_geometry_report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mtz');p.add_argument('--model',required=True);p.add_argument('--sequence-report',required=True)
    p.add_argument('--sequence',required=True);p.add_argument('--profile',required=True)
    p.add_argument('--residue-profile',required=True);p.add_argument('--rotamers',required=True)
    p.add_argument('--output',required=True);p.add_argument('--loop-starts',type=int,default=48)
    p.add_argument('--seed',type=int,default=0)
    p.add_argument('--map-amplitude',default='FWT')
    p.add_argument('--map-phase',default='PHWT')
    args=p.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    sources=[fingerprint(f) for f in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
    density=Density(args.mtz,args.map_amplitude,args.map_phase);density.profile=DensityProfile.load(args.profile)
    profile=ResidueProfile.load(args.residue_profile);library=load_library(args.rotamers)
    records=read_sequences(args.sequence)
    previous=json.loads(Path(args.sequence_report).read_text())
    if previous['sequence_input']['sha256']!=fingerprint(args.sequence)['sha256']:
        raise ValueError('Sequence report belongs to another sequence')
    structure=gemmi.read_structure(args.model)
    bindings=sequence_bindings(structure,previous['sequence_assignment'],records)
    names=[c.name for c in structure[0]]
    pruning=prune_weak_fragments(structure,density)
    bindings=[b for name,b in zip(names,bindings) if name in {c.name for c in structure[0]}]
    structure.make_mmcif_document().write_file(str(out/'pruned.cif'))
    chains=read_backbone(out/'pruned.cif',density)
    chains,loops=build_short_loops(chains,density,records,bindings,profile,starts=args.loop_starts,seed=args.seed)
    write_backbone(out/'backbone.cif',chains,density)
    sequence=assign_sequence(chains,density,profile,records,seed=args.seed)
    annotate_structure(out/'backbone.cif',out/'assigned_backbone.cif',sequence,records)
    rebuilt=gemmi.read_structure(str(out/'assigned_backbone.cif'))
    transfer=transfer_properties(structure,rebuilt)
    sidechains=build_sidechains(rebuilt,density,library)
    geometry=atom_geometry_report(rebuilt,library)
    rebuilt.make_mmcif_document().write_file(str(out/'model.cif'))
    report={'version':__version__,'engine':'native-phasebuild-post-refinement-rebuild',
            'parameters':vars(args),'input':density.source,'model':fingerprint(args.model),
            'sequence_input':fingerprint(args.sequence),'sequence_report':fingerprint(args.sequence_report),
            'profile':fingerprint(args.profile),'residue_profile':fingerprint(args.residue_profile),
            'rotamers':fingerprint(args.rotamers),'sources':sources,'runner':fingerprint(__file__),
            'source_changed_during_run':any(fingerprint(s['path'])!=s for s in sources),
            'pruning':pruning,'loops':loops,'property_transfer':transfer,'sequence_assignment':sequence,
            'sidechains':sidechains,'geometry':geometry_report(chains,density),'atom_geometry':geometry,
            'seconds':time.perf_counter()-started,
            'scope':'Uses the supplied map and refined coordinates/B factors. No Rfree, target reference, '
                    'or Buccaneer output is used by this rebuilding stage.'}
    (out/'build_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ['pruning','loops','property_transfer','seconds']},indent=2))


if __name__=='__main__':
    main()
