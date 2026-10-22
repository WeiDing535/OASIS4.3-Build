"""Autonomous short-loop hypotheses and density-supported sidechains."""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import gemmi
from phasebuild import __version__
from phasebuild.data import Density,fingerprint
from phasebuild.model import read_backbone
from phasebuild.profile import DensityProfile
from phasebuild.residue_profile import ResidueProfile
from phasebuild.sequencing import read_sequences,assign_sequence,annotate_structure
from phasebuild.rotamers import load_library
from phasebuild.sidechains import build_sidechains
from phasebuild.trace import write_backbone
from phasebuild.validation import geometry_report
from phasebuild.atom_validation import atom_geometry_report


def main():
    p=argparse.ArgumentParser()
    p.add_argument('mtz')
    p.add_argument('--model',required=True)
    p.add_argument('--profile',required=True)
    p.add_argument('--rotamers',required=True)
    p.add_argument('--sequence',required=True)
    p.add_argument('--residue-profile',required=True)
    p.add_argument('--sequence-report',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--skip-loops',action='store_true')
    p.add_argument('--loop-starts',type=int,default=48)
    p.add_argument('--seed',type=int,default=0)
    p.add_argument('--map-amplitude',default='FWT')
    p.add_argument('--map-phase',default='PHWT')
    args=p.parse_args()
    out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    sources=[fingerprint(f) for f in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
    density=Density(args.mtz,args.map_amplitude,args.map_phase);density.profile=DensityProfile.load(args.profile)
    target=ResidueProfile.load(args.residue_profile)
    records=read_sequences(args.sequence)
    chains=read_backbone(args.model,density)
    loop_report=None
    if not args.skip_loops:
        from phasebuild.loops import sequence_bindings,build_short_loops
        previous=json.loads(Path(args.sequence_report).read_text())
        if previous['sequence_input']['sha256']!=fingerprint(args.sequence)['sha256']:
            raise ValueError('Sequence report belongs to another sequence')
        bindings=sequence_bindings(gemmi.read_structure(args.model),previous['sequence_assignment'],records)
        chains,loop_report=build_short_loops(chains,density,records,bindings,target,starts=args.loop_starts,seed=args.seed,
                                            progress=lambda s:print(s,flush=True))
    write_backbone(out/'backbone.cif',chains,density)
    sequence=assign_sequence(chains,density,target,records,seed=args.seed)
    annotate_structure(out/'backbone.cif',out/'assigned_backbone.cif',sequence,records)
    structure=gemmi.read_structure(str(out/'assigned_backbone.cif'))
    library=load_library(args.rotamers)
    sidechains=build_sidechains(structure,density,library,progress=lambda s:print(s,flush=True))
    atom_geometry=atom_geometry_report(structure,library)
    structure.make_mmcif_document().write_file(str(out/'model.cif'))
    report={'version':__version__,'engine':'native-phasebuild-loops-and-sidechains',
            'parameters':vars(args),'input':density.source,'model':fingerprint(args.model),
            'profile':fingerprint(args.profile),'rotamers':fingerprint(args.rotamers),
            'sequence_input':fingerprint(args.sequence),'residue_profile':fingerprint(args.residue_profile),
            'sequence_report':fingerprint(args.sequence_report),'runner':fingerprint(__file__),
            'sources':sources,'source_changed_during_run':any(fingerprint(f['path'])!=f for f in sources),
            'loops':loop_report,'sequence_assignment':sequence,'sidechains':sidechains,
            'geometry':geometry_report(chains,density),'atom_geometry':atom_geometry,
            'total_seconds':time.perf_counter()-started,
            'scope':'No target reference or external building engine. Library derived from unrelated 1TQW. '
                    'Density fitting only; no refinement against raw observations or Rfree calculation.'}
    (out/'build_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'seconds':report['total_seconds'],'loops':loop_report,
                      'sidechain_atoms':sidechains['sidechain_atoms_built'],
                      'complete_sidechains':sidechains['complete_sidechains']},indent=2))


if __name__=='__main__':
    main()
