"""Grow missing termini from sequence-registered anchors, then fit sidechains."""
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
from phasebuild.loops import sequence_bindings
from phasebuild.sequencing import read_sequences,annotate_structure
from phasebuild.extension import extend_anchors
from phasebuild.trace import write_backbone
from phasebuild.rebuilding import transfer_properties
from phasebuild.sidechains import build_sidechains
from phasebuild.validation import geometry_report
from phasebuild.atom_validation import atom_geometry_report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mtz');p.add_argument('--model',required=True);p.add_argument('--sequence-report',required=True)
    p.add_argument('--sequence',required=True);p.add_argument('--profile',required=True)
    p.add_argument('--residue-profile',required=True);p.add_argument('--rotamers',required=True)
    p.add_argument('--output',required=True);p.add_argument('--maximum-extension',type=int,default=12)
    p.add_argument('--beam-width',type=int,default=6)
    p.add_argument('--map-amplitude',default='FWT')
    p.add_argument('--map-phase',default='PHWT')
    args=p.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter()
    sources=[fingerprint(f) for f in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
    density=Density(args.mtz,args.map_amplitude,args.map_phase);density.profile=DensityProfile.load(args.profile)
    target=ResidueProfile.load(args.residue_profile);library=load_library(args.rotamers)
    records=read_sequences(args.sequence);previous=json.loads(Path(args.sequence_report).read_text())
    if previous['sequence_input']['sha256']!=fingerprint(args.sequence)['sha256']:raise ValueError('Sequence differs')
    structure=gemmi.read_structure(args.model)
    bindings=sequence_bindings(structure,previous['sequence_assignment'],records)
    chains=read_backbone(args.model,density)
    chains,bindings,growth=extend_anchors(chains,bindings,records,density,target,
                                        maximum=args.maximum_extension,beam_width=args.beam_width,
                                        progress=lambda message:print(message,flush=True))
    write_backbone(out/'backbone.cif',chains,density)
    decisions=[]
    for ci,(chain,binding) in enumerate(zip(chains,bindings)):
        decision={'fragment_index':ci,'length':len(chain),'accepted':binding is not None}
        if binding is not None:
            decision.update(sequence_record=binding['record'],start=binding['start']+1,end=binding['end'],
                            reason='inherited sequence anchor plus map-supported terminal growth')
        decisions.append(decision)
    assignment={'decisions':decisions,'assigned_residues':sum(d['length'] for d in decisions if d['accepted']),
                'scope':'Registration inherited from audited anchors; new residues not an independent sequence-placement test.'}
    annotate_structure(out/'backbone.cif',out/'assigned_backbone.cif',assignment,records)
    grown=gemmi.read_structure(str(out/'assigned_backbone.cif'))
    transfer=transfer_properties(structure,grown)
    sidechains=build_sidechains(grown,density,library,progress=lambda message:print(message,flush=True))
    atom_geometry=atom_geometry_report(grown,library)
    grown.make_mmcif_document().write_file(str(out/'model.cif'))
    report={'version':__version__,'engine':'native-sequence-anchored-extension','parameters':vars(args),
            'input':density.source,'model':fingerprint(args.model),'sequence_input':fingerprint(args.sequence),
            'sequence_report':fingerprint(args.sequence_report),'profile':fingerprint(args.profile),
            'residue_profile':fingerprint(args.residue_profile),'rotamers':fingerprint(args.rotamers),
            'sources':sources,'runner':fingerprint(__file__),
            'source_changed_during_run':any(fingerprint(s['path'])!=s for s in sources),
            'extension':growth,'sequence_assignment':assignment,'property_transfer':transfer,
            'sidechains':sidechains,'geometry':geometry_report(chains,density),'atom_geometry':atom_geometry,
            'seconds':time.perf_counter()-started,
            'scope':'No target reference, Buccaneer coordinates or Rfree used by growth. '
                    'Input maps may carry prior model/free-set information.'}
    (out/'build_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'added_residues':growth['added_residues'],'events':growth['events'],
                      'seconds':report['seconds'],'atom_geometry':atom_geometry},indent=2))


if __name__=='__main__':
    main()
