"""Autonomous post-tracing improvement; no evaluation reference argument."""
import argparse
import json
import sys
import time
import shutil
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
import numpy as np
from phasebuild.data import Density, fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.model import read_backbone
from phasebuild.completion import complete_internal_gaps
from phasebuild.validation import geometry_report
from phasebuild.trace import write_backbone
from phasebuild import __version__


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mtz')
    p.add_argument('--model',required=True)
    p.add_argument('--profile',required=True)
    p.add_argument('--candidates',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--max-gap',type=int,default=2)
    p.add_argument('--sequence')
    p.add_argument('--residue-profile')
    p.add_argument('--sequence-shuffles',type=int,default=499)
    p.add_argument('--seed',type=int,default=0)
    p.add_argument('--polish-terminals',action='store_true')
    p.add_argument('--map-amplitude',default='FWT')
    p.add_argument('--map-phase',default='PHWT')
    args = p.parse_args()
    if bool(args.sequence)!=bool(args.residue_profile):
        p.error('--sequence and --residue-profile must be provided together')
    out = Path(args.output)
    out.mkdir(parents=True,exist_ok=False)
    started = time.perf_counter()
    source = [fingerprint(f) for f in sorted((Path(__file__).parent/'src/phasebuild').glob('*.py'))]
    density = Density(args.mtz,args.map_amplitude,args.map_phase)
    density.profile = DensityProfile.load(args.profile)
    chains = read_backbone(args.model,density)
    before = geometry_report(chains,density)
    with np.load(args.candidates,allow_pickle=False) as z:
        if json.loads(str(z['metadata']))['input']['sha256'] != density.source['sha256']:
            raise ValueError('Candidate archive belongs to another MTZ')
        parents,start = [],0
        for n in z['grown_lengths']:
            parents.append(z['grown'][start:start+n].copy())
            start += n
    chains,completion = complete_internal_gaps(chains,parents,density,args.max_gap)
    polish_report=None
    if args.polish_terminals:
        from phasebuild.polish import fit_terminals
        chains,polish_report=fit_terminals(chains,density)
    write_backbone(out/'backbone.cif',chains,density)
    sequence_report = None
    if args.sequence:
        from phasebuild.residue_profile import ResidueProfile
        from phasebuild.sequencing import read_sequences, assign_sequence, annotate_structure
        records = read_sequences(args.sequence)
        target = ResidueProfile.load(args.residue_profile)
        sequence_report = assign_sequence(chains,density,target,records,shuffles=args.sequence_shuffles,seed=args.seed)
        annotate_structure(out/'backbone.cif',out/'model.cif',sequence_report,records)
    else:
        shutil.copyfile(out/'backbone.cif',out/'model.cif')
    report = {'version':__version__,'engine':'native-phasebuild-post-tracing','input':density.source,
              'model':fingerprint(args.model),'profile':fingerprint(args.profile),
              'candidates':fingerprint(args.candidates),'runner':fingerprint(__file__),
              'source_files':source,'source_changed_during_run':any(fingerprint(f['path'])!=f for f in source),
              'parameters':vars(args),'completion':completion,'geometry_before':before,
              'geometry_after':geometry_report(chains,density),'total_seconds':time.perf_counter()-started,
              'sequence_assignment':sequence_report,
              'terminal_polish':polish_report,
              'sequence_input':fingerprint(args.sequence) if args.sequence else None,
              'residue_profile':fingerprint(args.residue_profile) if args.residue_profile else None,
              'scope':'Autonomous map/geometry improvement. No evaluation reference is read.'}
    (out/'build_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
