import argparse
import json
import platform
import time
from pathlib import Path
import gemmi
import numpy as np
import scipy
from . import __version__
from .data import Density, audit, fingerprint


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def main():
    parser = argparse.ArgumentParser(description='PhaseBuild experimental backbone prototype')
    sub = parser.add_subparsers(dest='command', required=True)
    a = sub.add_parser('audit')
    a.add_argument('mtz')
    a.add_argument('--output')
    t = sub.add_parser('trace')
    t.add_argument('mtz')
    t.add_argument('--amplitude', default='FWT')
    t.add_argument('--phase', default='PHWT')
    t.add_argument('--output', required=True)
    t.add_argument('--rotations', type=int, default=512)
    t.add_argument('--refine-count', type=int, default=96)
    t.add_argument('--max-length', type=int, default=30)
    t.add_argument('--max-residues', type=int, default=65)
    t.add_argument('--beam-width', type=int, default=4)
    t.add_argument('--min-score', type=float, default=.9)
    t.add_argument('--seed', type=int, default=0)
    t.add_argument('--profile', help='Local density model trained on an unrelated structure')
    t.add_argument('--selection-cutoff', type=float, help='Enable global profile-based fragment selection')
    training = sub.add_parser('train-target')
    training.add_argument('--mtz', required=True)
    training.add_argument('--structure', required=True)
    training.add_argument('--amplitude', default='FC_BEST.F_phi.F')
    training.add_argument('--phase', default='FC_BEST.F_phi.phi')
    training.add_argument('--resolution', type=float, default=1.9)
    training.add_argument('--output', required=True)
    external = sub.add_parser('build', help='Explicitly labelled external building backend')
    external.add_argument('mtz')
    external.add_argument('--engine', choices=['buccaneer'], required=True)
    external.add_argument('--output', required=True)
    external.add_argument('--runtime', '--ccp4', dest='ccp4')
    external.add_argument('--cycles', type=int, default=3)
    external.add_argument('--timeout', type=int, default=600)
    selection = sub.add_parser('select', help='Re-score archived candidates without repeating the search')
    selection.add_argument('mtz')
    selection.add_argument('--profile', required=True)
    selection.add_argument('--candidates', required=True)
    selection.add_argument('--output', required=True)
    selection.add_argument('--cutoff', type=float, default=2.5)
    selection.add_argument('--max-residues', type=int, default=65)
    e = sub.add_parser('evaluate')
    e.add_argument('model')
    e.add_argument('--reference', required=True)
    e.add_argument('--cutoff', type=float, default=1.)
    e.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.command == 'audit':
        result = audit(args.mtz)
        if args.output:
            save(args.output, result)
        else:
            print(json.dumps(result, indent=2))
    elif args.command == 'select':
        from .profile import DensityProfile
        from .selection import select_fragments
        from .trace import write_backbone
        started = time.perf_counter()
        density = Density(args.mtz)
        density.profile = DensityProfile.load(args.profile)
        with np.load(args.candidates, allow_pickle=False) as data:
            metadata = json.loads(str(data['metadata']))
            if metadata['input']['sha256'] != density.source['sha256']:
                raise ValueError('Candidate archive was generated from a different MTZ')
            candidates = list(data['seeds'])
            start = 0
            for length in data['grown_lengths']:
                candidates.append(data['grown'][start:start+length])
                start += length
        chains, stats = select_fragments(density, candidates, args.max_residues, args.cutoff)
        out = Path(args.output)
        out.mkdir(parents=True, exist_ok=False)
        write_backbone(out / 'backbone.cif', chains, density)
        stats.update(engine='native-phasebuild-cached-selection', input=fingerprint(args.mtz),
                     profile=fingerprint(args.profile), candidates=fingerprint(args.candidates),
                     parameters=vars(args), total_seconds=time.perf_counter()-started)
        save(out / 'build_report.json', stats)
        print(json.dumps(stats, indent=2))
    elif args.command == 'build':
        from .backends import buccaneer
        print(json.dumps(buccaneer(args.mtz, args.output, args.ccp4, args.cycles, args.timeout), indent=2))
    elif args.command == 'train-target':
        from .profile import train
        if Path(args.output).exists():
            parser.error('Training output already exists')
        density = Density(args.mtz, args.amplitude, args.phase, resolution=args.resolution)
        profile = train(density, args.structure, args.mtz)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        profile.save(args.output)
        save(str(args.output) + '.json', profile.metadata)
        print(json.dumps(profile.metadata, indent=2))
    elif args.command == 'trace':
        from .trace import trace, write_backbone
        for name in ['rotations', 'refine_count', 'max_length', 'max_residues', 'beam_width']:
            if getattr(args, name) < 1:
                parser.error('--' + name.replace('_', '-') + ' must be positive')
        if args.max_length < 3 or args.max_residues < 3:
            parser.error('Length and residue budgets must be at least 3')
        if args.selection_cutoff is not None and not args.profile:
            parser.error('--selection-cutoff requires --profile')
        out = Path(args.output)
        # Never mix artifacts from different runs or overwrite existing results.
        out.mkdir(parents=True, exist_ok=False)
        started = time.perf_counter()
        density = Density(args.mtz, args.amplitude, args.phase)
        if args.profile:
            from .profile import DensityProfile
            density.profile = DensityProfile.load(args.profile)
        density.write_map(out / 'map_sigma.ccp4')
        source_files = [fingerprint(p) for p in sorted(Path(__file__).parent.glob('*.py'))]
        chains, stats = trace(density, args.rotations, args.refine_count, args.max_length,
                              args.beam_width, args.max_residues, args.min_score, args.seed,
                              progress=lambda msg: print(msg, flush=True),
                              candidate_output=out / 'candidates.npz', selection_cutoff=args.selection_cutoff)
        write_backbone(out / 'backbone.cif', chains, density)
        stats.update(version=__version__, input=fingerprint(args.mtz), parameters=vars(args),
                     engine='native-phasebuild',
                     profile=fingerprint(args.profile) if args.profile else None,
                     source_files=source_files,
                     source_changed_during_run=any(fingerprint(p['path']) != p for p in source_files),
                     map_shape=list(density.shape), map_mean_before_normalization=density.mean,
                     map_sigma_before_normalization=density.sigma,
                     total_seconds=time.perf_counter() - started,
                     environment={'python': platform.python_version(), 'gemmi': gemmi.__version__,
                                  'numpy': np.__version__, 'scipy': scipy.__version__})
        save(out / 'build_report.json', stats)
        print(json.dumps(stats, indent=2))
    else:
        from .evaluate import evaluate
        result = evaluate(args.model, args.reference, args.cutoff)
        save(args.output, result)
        print(json.dumps({k: v for k, v in result.items() if k != 'matches'}, indent=2))
