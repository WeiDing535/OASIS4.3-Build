"""Opt-in fixed-protocol P+ exploration. Existing single-run entry is unchanged."""
import argparse
import json
import os
from pathlib import Path
import sys

for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS']:
    os.environ[key] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, help='Pipeline config JSON, not a pipeline report')
    parser.add_argument('--gui-config', type=Path, help='Saved GUI job JSON')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--workers', type=int, choices=[1, 2], default=2)
    parser.add_argument('--protocols', nargs='+', help='Subset of the four predefined protocols')
    parser.add_argument('--check', action='store_true', help='Read-only preflight; creates no output')
    parser.add_argument('--child-config', type=Path, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if args.child_config:
        from phasebuild.sad_pipeline import run_pipeline
        report = run_pipeline(json.loads(args.child_config.read_text()))
        return 0 if report['status'] == 'completed' else 1
    if bool(args.config) == bool(args.gui_config) or not args.output:
        parser.error('Exactly one of --config/--gui-config and --output are required')
    from phasebuild.portfolio import plan, run
    config = json.loads((args.gui_config or args.config).read_text())
    if args.gui_config:
        config = dict(config['config'])
        for key in ['mode', 'explore', 'explore_workers', 'protocol']:
            config.pop(key, None)
    result = (plan(config, args.output, args.protocols) if args.check else
              run(config, args.output, workers=args.workers, names=args.protocols))
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0 if args.check or result['status'] == 'completed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
