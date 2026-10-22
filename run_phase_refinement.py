"""Optional phase-aware REFMAC branch; original experimental priors required."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild.phase_refinement import merge_experimental_phases, refmac_with_experimental_prior


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native', required=True)
    parser.add_argument('--experimental-phases', required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--cycles', type=int, default=10)
    parser.add_argument('--representation', choices=['hl', 'phifom'], default='hl')
    parser.add_argument('--amplitude-only-control', action='store_true')
    parser.add_argument('--runtime', '--ccp4', dest='runtime')
    a = parser.parse_args()
    folder = Path(a.output)
    folder.mkdir(parents=True, exist_ok=False)
    merged = folder / 'native_with_original_prior.mtz'
    merge_experimental_phases(a.native, a.experimental_phases, merged, a.representation)
    result = refmac_with_experimental_prior(merged, a.model, folder / 'refinement',
        use_prior=not a.amplitude_only_control, cycles=a.cycles, representation=a.representation, ccp4=a.runtime)
    print(json.dumps(result['statistics']['final']), flush=True)


if __name__ == '__main__':
    main()
