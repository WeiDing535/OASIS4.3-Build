"""Opt-in removal of extra protein fragments; refuses ambiguous copy/sequence cases."""
import argparse
import json
from pathlib import Path
import sys
import gemmi

sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild.data import fingerprint
from phasebuild.model_cleanup import cleanup_single_copy


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', required=True)
    parser.add_argument('--sequence', required=True)
    parser.add_argument('--asu-copies', required=True, type=int)
    parser.add_argument('--output', required=True, help='New directory; existing directories are refused')
    args = parser.parse_args(argv)
    text = Path(args.sequence).read_text()
    if sum(line.startswith('>') for line in text.splitlines()) > 1:
        raise ValueError('Exactly one sequence record is supported')
    sequence = ''.join(line.strip() for line in text.splitlines() if not line.startswith('>'))
    cleaned, report = cleanup_single_copy(gemmi.read_structure(args.model), sequence,
                                        asu_copies=args.asu_copies)
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    output = out / 'model.cif'
    cleaned.make_mmcif_document().write_file(str(output))
    report.update(input_model=fingerprint(args.model), input_sequence=fingerprint(args.sequence),
                  output_model=fingerprint(output))
    (out / 'cleanup_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'output': str(output),
                      'removed_protein_residues': report['removed_protein_residues']}))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
