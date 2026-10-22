"""Report and optionally hash-verify the portable PhaseBuild installation."""
import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'src'))


def verify_manifest(root, manifest):
    records = json.loads(manifest.read_text())['files']
    issues = []
    for record in records:
        relative = record['relative_path']
        path = (root / relative).resolve()
        try:
            path.relative_to(root.resolve())
        except ValueError:
            issues.append('Path escapes bundle: ' + relative)
            continue
        expected = record.get('sha256')
        if expected is None and 'copied' in record:
            expected = record['copied']['sha256']
        if not path.is_file():
            issues.append('Missing: ' + relative)
        elif hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            issues.append('Hash differs: ' + relative)
    return {'checked_files': len(records), 'issues': issues, 'passed': not issues,
            'manifest_sha256': hashlib.sha256(manifest.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime')
    parser.add_argument('--verify', action='store_true', help='Verify every packaged file hash')
    parser.add_argument('--output')
    args = parser.parse_args()
    from phasebuild.runtime import resolve_runtime, runtime_info
    import gemmi
    import numpy
    import scipy
    runtime = resolve_runtime(args.runtime, required=('cbuccaneer', 'refmac5'))
    result = {'status': 'success', 'platform': platform.platform(), 'architecture': platform.machine(),
              'python': sys.executable, 'python_prefix': sys.prefix, 'runtime': runtime_info(runtime),
              'packages': {name: {'version': module.__version__, 'path': module.__file__}
                           for name, module in [('gemmi', gemmi), ('numpy', numpy), ('scipy', scipy)]},
              'scope': 'Local evaluation build; external Buccaneer and REFMAC contributions retained',
              'public_redistribution_ready': False}
    if args.verify:
        bundle = ROOT.parent
        result['verification'] = {}
        for name, folder, manifest_name in [('engines', runtime, 'runtime_manifest.json'),
                                             ('python', bundle / 'python', 'python_manifest.json'),
                                             ('app', ROOT, 'app_manifest.json')]:
            manifest = folder / manifest_name
            if not manifest.is_file():
                result['verification'][name] = {'passed': False, 'issues': ['Manifest missing: ' + str(manifest)]}
            else:
                result['verification'][name] = verify_manifest(folder, manifest)
        if not all(row['passed'] for row in result['verification'].values()):
            result['status'] = 'failed'
    text = json.dumps(result, indent=2, ensure_ascii=False) + '\n'
    if args.output:
        Path(args.output).write_text(text)
    print(text, end='')
    if result['status'] != 'success':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
