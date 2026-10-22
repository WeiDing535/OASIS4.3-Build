"""Assemble application sources and relative launchers; no public publication."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = sorted(path.name for path in ROOT.glob('*.py'))
MODELS = ['1tqw_1p9_symmetric.npz', 'residue_1tqw_1p9.npz', 'rotamers_1tqw_v08.json']
LAUNCHER = '''#!/bin/sh
set -eu
base=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
unset PYTHONHOME PYTHONPATH CCP4 CLIBD CLIBD_MON CBIN CEXAM CINCL CLIB
unset DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH DYLD_FRAMEWORK_PATH DYLD_INSERT_LIBRARIES
export PHASEBUILD_STANDALONE=1 PHASEBUILD_RUNTIME="$base/runtime"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$base/app/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mode=${1:---help}
if [ "$#" -gt 0 ]; then shift; fi
case "$mode" in
  doctor) script=run_runtime_doctor.py ;;
  hybrid-hl) script=run_phasebuild_hybrid_hl.py ;;
  build-hl) script=run_build_from_phases.py ;;
  hybrid) script=run_phasebuild_hybrid.py ;;
  native) script=run_sad_development.py ;;
  refine) script=run_refine.py ;;
  rfactor) script=run_rfactor.py ;;
  phase-refine) script=run_phase_refinement.py ;;
  gui) exec "$base/python/bin/python3" "$base/app/run_gui.py" --project "$base/app" "$@" ;;
  audit|evaluate) exec "$base/python/bin/python3" "$base/app/run_phasebuild.py" "$mode" "$@" ;;
  --help|-h|help)
    cat <<'HELP'
PhaseBuild local portable preview (macOS x86_64)
  ./phasebuild gui                  Open the local GUI
  ./phasebuild doctor [--verify]    Inspect/verify bundled components
  ./phasebuild hybrid-hl --help     PhaseBuild seed + Buccaneer + original HL refinement
  ./phasebuild build-hl --help      Experimental phases, no initial model + original HL refinement
  ./phasebuild native --help        Experimental native PhaseBuild workflow
  ./phasebuild refine --help        REFMAC refinement
  ./phasebuild rfactor --help       Independent fixed-coordinate R factors
  ./phasebuild evaluate --help      Separate post-hoc structure evaluation
External Buccaneer/REFMAC components are explicitly included. This preview is
for local evaluation, with public redistribution conditions still unresolved.
HELP
    exit 0 ;;
  *) echo "Unknown command: $mode" >&2; exit 2 ;;
esac
exec "$base/python/bin/python3" "$base/app/$script" "$@"
'''


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default=str(ROOT / 'dist/local_macos_x86_64'))
    args = parser.parse_args(argv)
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    app = out / 'app'
    app.mkdir(exist_ok=False)
    shutil.copytree(ROOT / 'src', app / 'src', ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for name in SCRIPTS:
        shutil.copy2(ROOT / name, app / name)
    (app / 'models').mkdir()
    for name in MODELS:
        shutil.copy2(ROOT / 'models' / name, app / 'models' / name)
    # Preserve the existing first-round native command's bank discovery path;
    # only these two unrelated trained profiles are shipped, not research cases.
    bank = app / 'runs/sad_four_optimization/models'
    bank.mkdir(parents=True)
    for name in ['backbone_2p8.npz', 'residue_2p8.npz']:
        shutil.copy2(ROOT / 'runs/sad_four_optimization/models' / name, bank / name)
    # Ship source documentation and license alongside every application copy.
    for name in ['LICENSE', 'AUTHORS.md', 'README.md', 'README_EN.md', 'GUI_GUIDE.md',
                 'SAD_RESEARCH_GUIDE.md', 'WEAK_PHASE_DEVELOPMENT.md', 'pyproject.toml']:
        shutil.copy2(ROOT / name, app / name)
    shutil.copytree(ROOT / 'packaging', app / 'packaging',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    files = [{'relative_path': str(path.relative_to(app)), 'sha256': digest(path),
              'bytes': path.stat().st_size} for path in sorted(app.rglob('*')) if path.is_file()]
    (app / 'app_manifest.json').write_text(json.dumps({'distribution_status': 'local-evaluation',
        'scientific_protocol': 'Opt-in integrated: optional initial seed, external Buccaneer3/2/2/2/2 and fixed originalHL REFMAC5x10',
        'files': files}, indent=2) + '\n')
    (out / 'phasebuild').write_text(LAUNCHER)
    (out / 'phasebuild').chmod(0o755)
    (out / 'PhaseBuild.command').write_text('#!/bin/sh\nset -eu\nbase=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)\nexec "$base/phasebuild" gui "$@"\n')
    (out / 'PhaseBuild.command').chmod(0o755)
    print(json.dumps({'bundle': str(out), 'app_files': len(files),
                      'ready_for_local_validation': (out / 'runtime/runtime_manifest.json').is_file()
                        and (out / 'python/python_manifest.json').is_file()}, indent=2))


if __name__ == '__main__':
    main()
