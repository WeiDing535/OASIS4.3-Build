#!/bin/bash
# OASIS4.3 GUI launcher (public release, no local hardcoded paths)
# This version does not bundle a CCP4 runtime; the external engine path is entered in the GUI advanced parameters.
# The Python interpreter is probed in the following order; the first one that is "executable and has numpy/gemmi/scipy installed" is chosen:
#   1) Environment variable PHASEBUILD_PYTHON
#   2) Bundled .conda-unified/bin/python (if included at packaging time)
#   3) Bundled .local_python file (created by local users, contains the interpreter absolute path; gitignored)
#   4) Currently active conda environment CONDA_PREFIX
#   5) python3.12 in PATH
PROJECT_DIR="$(cd -- "$(dirname -- "$0")" && pwd)"
cd -- "$PROJECT_DIR" || exit 1

export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1

# Runtime environment: CCP4 is specified by the user in the GUI; the launcher does not bind to a bundled path.

local_py=""
if [ -f "$PROJECT_DIR/.local_python" ]; then
  local_py="$(head -n 1 "$PROJECT_DIR/.local_python" | tr -d '\r')"
fi
path_py="$(command -v python3.12 || true)"
conda_py=""
if [ -n "${CONDA_PREFIX:-}" ]; then
  conda_py="$CONDA_PREFIX/bin/python"
fi

# Probing: the interpreter must be executable and able to import all scientific computing dependencies.
usable() {
  local py="$1"
  [ -n "$py" ] && [ -x "$py" ] && "$py" -c 'import numpy, gemmi, scipy' >/dev/null 2>&1
}

runtime=""
for candidate in \
  "${PHASEBUILD_PYTHON:-}" \
  "$PROJECT_DIR/.conda-unified/bin/python" \
  "$local_py" \
  "$conda_py" \
  "$path_py"; do
  if usable "$candidate"; then
    runtime="$candidate"
    break
  fi
  if [ -n "$candidate" ] && [ -x "$candidate" ]; then
    echo "  Skipping $candidate (runnable but missing numpy/gemmi/scipy dependencies)"
  fi
done

if [ -z "$runtime" ]; then
  echo 'No Python 3.12 environment with complete dependencies found. Choose one of the following:'
  echo '  1) First conda activate oasis43-unified (or a 3.12 environment with numpy/gemmi/scipy), then launch;'
  echo '  2) Or set the environment variable PHASEBUILD_PYTHON to the interpreter of that environment;'
  echo '  3) Or create a .local_python file in the package root directory and write the interpreter absolute path into it.'
  echo '  Dependencies: numpy 1.26.x, gemmi 0.7.5, scipy 1.13.x (see environment.yml).'
  read -r -p 'Press Enter to exit...'
  exit 1
fi
export PHASEBUILD_PYTHON="$runtime"

echo "OASIS4.3 launching..."
echo "  Python : $runtime"
echo "  Engines: user-supplied CCP4 (phaser / cparrot / buccaneer / refmac5)"
"$runtime" -B run_desktop.py --project "$PROJECT_DIR" "$@"
status=$?
if [ "$status" -ne 0 ]; then
  echo 'Launch failed; please read the messages above.'
  read -r -p 'Press Enter to exit...'
fi
exit "$status"
