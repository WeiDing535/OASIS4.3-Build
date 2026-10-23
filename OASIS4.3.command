#!/bin/bash
cd -- "$(dirname -- "$0")" || exit 1
export PYTHONDONTWRITEBYTECODE=1
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
runtime="${PHASEBUILD_PYTHON:-}"
if [ -z "$runtime" ] && [ -x "$PWD/.conda-unified/bin/python" ]; then
  runtime="$PWD/.conda-unified/bin/python"
fi
if [ -z "$runtime" ] && [ -n "${CONDA_PREFIX:-}" ]; then
  runtime="$CONDA_PREFIX/bin/python"
fi
if [ -z "$runtime" ]; then
  runtime="$(command -v python3.12 || true)"
fi
if [ ! -x "$runtime" ]; then
  echo '找不到运行环境：'"$runtime"
  read -r -p '按回车退出…'
  exit 1
fi
"$runtime" -B run_gui.py --project "$PWD" "$@"
status=$?
if [ "$status" -ne 0 ]; then
  read -r -p '启动失败，请阅读上方提示。按回车退出…'
fi
exit "$status"
