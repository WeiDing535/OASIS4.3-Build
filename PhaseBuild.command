#!/bin/bash
cd -- "$(dirname -- "$0")" || exit 1
export PYTHONDONTWRITEBYTECODE=1
export PHASEBUILD_RUNTIME="${PHASEBUILD_RUNTIME:-/Applications/ccp4-9}"
runtime="${PHASEBUILD_PYTHON:-/Applications/ccp4-9/bin/ccp4-python}"
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
