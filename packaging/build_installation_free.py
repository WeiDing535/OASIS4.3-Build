"""Assemble a focused local bundle with zero separate suite installations.

Copies existing authorized runtime components; installs nothing. Never copies
Phenix. The included engine manifest retains its redistribution restrictions.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib

ROOT=Path(__file__).resolve().parents[1]
LAUNCHER='''#!/bin/sh
set -eu
base=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
unset PYTHONHOME PYTHONPATH CCP4 CLIBD CLIBD_MON CBIN CEXAM CINCL CLIB
unset LD_LIBRARY_PATH LD_PRELOAD LD_AUDIT
unset DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH DYLD_FRAMEWORK_PATH DYLD_INSERT_LIBRARIES
unset PHENIX PHENIX_ROOT PHENIX_VERSION LIBTBX_BUILD
export OASIS_INSTALLATION_FREE=1 PHASEBUILD_STANDALONE=1 PHASEBUILD_RUNTIME="$base/runtime"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1 PYTHONPATH="$base/app/src"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
mode=${1:-gui}
if [ "$#" -gt 0 ]; then shift; fi
case "$mode" in
 gui) exec "$base/python/bin/python3" -B "$base/app/run_desktop.py" --project "$base/app" "$@" ;;
 sad-to-model) exec "$base/python/bin/python3" -B "$base/app/run_sad_to_model.py" "$@" ;;
 portfolio) exec "$base/python/bin/python3" -B "$base/app/run_sad_portfolio.py" "$@" ;;
 doctor) exec "$base/python/bin/python3" -B "$base/app/portable_doctor.py" "$base" ;;
 *) echo 'Usage: ./oasis [gui|sad-to-model|portfolio|doctor]' >&2; exit 2 ;;
esac
'''
DOCTOR='''import sys,json,os
from pathlib import Path
import numpy,scipy,gemmi
from phasebuild.runtime import resolve_runtime
base=Path(sys.argv[1]).resolve()
paths={'python':sys.executable,'numpy':numpy.__file__,'scipy':scipy.__file__,'gemmi':gemmi.__file__}
for name,path in paths.items():Path(path).resolve().relative_to(base/'python')
assert sys.version_info[:2]==(3,12)
assert resolve_runtime(required=('cbuccaneer','refmac5'))==base/'runtime'
assert os.environ.get('OASIS_INSTALLATION_FREE')=='1'
from phasebuild.oasis_workflow import defaults
cfg=defaults(base/'app');assert cfg['method'] in ('pplus_internal','pplus_bundled') and not cfg['phenix_bin']
doc=json.loads((base/'runtime/runtime_manifest.json').read_text())
names={Path(e['relative_path']).name for e in doc.get('dependency_graph',[])}
labels={'phaser':'Phaser','cparrot':'Parrot','cbuccaneer':'Buccaneer','refmac5':'REFMAC'}
engines=[labels[n] for n in labels if n in names]
print(json.dumps({'status':'passed','paths':paths,'python_version':sys.version,'default_method':cfg['method'],'requires_separate_Phenix':False,'requires_separate_CCP4':False,'bundled_external_engines':engines,'public_redistribution_ready':False},indent=2))
'''


ENGINE_LABELS={'phaser':'Phaser','cparrot':'Parrot','cbuccaneer':'Buccaneer','refmac5':'REFMAC'}


def bundled_engines(root):
    doc=json.loads((Path(root)/'runtime'/'runtime_manifest.json').read_text())
    names={Path(entry['relative_path']).name for entry in doc.get('dependency_graph',[])}
    return [name for name in ENGINE_LABELS if name in names]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def manifest(root, *, exclude=()):
    root=Path(root).resolve();rows=[]
    for p in sorted(root.rglob('*')):
        if str(p.relative_to(root)) in exclude:continue
        if p.is_symlink():
            try:p.resolve(strict=True).relative_to(root)
            except (ValueError,OSError):raise ValueError('Broken/escaping copied symlink: '+str(p))
        if p.is_file():rows.append({'relative_path':str(p.relative_to(root)),'sha256':sha(p)})
    return rows


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--python-source',type=Path,required=True);p.add_argument('--runtime-source',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(argv);out=a.output.resolve()
    if out.exists():raise FileExistsError('Refusing existing output: '+str(out))
    for src in [a.python_source,a.runtime_source]:
        if src.resolve()==out or src.resolve() in out.parents:raise ValueError('Output must not be inside a source runtime')
    probe=subprocess.run([str(a.python_source/'bin/python3'),'-B','-c','import sys,numpy,scipy,gemmi;assert sys.version_info[:2]==(3,12)'],capture_output=True,text=True)
    if probe.returncode:raise ValueError('Source Python must already contain Python3.12/numpy/scipy/gemmi: '+probe.stderr)
    spec=importlib.util.spec_from_file_location('focused_packager',ROOT/'package_oasis.py');pack=importlib.util.module_from_spec(spec);spec.loader.exec_module(pack)
    engine_manifest=json.loads((a.runtime_source/'runtime_manifest.json').read_text())
    for row in engine_manifest['files']:
        if sha(a.runtime_source/row['relative_path'])!=row['sha256']:raise ValueError('Damaged source engine runtime')
    out.mkdir(parents=True,exist_ok=False);app=out/'app';app.mkdir()
    for source in pack.source_files():
        dst=app/source.relative_to(ROOT);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,dst)
    shutil.copy2(ROOT/'pyproject.toml',app/'pyproject.toml')
    (app/'portable_doctor.py').write_text(DOCTOR)
    print('Copy scientific Python 3.12 (existing environment; no installation)',flush=True)
    shutil.copytree(a.python_source,out/'python',symlinks=True,ignore=shutil.ignore_patterns('__pycache__','*.pyc','.DS_Store','conda-meta'))
    print('Copy existing authorized local engine runtime',flush=True)
    shutil.copytree(a.runtime_source,out/'runtime',symlinks=True)
    (out/'oasis').write_text(LAUNCHER);(out/'oasis').chmod(0o755)
    (out/'OASIS4.3.command').write_text('#!/bin/sh\nset -eu\nbase=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)\nexec "$base/oasis" gui "$@"\n');(out/'OASIS4.3.command').chmod(0o755)
    print('Check copied links, hash files and verify clean-environment doctor',flush=True)
    pyfiles=manifest(out/'python',exclude={'python_manifest.json'});appfiles=manifest(app,exclude={'app_manifest.json'});manifest(out/'runtime')
    (out/'python/python_manifest.json').write_text(json.dumps({'files':pyfiles,'version':'3.12','source':'Existing local scientific environment; not downloaded'},indent=2)+'\n')
    (app/'app_manifest.json').write_text(json.dumps({'files':appfiles},indent=2)+'\n')
    env={'PATH':'/usr/bin:/bin:/usr/sbin:/sbin','HOME':str(out.parent/(out.name+'-clean-home')),'LANG':'C'}
    run=subprocess.run([str(out/'oasis'),'doctor'],env=env,capture_output=True,text=True,timeout=120)
    if run.returncode:raise RuntimeError('Copied runtime doctor failed: '+run.stdout+run.stderr)
    # The bundle is meant to be moved and downloaded; absolute build paths in the
    # doctor report would be wrong on the receiving machine.
    doctor=json.loads(run.stdout)
    doctor['paths']={k:str(Path(v).relative_to(out)) if str(v).startswith(str(out)) else v
                     for k,v in doctor.get('paths',{}).items()}
    (out/'DOCTOR.json').write_text(json.dumps(doctor,indent=2)+'\n')
    engines=bundled_engines(out);parrot='cparrot' in engines
    scope=('Bundled Phaser SAD -> P+ -> Parrot density modification -> Buccaneer/REFMAC. '
        'The RESOLVE/CCTBX density runtime is not part of this assembly.') if parrot else \
        'Internal approximate likelihood and DM. External-branch SEC17 Rfree27.42% is not claimed for this route.'
    report={'version':tomllib.loads((ROOT/'pyproject.toml').read_text())['tool']['oasis43']['version'],'platform':'macOS Intel x86_64','status':'assembled',
        'default_method':'pplus_bundled' if {'phaser','cparrot'}<=set(engines) else 'pplus_internal',
        'requires_separate_Phenix_installation':False,'requires_separate_CCP4_installation':False,'Phenix_bundled':False,
        'bundled_external_engines':[ENGINE_LABELS[n] for n in ENGINE_LABELS if n in engines],
        'bundled_engine_ids':[n for n in ENGINE_LABELS if n in engines],
        'scope':'Local evaluation with existing authorized engine binaries. Public redistribution conditions unresolved.',
        'public_redistribution_ready':False,'dependencies_installed':False,'python_files':len(pyfiles),'app_files':len(appfiles),
        'engine_manifest_SHA256':sha(out/'runtime/runtime_manifest.json'),'doctor':doctor,
        'scientific_scope':scope}
    (out/'bundle.json').write_text(json.dumps(report,indent=2)+'\n')
    engine_text='、'.join(ENGINE_LABELS[n] for n in ENGINE_LABELS if n in engines)
    dm_text='Parrot 密度修饰' if parrot else '内置 RESOLVE/MLHL 密度修饰'
    (out/'使用说明.txt').write_text('OASIS4.3 本机免安装套件评估包（macOS Intel）\n\n双击 OASIS4.3.command；保持整个目录完整。无需另装 Python、Phenix 或整套 CCP4。\n随包引擎：'+engine_text+'；密度修饰：'+dm_text+'。内部 SAD/P+ 使用随包 NumPy/SciPy/Gemmi。不是全部算法由 OASIS 原创。\nPhenix 对照模式在此包中禁止调用。数据文件由用户提供；运行写入新目录。\n这是由本机已有授权组件构建的本机评估包；公开再分发条件尚未落实，不是公开发行包。\n')
    print(json.dumps({'bundle':str(out),'status':report['status'],'doctor':'passed'},indent=2))

if __name__=='__main__':main()
