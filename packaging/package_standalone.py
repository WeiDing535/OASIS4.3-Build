"""Assemble a fully self-contained OASIS4.3 folder from verified components.

The result runs without any separately installed CCP4 or Phenix: the bundled
runtime carries Phaser, Buccaneer, REFMAC and Parrot, and the bundled Python
carries NumPy, SciPy and Gemmi. Sources are never modified and the destination
must be a new path.
"""
import argparse
import hashlib
import json
import os
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
ENGINES=('phaser','cbuccaneer','refmac5','cparrot')
HASH_CHUNK=1<<20


def digest(path):
    handle=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(HASH_CHUNK), b''):
            handle.update(block)
    return handle.hexdigest()


def verify_manifest(root, name='runtime_manifest.json'):
    root=Path(root);manifest=root/name
    if not manifest.is_file():return {'manifest':None,'files_checked':0}
    rows=json.loads(manifest.read_text())['files']
    for row in rows:
        path=root/row['relative_path']
        if not path.is_file():raise ValueError('Missing runtime file: '+row['relative_path'])
        if digest(path)!=row['sha256']:raise ValueError('Runtime file differs: '+row['relative_path'])
    return {'manifest':digest(manifest),'files_checked':len(rows)}


def app_files(root=ROOT):
    """Reuse the focused source selection so the bundle ships only the product."""
    import importlib.util
    spec=importlib.util.spec_from_file_location('oasis_packager',ROOT/'package_oasis.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module.source_files(root)


def install_parrot(runtime, binary):
    """Add cparrot to the bundled runtime and record it in the manifest."""
    runtime=Path(runtime);target=runtime/'bin'/'cparrot'
    source=Path(binary) if binary else None
    if source is None:
        for candidate in (runtime/'bin'/'cparrot',Path('/Applications/ccp4-9/bin/cparrot')):
            if candidate.is_file():source=candidate;break
    if source is None:raise ValueError('No cparrot binary available; pass --parrot-binary')
    if not target.is_file() or digest(target)!=digest(source):
        shutil.copy2(source,target)
        os.chmod(target,os.stat(target).st_mode|stat.S_IXUSR|stat.S_IXGRP|stat.S_IXOTH)
    manifest_path=runtime/'runtime_manifest.json'
    manifest=json.loads(manifest_path.read_text())
    relative='bin/cparrot'
    entry=digest(target)
    # Self-referential bookkeeping is declared excluded from the inventory; keeping
    # a hash of the checksum file inside the manifest cannot stay consistent.
    excluded=set(manifest.get('inventory_excludes') or ())|{'runtime_manifest.json','runtime_manifest.sha256'}
    manifest['files']=[row for row in manifest['files']
                       if row['relative_path']!=relative and row['relative_path'] not in excluded]
    manifest['files'].append({'relative_path':relative,'sha256':entry})
    manifest['files'].sort(key=lambda row:row['relative_path'])
    engines=list(manifest.get('bundled_external_engines') or [])
    if 'Parrot' not in engines:engines.append('Parrot')
    manifest['bundled_external_engines']=engines
    manifest['copied_file_count']=int(manifest.get('copied_file_count',0))+1
    manifest['inventory_file_count']=int(manifest.get('inventory_file_count',0))+1
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    checksum=manifest_path.with_suffix('.sha256')
    if checksum.exists():
        checksum.write_text('%s  runtime_manifest.json\n'%digest(manifest_path))
    return {'binary':str(target),'sha256':entry,'manifest':digest(manifest_path)}


def write_launcher(output):
    launcher=output/'OASIS4.3.command'
    launcher.write_text('''#!/bin/bash
# OASIS4.3 self-contained launcher: no CCP4, Phenix or Python installation needed.
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
export OASIS4_ROOT="$HERE"
export PYTHONHOME="$HERE/python"
export PYTHONPATH="$HERE/app/src"
export PATH="$HERE/python/bin:$HERE/runtime/bin:/usr/bin:/bin:/usr/sbin:/sbin"
export CCP4="$HERE/runtime"
export CLIBD="$HERE/runtime/lib/data"
export CLIBD_MON="$HERE/runtime/lib/data/monomers/"
export SYMINFO="$HERE/runtime/lib/data/syminfo.lib"
export CINCL="$HERE/runtime/include"
export CCP4_SCR="${TMPDIR:-/tmp}/oasis43-scratch"
mkdir -p "$CCP4_SCR"
if [ "$#" -gt 0 ]; then
  exec "$HERE/python/bin/python3" -u "$HERE/app/$1" "${@:2}"
fi
exec "$HERE/python/bin/python3" -u "$HERE/app/run_gui.py" --project "$HERE"
''')
    os.chmod(launcher,0o755)
    return launcher


USAGE='''OASIS4.3 独立安装包

解压整个目录到可写位置，双击 OASIS4.3.command（macOS）即可启动界面。
命令行用法：./OASIS4.3.command run_sad_to_model.py --help

本包自带运行所需全部组件，不需要另外安装 CCP4、Phenix、Python、NumPy、SciPy 或 Gemmi：
  runtime/bin   Phaser、Buccaneer、REFMAC、Parrot 及所需动态库
  runtime/lib/data  对称信息、单体字典、参考直方图
  python/      Python 3.12 + NumPy/SciPy/Gemmi
  app/         OASIS4.3 程序与 src/

推荐的单波长反常散射（SAD）全流程：
  ./OASIS4.3.command run_sad_to_model.py --sad SAD.mtz --native NATIVE.mtz \\
     --sequence SEQ.fasta --heavy HEAVY.pdb --wavelength 0.97918 --elements Se \\
     --method pplus_bundled --dm-engine parrot --feedback-cycles 2 \\
     --dm-solvent 0.5 --building-recipe pipeline --output NEW_DIR

说明：
- 输出总是写入新目录，不覆盖已有结果。
- 每个案例的关系上限（--max-per-target）过大时程序会提示内存预算不足，按提示调小即可。
- 结果取舍只使用工作集指标（工作集 R 因子）；R-free 仅用于最终独立评估。
- 这是本机评估包：外部引擎 Buccaneer/REFMAC/Phaser/Parrot 的公开再分发条件尚未落实，
  不代表可以公开发布。详见 runtime/THIRD_PARTY_NOTICES.md 与 runtime/share/licenses。
'''


def verify_bundle(output):
    """Run the bundle with a minimal environment that cannot see host CCP4/Phenix."""
    output=Path(output).resolve()
    env={'PATH':'/usr/bin:/bin:/usr/sbin:/sbin','HOME':str(output/'.verify-home'),
         'LANG':'C','LC_ALL':'C'}
    report={}
    python=output/'python'/'bin'/'python3'
    if not python.is_file():raise ValueError('Bundled Python is missing')
    probe=subprocess.run([str(python),'-c',
        "import sys, numpy, scipy, gemmi; print(sys.version.split()[0], numpy.__version__, scipy.__version__, gemmi.__version__)"],
        env=env,capture_output=True,text=True,timeout=180)
    if probe.returncode:raise RuntimeError('Bundled Python cannot import its packages: '+probe.stderr)
    report['python']=probe.stdout.strip()
    runtime=output/'runtime'
    for name in ENGINES:
        binary=runtime/'bin'/name
        if not binary.is_file():raise ValueError('Bundled engine missing: '+name)
        if not os.access(binary,os.X_OK):raise ValueError('Bundled engine is not executable: '+name)
    report['engines']={name:digest(runtime/'bin'/name) for name in ENGINES}
    # Engine startup: PARROT with no input must still initialise its libraries.
    parrot=subprocess.run([str(runtime/'bin'/'cparrot')],env={**env,'CCP4':str(runtime),
        'CLIBD':str(runtime/'lib/data'),'CCP4_SCR':str(output/'.verify-home'/'scratch')},
        capture_output=True,text=True,timeout=180)
    banner=(parrot.stdout+parrot.stderr).lower()
    report['parrot_startup']='cparrot' in banner and 'ccp4' in banner
    if not report['parrot_startup']:raise RuntimeError('Bundled Parrot did not initialise')
    # End-to-end preflight of the pipeline with the bundled engines only.
    check=subprocess.run([str(python),'-u',str(output/'app'/'run_sad_to_model.py'),'--help'],
        env=env,capture_output=True,text=True,timeout=180)
    report['cli_import']=check.returncode==0 and 'dm-engine' in check.stdout
    if not report['cli_import']:raise RuntimeError('Bundled CLI failed to start: '+(check.stderr or check.stdout)[:800])
    # The Parrot adapter must resolve the bundled binary, never a host installation.
    locate=subprocess.run([str(python),'-c',
        "import sys;sys.path.insert(0,sys.argv[1]);from phasebuild import parrot_dm;"
        "print(parrot_dm.locate(None,sys.argv[2]))",
        str(output/'app'/'src'),str(output/'runtime')],env=env,capture_output=True,text=True,timeout=120)
    resolved=locate.stdout.strip()
    report['parrot_resolved']=resolved
    if Path(resolved).resolve()!= (output/'runtime'/'bin'/'cparrot').resolve():
        raise RuntimeError('Parrot resolved outside the bundle: '+resolved)
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--app-source',type=Path,default=ROOT)
    parser.add_argument('--python-source',type=Path,required=True)
    parser.add_argument('--runtime-source',type=Path,required=True)
    parser.add_argument('--parrot-binary',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    archives=parser.add_mutually_exclusive_group()
    archives.add_argument('--zip',type=Path)
    archives.add_argument('--tar',type=Path)
    args=parser.parse_args(argv)
    output=args.output.resolve()
    if output.exists():raise FileExistsError('Refusing to overwrite: '+str(output))
    report={'product':'OASIS4.3','platform':'%s-%s'%(platform.system().lower(),platform.machine().lower()),
            'self_contained':True,'requires_separate_ccp4_installation':False,
            'requires_separate_phenix_installation':False,'dependencies_installed':False,
            'distribution_status':'local-evaluation','public_redistribution_ready':False,
            'bundled_engines':list(ENGINES),'validation':{}}
    output.mkdir(parents=True)
    print('Verify sources',flush=True)
    report['validation']['python_source']=verify_manifest(args.python_source.resolve(),'python_manifest.json')
    report['validation']['runtime_source']=verify_manifest(args.runtime_source.resolve(),'runtime_manifest.json')
    for label,source in [('python',args.python_source),('runtime',args.runtime_source)]:
        print('Copy '+label,flush=True)
        shutil.copytree(Path(source).resolve(),output/label,symlinks=True)
    print('Add Parrot to the bundled runtime',flush=True)
    report['validation']['parrot']=install_parrot(output/'runtime',args.parrot_binary)
    print('Copy application sources',flush=True)
    app=output/'app';app.mkdir()
    source_root=args.app_source.resolve()
    files=app_files(source_root)
    for path in files:
        relative=Path(path).relative_to(source_root)
        target=app/relative
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(path,target)
    report['app_files']=len(files)
    write_launcher(output)
    (output/'使用说明.txt').write_text(USAGE)
    if (source_root/'LICENSE').is_file():shutil.copy2(source_root/'LICENSE',output/'LICENSE')
    print('Verify bundle',flush=True)
    report['validation']['bundle']=verify_bundle(output)
    report['status']='assembled-and-verified'
    (output/'bundle.json').write_text(json.dumps(report,indent=2,ensure_ascii=False)+'\n')
    shutil.rmtree(output/'.verify-home',ignore_errors=True)
    print(json.dumps(report,indent=1)[:1500])
    if args.zip:
        archive=Path(args.zip)
        if archive.exists():raise FileExistsError('Refusing to overwrite archive')
        with zipfile.ZipFile(archive,'x',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
            for path in sorted(output.rglob('*')):
                z.write(path,str(output.name)+'/'+str(path.relative_to(output)))
        _finish_archive(archive)
    if args.tar:
        archive=Path(args.tar)
        if archive.exists():raise FileExistsError('Refusing to overwrite archive')
        # Stream through the system tar: Python's tarfile keeps the whole gzip
        # stream open in one process and produced truncated archives on this host.
        completed=subprocess.run(['tar','-czf',str(archive),'-C',str(output.parent),output.name],
                                 capture_output=True,text=True)
        if completed.returncode:
            raise RuntimeError('tar failed: '+(completed.stderr or '')[:500])
        _finish_archive(archive)


def _finish_archive(archive):
    """Verify the archive really decompresses before publishing its checksum."""
    archive=Path(archive)
    with archive.open('rb') as handle:
        head=handle.read(2)
    if head!=b'\x1f\x8b':
        raise RuntimeError('Archive is not gzip: '+str(archive))
    if archive.name.endswith(('.tar.gz','.tgz')):
        test=subprocess.run(['gzip','-t',str(archive)],capture_output=True,text=True)
        if test.returncode:
            raise RuntimeError('Archive is truncated: '+(test.stderr or '').strip()[:300])
        listing=subprocess.run(['tar','-tzf',str(archive)],capture_output=True,text=True)
        if listing.returncode:raise RuntimeError('Archive listing failed')
        entries=len(listing.stdout.splitlines())
    else:
        with zipfile.ZipFile(archive) as z:
            if z.testzip() is not None:raise RuntimeError('ZIP integrity check failed')
            entries=len(z.namelist())
    Path(str(archive)+'.sha256').write_text(digest(archive)+'  '+archive.name+'\n')
    print('archive %s entries=%d bytes=%d'%(archive,entries,archive.stat().st_size),flush=True)


if __name__=='__main__':
    main()
