"""Assemble a fully self-contained OASIS4.3 folder from verified components.

The result runs without any separately installed CCP4 or Phenix: the bundled
runtime carries Phaser, Buccaneer and REFMAC, density modification is the
in-tree oasis-dm v2, and the bundled Python carries NumPy, SciPy and Gemmi. Sources are never modified and the destination
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
# Density modification is oasis-dm v2 and ships in-tree, so the runtime needs no
# external density-modification program at all.
ENGINES=('phaser','cbuccaneer','refmac5')
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
# Fail loudly on a missing bundled file instead of silently using a local CCP4.
export PHASEBUILD_STANDALONE=1
export PHASEBUILD_RUNTIME="$HERE/runtime"
mkdir -p "$CCP4_SCR"
if [ "$#" -gt 0 ]; then
  exec "$HERE/python/bin/python3" -u "$HERE/app/$1" "${@:2}"
fi
exec "$HERE/python/bin/python3" -u "$HERE/app/run_desktop.py" --project "$HERE/app"
''')
    os.chmod(launcher,0o755)
    return launcher


USAGE='''OASIS4.3 独立安装包

解压整个目录到可写位置，双击 OASIS4.3.command（macOS）即可启动界面。
命令行用法：./OASIS4.3.command run_sad_to_model.py --help

本包自带运行所需全部组件，不需要另外安装 CCP4、Phenix、Python、NumPy、SciPy 或 Gemmi：
  runtime/bin   Phaser、Buccaneer、REFMAC 及所需动态库
  runtime/lib/data  对称信息、单体字典、参考直方图
  python/      Python 3.12 + NumPy/SciPy/Gemmi
  app/         OASIS4.3 程序与 src/

推荐的单波长反常散射（SAD）全流程：
  ./OASIS4.3.command run_sad_to_model.py --sad SAD.mtz --native NATIVE.mtz \\
     --sequence SEQ.fasta --heavy HEAVY.pdb --wavelength 0.97918 --elements Se \\
     --method pplus_bundled --dm-engine oasis-dm --feedback-cycles 2 \\
     --dm-solvent 0.5 --building-recipe pipeline --output NEW_DIR

ASU 拷贝数默认自动判定（--copies auto）：先读 RCSB FASTA 头里的 Chains A, B, C，
没有该注释时取 Matthews 溶剂含量 20-65% 的最小拷贝数；也可用 --copies 3 显式指定。

说明：
- 输出总是写入新目录，不覆盖已有结果。
- 每个案例的关系上限（--max-per-target）过大时程序会提示内存预算不足，按提示调小即可。
- 结果取舍只使用工作集指标（工作集 R 因子）；R-free 仅用于最终独立评估。
- 密度修饰由随程序提供的 oasis-dm v2 完成，包内不需要任何外部密度修饰程序。
- 这是本机评估包：外部引擎 Buccaneer/REFMAC/Phaser 的公开再分发条件尚未落实，
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
    present=[name for name in ENGINES if (runtime/'bin'/name).is_file()]
    report['engines']={name:digest(runtime/'bin'/name) for name in present}
    # End-to-end preflight of the pipeline with the bundled engines only.
    check=subprocess.run([str(python),'-u',str(output/'app'/'run_sad_to_model.py'),'--help'],
        env=env,capture_output=True,text=True,timeout=180)
    report['cli_import']=check.returncode==0 and 'dm-engine' in check.stdout
    if not report['cli_import']:raise RuntimeError('Bundled CLI failed to start: '+(check.stderr or check.stdout)[:800])
    # The desktop window must import inside the bundle: Tkinter ships with the
    # bundled Python, so the GUI needs no extra installation.
    desktop=subprocess.run([str(python),'-c',
        "import sys;sys.path.insert(0,sys.argv[1]);import tkinter;from phasebuild import desktop;"
        "print(tkinter.TkVersion,desktop.App.__name__)",str(output/'app'/'src')],
        env=env,capture_output=True,text=True,timeout=120)
    if desktop.returncode:raise RuntimeError('Bundled desktop UI cannot start: '+desktop.stderr[:800])
    report['desktop']=desktop.stdout.strip()
    launcher_help=subprocess.run([str(python),'-u',str(output/'app'/'run_desktop.py'),'--help'],
        env=env,capture_output=True,text=True,timeout=120)
    report['desktop_cli']=launcher_help.returncode==0 and 'project' in launcher_help.stdout
    if not report['desktop_cli']:raise RuntimeError('Bundled desktop launcher failed: '+(launcher_help.stderr or launcher_help.stdout)[:800])
    # The window must find the pipeline scripts from the bundle root as well.
    project_check=subprocess.run([str(python),'-c',
        "import sys;sys.path.insert(0,sys.argv[1]);from phasebuild.desktop import resolve_project;"
        "p=resolve_project(sys.argv[2]);assert p and (p/'run_sad_to_model.py').is_file(),'cannot resolve project';print(p)",
        str(output/'app/src'),str(output)],cwd=str(output),env=env,capture_output=True,text=True,timeout=120)
    report['desktop_project']=project_check.returncode==0
    if not report['desktop_project']:raise RuntimeError('Bundled desktop cannot locate its scripts: '+(project_check.stderr or project_check.stdout)[:800])
    # Density modification must resolve to the in-tree oasis-dm, never a host copy.
    locate=subprocess.run([str(python),'-c',
        "import sys;sys.path.insert(0,sys.argv[1]);from phasebuild import oasis_dm;"
        "print(oasis_dm.locate())",
        str(output/'app'/'src')],env=env,capture_output=True,text=True,timeout=120)
    resolved=locate.stdout.strip()
    report['oasis_dm_resolved']=resolved
    if not resolved or not Path(resolved).is_file():
        raise RuntimeError('Bundled oasis-dm did not resolve: '+(resolved or locate.stderr[:400]))
    return report


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--app-source',type=Path,default=ROOT)
    parser.add_argument('--python-source',type=Path,required=True)
    parser.add_argument('--runtime-source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--keep-optional',action='store_true',
                        help='keep optional RESOLVE/plotting/validation payloads (adds ~1.4 GB)')
    parser.add_argument('--keep-all-monomers',action='store_true',
                        help='keep every ligand dictionary (adds ~660 MB)')
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
            'bundled_engines':list(ENGINES),
            'validation':{}}
    output.mkdir(parents=True)
    print('Verify sources',flush=True)
    report['validation']['python_source']=verify_manifest(args.python_source.resolve(),'python_manifest.json')
    report['validation']['runtime_source']=verify_manifest(args.runtime_source.resolve(),'runtime_manifest.json')
    for label,source in [('python',args.python_source),('runtime',args.runtime_source)]:
        print('Copy '+label,flush=True)
        shutil.copytree(Path(source).resolve(),output/label,symlinks=True)
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
    if not args.keep_optional:
        print('Prune optional components',flush=True)
        import importlib.util
        spec=importlib.util.spec_from_file_location('oasis_prune',Path(__file__).with_name('prune_bundle.py'))
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        removed=module.prune(output)
        if not args.keep_all_monomers:removed+=module.prune_monomers(output)
        removed+=module.prune_late(output)
        saved,linked=module.dedupe(output)
        report['deduplicated_links']=linked
        for name in ('runtime/runtime_manifest.json','python/python_manifest.json'):
            kept,dropped=module.refresh_manifest(output/Path(name).parent,Path(name).name)
            report['validation'].setdefault('manifests',{})[name]={'rows_kept':kept,'rows_dropped':dropped}
        report['pruned_bytes']=sum(item['bytes'] for item in removed)+saved
        report['pruned_items']=len(removed)
        print('Pruned %d items, %.0f MB'%(len(removed),report['pruned_bytes']/1e6),flush=True)
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
