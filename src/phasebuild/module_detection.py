"""Which optional CCP4 modules this machine can really run.

OASIS4.3 always runs its own P+ phasing (module 1).  Parrot (module 2) and
Buccaneer/REFMAC (module 3) are CCP4 programs, so the window must not offer a
stage the machine cannot execute.  This module answers one question — "are the
executables there?" — and never raises: a broken or absent CCP4 installation is
a normal answer, not an error.

Search order: the run directory the form already carries (the packaged runtime),
then the usual CCP4 environment variables, then ``PATH``.
"""
import os
from pathlib import Path

# key, title, executables any one of which proves the module is installed
MODULES=(
    ('parrot','Parrot 密度修饰',('cparrot','parrot')),
    ('build','Buccaneer 建模 + REFMAC 精修',('cbuccaneer','buccaneer')),
)
# Buccaneer and REFMAC travel together: the builder writes the model, REFMAC
# refines it, and a machine with only one of them cannot finish the stage.
BUILD_PARTNER=('refmac5',)

MISSING_HINT='未检测到 CCP4：只做 P+ 定相，输出相位/FOM/MTZ。需要密度修饰或建模时，'
MISSING_TAIL='请安装 CCP4，或设置 CCP4/CBIN 环境变量、在“Buccaneer / REFMAC 运行目录”里填 bin 目录，再点“重新检测”。'


def _candidate_dirs(runtime=None,environ=None):
    """Directories that could hold the CCP4 executables, in search order."""
    env=os.environ if environ is None else environ
    dirs=[]
    def add(value):
        if not value:return
        path=Path(value).expanduser()
        if path.is_dir() and path not in dirs:dirs.append(path)
    if runtime:
        add(Path(runtime).expanduser()/'bin')
        add(runtime)
    for name in ('CBIN','CCP4_BIN','OASIS_CCP4_BIN'):
        add(env.get(name))
    for name in ('CCP4','CCP4_ROOT'):
        root=env.get(name)
        if root:
            add(Path(root)/'bin')
            add(root)
    # A double-clicked GUI inherits almost no environment, so also look where
    # the installers actually put CCP4 on this platform.
    for pattern in ('/Applications/ccp4-*/bin','/opt/ccp4*/bin','/usr/local/ccp4*/bin'):
        for path in sorted(Path('/').glob(pattern.lstrip('/'))):
            add(path)
    for entry in str(env.get('PATH') or '').split(os.pathsep):
        add(entry)
    return dirs


def _find(names,dirs):
    for folder in dirs:
        for name in names:
            path=folder/name
            if path.is_file() and os.access(path,os.X_OK):
                return str(path)
    return ''


def detect(runtime=None,environ=None):
    """Report module availability; every field is plain JSON-safe data."""
    dirs=_candidate_dirs(runtime,environ)
    found={key:_find(names,dirs) for key,_,names in MODULES}
    refmac=_find(BUILD_PARTNER,dirs)
    modules={}
    for key,title,names in MODULES:
        available=bool(found[key]) and (key!='build' or bool(refmac))
        missing=[]
        if not found[key]:missing.append('/'.join(names))
        if key=='build' and not refmac:missing.append('refmac5')
        modules[key]=dict(title=title,available=available,path=found[key],
                          reason='' if available else '未找到可执行文件：'+'、'.join(missing))
    if refmac:modules['build']['refmac']=refmac
    ready=[modules[key]['title'] for key,_,_ in MODULES if modules[key]['available']]
    if len(ready)==2:
        summary='已检测到 CCP4：Parrot 与 Buccaneer/REFMAC 默认启用（可在上方取消勾选）。'
    elif ready:
        summary='只检测到 CCP4 的一部分：默认启用 %s；%s 未找到，保持关闭。'%(
            '、'.join(ready),'、'.join(modules[key]['title'] for key,_,_ in MODULES if not modules[key]['available']))
    else:
        summary=MISSING_HINT+MISSING_TAIL
    return dict(available=bool(ready),modules=modules,summary=summary,
                ccp4_ready=all(modules[key]['available'] for key,_,_ in MODULES),
                searched=[str(path) for path in dirs])


def summary(config=None):
    """Convenience wrapper for the window: detect from a configuration dict."""
    config=config or {}
    return detect(config.get('runtime') or None)
