"""Which optional CCP4 modules this machine can really run.

OASIS4.3 always runs P+ phasing together with oasis-dm density modification
(module 1): oasis-dm ships in-tree, so that half needs no detection at all.
Buccaneer/REFMAC (module 2) are CCP4 programs, so the window must not offer a
stage the machine cannot execute.  This module answers one question — "are the
executables there?" — and never raises: a broken or absent CCP4 installation is
a normal answer, not an error.

Search order: the location the user typed (a CCP4 root, its ``bin``, or a
``ccp4.setup``/``bashrc`` file to read ``CCP4``/``CBIN`` out of), then the run
directory the form already carries (the packaged runtime), then the usual CCP4
environment variables and install paths, then ``PATH``.
"""
import os
import re
from pathlib import Path

# key, title, executables any one of which proves the module is installed
MODULES=(
    ('build','Buccaneer 建模 + REFMAC 精修',('cbuccaneer','buccaneer')),
)
# Density modification is not detected: oasis-dm v2 is part of the program and is
# always available.  It stays in the report so the window can show one list.
IN_TREE=(('dm','oasis-dm v2 密度修饰'),)
# Buccaneer and REFMAC travel together: the builder writes the model, REFMAC
# refines it, and a machine with only one of them cannot finish the stage.
BUILD_PARTNER=('refmac5',)
# Where the installers put CCP4 on this platform.  A double-clicked GUI inherits
# almost no environment, so these are searched explicitly.
INSTALL_PATTERNS=('/Applications/ccp4-*/bin','/opt/ccp4*/bin','/usr/local/ccp4*/bin')

MISSING_HINT='未检测到 CCP4：只做模块 1（P+ 定相 + oasis-dm 密度修饰），输出相位/FOM/MTZ。'
MISSING_TAIL=('需要 Buccaneer 建模与 REFMAC 精修时请安装 CCP4；已装但不在默认位置时，'
              '在“CCP4 位置”里填 CCP4 根目录（含 bin/）或 ccp4.setup / bashrc 的路径，再点“重新检测”。')

# CCP4=/opt/ccp4-8, export CBIN=/opt/ccp4-8/bin, CCP4_ROOT='/Applications/ccp4-9'
_ASSIGNMENT=re.compile(r'''(?:^|[;\s])(?:export\s+)?(CCP4|CCP4_ROOT|CBIN|CCP4_BIN)\s*=\s*("[^"]*"|'[^']*'|[^\s;]+)''')


def _script_dirs(path):
    """CCP4 locations named in a setup script the user pointed us at."""
    try:
        text=path.read_text(errors='ignore')
    except OSError:
        return []
    home=str(Path.home())
    roots=[]
    for match in _ASSIGNMENT.finditer(text):
        value=match.group(2).strip('"\'').replace('${HOME}',home).replace('$HOME',home)
        if not value.startswith('/'):
            continue
        candidate=Path(value).expanduser()
        # CBIN points at <root>/bin; CCP4/CCP4_ROOT point at <root>.
        roots.extend([candidate.parent,candidate] if candidate.name=='bin' else [candidate,candidate/'bin'])
    return roots


def location_dirs(value):
    """Candidate directories for whatever the user typed, plus a short note."""
    if not value:
        return [],[]
    path=Path(value).expanduser()
    if path.is_file():
        roots=_script_dirs(path)
        if roots:
            return roots,['已从 %s 解析出 CCP4 位置'%path.name]
        return [],[path.name+' 里没有找到 CCP4 / CBIN 赋值']
    return ([path.parent,path] if path.name=='bin' else [path,path/'bin']),[]


def _candidate_dirs(runtime=None,environ=None,location=None):
    """Directories that could hold the CCP4 executables, in search order."""
    env=os.environ if environ is None else environ
    dirs=[]
    def add(value):
        if not value:return
        path=Path(value).expanduser()
        if path.is_dir() and path not in dirs:dirs.append(path)
    for path in location_dirs(location)[0]:
        add(path)
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
    for pattern in INSTALL_PATTERNS:
        for path in sorted(Path('/').glob(pattern.lstrip('/'))):
            add(path)
    for entry in str(env.get('PATH') or '').split(os.pathsep):
        add(entry)
    return dirs


def _install_root(path):
    """The installation root that owns an executable (strip a trailing bin/)."""
    parent=Path(path).parent
    return str(parent.parent if parent.name=='bin' else parent)


def _find(names,dirs):
    for folder in dirs:
        for name in names:
            path=folder/name
            if path.is_file() and os.access(path,os.X_OK):
                return str(path)
    return ''


def detect(runtime=None,environ=None,location=None):
    """Report module availability; every field is plain JSON-safe data."""
    typed,notes=location_dirs(location)
    dirs=_candidate_dirs(runtime,environ,location)
    found={key:_find(names,dirs) for key,_,names in MODULES}
    refmac=_find(BUILD_PARTNER,dirs)
    modules={key:dict(title=title,available=True,path='内置（随程序提供）',
                      reason='') for key,title in IN_TREE}
    for key,title,names in MODULES:
        available=bool(found[key]) and (key!='build' or bool(refmac))
        missing=[]
        if not found[key]:missing.append('/'.join(names))
        if key=='build' and not refmac:missing.append('refmac5')
        modules[key]=dict(title=title,available=available,path=found[key],
                          reason='' if available else '未找到可执行文件：'+'、'.join(missing))
    if refmac:modules['build']['refmac']=refmac
    ready=[modules[key]['title'] for key,_,_ in MODULES if modules[key]['available']]
    if ready:
        summary='已检测到 CCP4：Buccaneer/REFMAC 建模默认启用（可在上方取消勾选）。'
    else:
        summary=MISSING_HINT+MISSING_TAIL
    # Which installation actually answered.  ``found_at`` is where the
    # executables are; ``resolved`` is set only when that is the place the user
    # pointed at, so the window can tell "found it there" from "found it
    # somewhere else" and never send a run to a directory with no programs.
    found_at=''
    for path in found.values():
        if path:
            found_at=_install_root(path)
            break
    resolved=''
    if found_at:
        typed_roots=[str(item) for item in typed]
        if not typed or found_at in typed_roots:
            resolved=found_at
    return dict(available=bool(ready),modules=modules,summary=summary,
                ccp4_ready=all(modules[key]['available'] for key,_,_ in MODULES),
                location=dict(value=str(location or ''),resolved=resolved,found_at=found_at,
                              note='；'.join(notes)),
                searched=[str(path) for path in dirs])


def summary(config=None):
    """Convenience wrapper for the window: detect from a configuration dict."""
    config=config or {}
    runtime=config.get('runtime') or None
    location=config.get('ccp4_location') or None
    # A record may carry the setup script itself as its runtime; parse it.
    if not location and runtime and Path(runtime).expanduser().is_file():
        location=runtime
    return detect(runtime,location=location)
