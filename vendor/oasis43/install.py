#!/usr/bin/env python3
"""One-command local installer: python3.12 install.py [--no-launch]."""
from pathlib import Path
import argparse
import json
import os
import subprocess
import sys
import time
import tomllib


def probe_python(executable):
    """Check the standard library without trusting a copied launcher's path."""
    code='import encodings,json,sys; print(json.dumps(dict(version=list(sys.version_info[:2]),prefix=sys.prefix,base_prefix=sys.base_prefix)))'
    try:
        result=subprocess.run([str(executable),'-I','-c',code],capture_output=True,text=True,timeout=30)
        return json.loads(result.stdout) if result.returncode==0 else None
    except (OSError,ValueError,subprocess.TimeoutExpired):
        return None


def canonical_base_python():
    # uv/python-build-standalone can be reached through ~/.local/bin aliases.
    # Copying that executable with an alias directory as venv.home loses stdlib.
    base=Path(getattr(sys,'_base_executable',sys.executable)).resolve(strict=True)
    info=probe_python(base)
    if not info or info['version']!=[3,12]:
        raise RuntimeError(f'Cannot locate a working Python 3.12 standard library: {base}')
    return base


def ensure_environment(env):
    env=Path(env)
    python=env/('Scripts/python.exe' if os.name=='nt' else 'bin/python')
    info=probe_python(python) if python.is_file() else None
    if info and info['version']==[3,12] and Path(info['prefix']).resolve()==env.resolve():
        return python
    base=canonical_base_python()
    if env.exists():
        if env.is_symlink() or not (env/'pyvenv.cfg').is_file():
            raise RuntimeError(f'Refusing to replace an unrecognized environment directory: {env}; choose a new --prefix')
        backup=env.with_name(env.name+'.broken-'+str(time.time_ns()))
        env.rename(backup)
        print(f'Preserved unusable environment: {backup}',flush=True)
    print(f'Creating environment with base Python: {base}',flush=True)
    subprocess.run([str(base),'-I','-m','venv',str(env)],check=True)
    info=probe_python(python)
    if not info or info['version']!=[3,12] or Path(info['prefix']).resolve()!=env.resolve():
        raise RuntimeError(f'New environment failed its standard-library check: {env}')
    return python


def main():
    parser=argparse.ArgumentParser(description='Install OASIS4.3-GUI-0.1 into an isolated user environment; no Parrot required')
    parser.add_argument('--prefix',type=Path,default=Path.home()/'.oasis-sad')
    parser.add_argument('--no-launch',action='store_true',help='Install without opening the GUI')
    parser.add_argument('--cli-only',action='store_true',help='Install the command-line program without Qt')
    parser.add_argument('--launch-only',action='store_true',help='Launch an already installed GUI')
    parser.add_argument('--reference-source',type=Path,help='Existing pinned reference MTZ for offline reference setup')
    parser.add_argument('--index-url',default='https://pypi.org/simple',help='Python package index URL')
    args=parser.parse_args()
    if sys.version_info[:2]!=(3,12):parser.error('Use Python 3.12: python3.12 install.py (Windows: py -3.12 install.py)')
    root=Path(__file__).resolve().parent;prefix=args.prefix.expanduser().resolve();env=prefix/'venv'
    scripts=env/('Scripts' if os.name=='nt' else 'bin')
    python=scripts/('python.exe' if os.name=='nt' else 'python')
    gui=scripts/('oasis-gui.exe' if os.name=='nt' else 'oasis-gui')
    if not args.launch_only:
        version=tomllib.loads((root/'pyproject.toml').read_text())['project']['version']
        wheels=list((root/'dist').glob(f'oasis4_3_gui-{version}-*.whl'))
        if len(wheels)>1:parser.error(f'Multiple OASIS {version} wheels found in dist/')
        constraints=root/'constraints.txt'
        prefix.mkdir(parents=True,exist_ok=True)
        python=ensure_environment(env)
        subprocess.run([str(python),'-c','import sys; assert sys.version_info[:2]==(3,12), "Choose a fresh --prefix for Python 3.12"'],check=True)
        subprocess.run([str(python),'-m','pip','--isolated','install','--index-url',args.index_url,'--upgrade','pip'],check=True)
        requirement=str(wheels[0] if wheels else root)+('' if args.cli_only else '[gui]')
        subprocess.run([str(python),'-m','pip','--isolated','install','--index-url',args.index_url,'--force-reinstall','--constraint',str(constraints),requirement],check=True)
        subprocess.run([str(python),'-m','pip','check'],check=True)
        provision=[str(python),'-m','oasis_py','reference-fetch']
        if args.reference_source:provision+=['--source',str(args.reference_source.expanduser().resolve())]
        subprocess.run(provision,check=True)
        subprocess.run([str(python),'-c','from oasis_py.reference_data import read_reference; m=read_reference(); print("Default reference ready:", m.nreflections, "reflections")'],check=True)
        print('OASIS4.3-GUI-0.1 installation complete. Command-line executable:',scripts/('oasis-py.exe' if os.name=='nt' else 'oasis-py'))
        print('GUI launcher:',gui)
    if not args.no_launch and not args.cli_only:
        if not gui.is_file():parser.error('GUI not installed; rerun without --launch-only or --cli-only')
        subprocess.Popen([str(gui)],start_new_session=(os.name!='nt'))
    return 0

if __name__=='__main__':
    try:raise SystemExit(main())
    except RuntimeError as error:
        print(f'Installation/check failed: {error}',file=sys.stderr)
        raise SystemExit(1)
    except subprocess.CalledProcessError as error:
        print(f'Installation/check failed (exit {error.returncode}). See output above; correct the cause and rerun.',file=sys.stderr)
        raise SystemExit(error.returncode)
