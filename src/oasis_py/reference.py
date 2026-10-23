# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Isolated reference execution; this is an adapter, not a Python rewrite."""
from __future__ import annotations
from pathlib import Path
import os
import subprocess


def run_reference(executable,keywords,workdir,hklin=None,hklout='reference.mtz',
                  model=None,library_dir=None,timeout=600,program='oasis'):
    if program not in {'oasis','sapi'}:raise ValueError('program must be oasis or sapi')
    if program=='sapi' and model is not None:raise ValueError('Use KAR keywords for SAPI model input')
    if timeout<=0:raise ValueError('timeout must be positive')
    executable=Path(executable).resolve();workdir=Path(workdir).resolve()
    if not executable.is_file():raise ValueError(f'Executable not found: {executable}')
    workdir.mkdir(parents=True,exist_ok=True)
    if any(workdir.iterdir()):
        raise ValueError('Reference directory must be empty; legacy scratch files have fixed names')
    command=[str(executable)]
    if hklin is not None:
        command+=['HKLIN',str(Path(hklin).resolve())]
        if program=='oasis':command+=['HKLOUT',str(workdir/hklout)]
        else:command+=['SAPIPKS',str(workdir/'sites.pks')]
    if model is not None:command+=['XYZIN',str(Path(model).resolve())]
    env=os.environ.copy()
    scratch=workdir/'scratch';scratch.mkdir()
    env['CCP4_SCR']=str(scratch)
    if library_dir:
        env['DYLD_LIBRARY_PATH']=str(Path(library_dir).resolve())+os.pathsep+env.get('DYLD_LIBRARY_PATH','')
        env['LD_LIBRARY_PATH']=str(Path(library_dir).resolve())+os.pathsep+env.get('LD_LIBRARY_PATH','')
    keyword_text=Path(keywords).read_text()
    (workdir/'input.key').write_text(keyword_text)
    (workdir/'command.txt').write_text(repr(command)+'\n')
    with (workdir/'reference.log').open('w') as log:
        try:
            result=subprocess.run(command,input=keyword_text,text=True,stdout=log,stderr=subprocess.STDOUT,
                                  cwd=workdir,env=env,timeout=timeout,check=False)
        except subprocess.TimeoutExpired:
            raise RuntimeError(f'Reference exceeded {timeout}s and was terminated; see {workdir}/reference.log')
    if result.returncode:
        raise RuntimeError(f'Reference exit {result.returncode}; see {workdir}/reference.log')
    expected=workdir/(hklout if program=='oasis' else 'sites.pks')
    if hklin is not None and not expected.is_file():
        raise RuntimeError(f'Reference reported success without the expected output: {expected.name}')
    return {'returncode':result.returncode,'directory':str(workdir),
            'output':str(workdir/(hklout if program=='oasis' else 'sites.pks')) if hklin is not None else None}
