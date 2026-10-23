"""Regression coverage for standalone Python entered through a launcher alias."""
import importlib.util
from pathlib import Path
import os
import subprocess
import sys
import pytest

SOURCE=Path(__file__).resolve().parents[1]/'install.py'
spec=importlib.util.spec_from_file_location('oasis_local_installer',SOURCE)
installer=importlib.util.module_from_spec(spec);spec.loader.exec_module(installer)


def test_alias_entry_creates_a_working_environment(tmp_path):
    if os.name=='nt':pytest.skip('POSIX launcher symlink regression')
    alias=tmp_path/'launchers/python3.12';alias.parent.mkdir()
    alias.symlink_to(Path(sys._base_executable).resolve())
    env=tmp_path/'created'
    code=('import importlib.util; from pathlib import Path; '
          f's=importlib.util.spec_from_file_location("installer",{str(SOURCE)!r}); '
          'm=importlib.util.module_from_spec(s); s.loader.exec_module(m); '
          f'm.ensure_environment(Path({str(env)!r}))')
    subprocess.run([str(alias),'-I','-c',code],check=True,capture_output=True,text=True)
    info=installer.probe_python(env/'bin/python')
    assert info and Path(info['prefix']).resolve()==env.resolve()
    cfg=(env/'pyvenv.cfg').read_text()
    assert str(alias.parent) not in cfg


def test_broken_environment_is_backed_up_before_repair(tmp_path):
    env=tmp_path/'venv';env.mkdir();(env/'pyvenv.cfg').write_text('home = /install\n')
    (env/'preserve.txt').write_text('existing environment data')
    python=installer.ensure_environment(env)
    assert installer.probe_python(python)
    backups=list(tmp_path.glob('venv.broken-*'));assert len(backups)==1
    assert (backups[0]/'preserve.txt').read_text()=='existing environment data'


def test_unrecognized_directory_is_preserved(tmp_path):
    env=tmp_path/'venv';env.mkdir();(env/'keep.txt').write_text('untouched')
    with pytest.raises(RuntimeError,match='unrecognized'):installer.ensure_environment(env)
    assert (env/'keep.txt').read_text()=='untouched'


def test_working_environment_is_reused(tmp_path):
    env=tmp_path/'venv';python=installer.ensure_environment(env)
    (env/'keep.txt').write_text('keep installed packages')
    assert installer.ensure_environment(env)==python
    assert (env/'keep.txt').exists() and not list(tmp_path.glob('venv.broken-*'))
