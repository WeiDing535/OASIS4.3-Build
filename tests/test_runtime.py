from pathlib import Path
import pytest
from phasebuild.runtime import resolve_runtime, runtime_environment, default_runtime


def executable(root, name):
    path = root / 'bin' / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text('#!/bin/sh\nexit 0\n')
    path.chmod(0o755)
    return path


def test_explicit_missing_runtime_never_falls_back_to_installed_ccp4(tmp_path, monkeypatch):
    installed = tmp_path / 'installed'
    executable(installed, 'refmac5')
    monkeypatch.setenv('CCP4', str(installed))
    with pytest.raises(ValueError, match='does not exist'):
        resolve_runtime(tmp_path / 'missing', required=('refmac5',))


def test_portable_mode_does_not_discover_host_ccp4(monkeypatch, tmp_path):
    monkeypatch.delenv('PHASEBUILD_RUNTIME', raising=False)
    monkeypatch.setenv('CCP4', str(tmp_path))
    monkeypatch.setenv('PHASEBUILD_STANDALONE', '1')
    assert default_runtime() == ''


def test_environment_replaces_host_library_paths_with_selected_runtime(tmp_path):
    executable(tmp_path, 'refmac5')
    original = {'PATH':'/usr/bin:/bin','CCP4':'/old','CLIBD':'/old/data',
                'CCP4_SCR':'/old/tmp','DYLD_LIBRARY_PATH':'/old/lib',
                'PYTHONHOME':'/old/python','PYTHONPATH':'/old/packages','OMP_NUM_THREADS':'1'}
    result = runtime_environment(tmp_path, original)
    assert result['CCP4'] == str(tmp_path.resolve())
    assert result['CLIBD_MON'] == str(tmp_path / 'lib/data/monomers') + '/'
    assert result['CINCL'] == str(tmp_path / 'include')
    assert result['CCP4_SCR'] != '/old/tmp' and Path(result['CCP4_SCR']).is_dir()
    assert not set(['DYLD_LIBRARY_PATH','PYTHONHOME','PYTHONPATH']) & set(result)
    assert result['PATH'] == original['PATH'] and original['CCP4'] == '/old'


def test_clean_environment_creates_job_scratch_outside_runtime(tmp_path):
    root = tmp_path / 'runtime'; root.mkdir()
    scratch = tmp_path / 'job with spaces' / 'scratch'
    result = runtime_environment(root, {'PATH':'/usr/bin:/bin','LANG':'C'}, scratch=scratch)
    assert scratch.is_dir() and result['CCP4_SCR'] == str(scratch)
    assert result['CINCL'] == str(root / 'include')
    assert not (root / 'scratch').exists()


def test_required_runtime_files_must_be_executable(tmp_path):
    path = executable(tmp_path, 'cbuccaneer')
    assert resolve_runtime(tmp_path, required=('cbuccaneer',)) == tmp_path
    path.chmod(0o644)
    with pytest.raises(ValueError, match='executable unavailable'):
        resolve_runtime(tmp_path, required=('cbuccaneer',))
