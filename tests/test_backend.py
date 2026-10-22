import json
from types import SimpleNamespace
import pytest
from phasebuild import backends


def test_external_failure_is_recorded_and_not_reported_as_success(tmp_path, monkeypatch):
    ccp4 = tmp_path / 'ccp4'
    binary = ccp4 / 'bin/cbuccaneer'
    binary.parent.mkdir(parents=True)
    binary.write_text('fake binary')
    binary.chmod(0o755)
    lib = ccp4 / 'lib/data/reference_structures'
    lib.mkdir(parents=True)
    for extension in ['pdb','mtz']:
        (lib / ('reference-1tqw.' + extension)).write_text('fixture')
    monkeypatch.setattr(backends, 'audit', lambda _: {'input':{'path':'fixture'},
                        'columns':[{'label':x} for x in ['F','SIGF','FWT','PHWT','FOM']]})
    monkeypatch.setattr(backends.subprocess, 'run', lambda *args, **kwargs: SimpleNamespace(returncode=7))
    out = tmp_path / 'run'
    with pytest.raises(RuntimeError, match='Buccaneer failed'):
        backends.buccaneer('fixture.mtz', out, str(ccp4))
    report = json.loads((out / 'build_report.json').read_text())
    assert report['status'] == 'failed'
    assert report['returncode'] == 7
    assert not (out / 'backbone.cif').exists()
