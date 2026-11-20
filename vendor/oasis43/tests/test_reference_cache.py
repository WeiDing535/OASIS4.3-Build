import hashlib
import pytest
from oasis_py import reference_data as ref


def test_cached_reference_is_verified_and_does_not_download(monkeypatch):
    monkeypatch.setattr(ref.urllib.request,'urlopen',lambda *a,**k:pytest.fail('Unexpected network access'))
    path=ref.fetch_reference(source=ref.RESOURCE_DIR/'reference-1tqw.mtz')
    assert ref.digest(path)==ref.reference_info()['mtz_sha256']
    assert ref.fetch_reference()==path


def test_corrupt_download_never_replaces_existing_cache(tmp_path):
    path=ref.fetch_reference(source=ref.RESOURCE_DIR/'reference-1tqw.mtz')
    before=path.read_bytes();bad=tmp_path/'bad.mtz';bad.write_bytes(b'not a reference')
    with pytest.raises(ValueError,match='SHA256'):ref.fetch_reference(source=bad,force=True)
    assert path.read_bytes()==before


def test_unprovisioned_reference_fails_with_an_offline_instruction(monkeypatch,tmp_path):
    monkeypatch.setattr(ref,'RESOURCE_DIR',tmp_path)
    with pytest.raises(ValueError,match='reference-fetch'):ref.bundled_reference()


def test_corrupt_cache_is_not_silently_used():
    path=ref.fetch_reference(source=ref.RESOURCE_DIR/'reference-1tqw.mtz');path.write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='corrupt'):ref.fetch_reference()
