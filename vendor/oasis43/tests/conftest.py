"""Ordinary tests use synthetic reference data; never download benchmark targets."""
import json
import hashlib
import numpy as np
import gemmi
import pytest


@pytest.fixture(scope='session')
def synthetic_reference(tmp_path_factory):
    root=tmp_path_factory.mktemp('synthetic-reference')
    m=gemmi.Mtz(with_base=True);m.cell=gemmi.UnitCell(24,26,28,90,90,90);m.spacegroup=gemmi.SpaceGroup('P 1');m.add_dataset('synthetic')
    m.add_column('FC_BEST.F_phi.F','F');m.add_column('FC_BEST.F_phi.phi','P')
    h=gemmi.make_miller_array(m.cell,m.spacegroup,3.)
    f=100*np.exp(-m.cell.calculate_1_d2_array(h)*10)
    phi=(h@np.array([.13,.26,.39]))*360
    m.set_data(np.column_stack([h,f,phi]).astype('float32'));m.update_reso()
    path=root/'reference-1tqw.mtz';m.write_to_file(str(path))
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    (root/'reference-1tqw.json').write_text(json.dumps(dict(pdb_id='SYNTHETIC_TEST_ONLY',mtz_sha256=digest,
        source_sha256=digest,public_data_url='https://example.invalid/not-downloaded',coefficients_changed=False,
        resolution_limit=3.,reflections=m.nreflections)))
    return root


@pytest.fixture(autouse=True)
def isolate_reference(request,monkeypatch,synthetic_reference,tmp_path):
    if request.node.get_closest_marker('reference_data'):return
    import oasis_py.reference_data as module
    monkeypatch.setattr(module,'RESOURCE_DIR',synthetic_reference)
    monkeypatch.setenv('OASIS_REFERENCE_CACHE',str(tmp_path/'reference-cache'))
