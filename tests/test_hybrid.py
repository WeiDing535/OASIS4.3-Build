from types import SimpleNamespace
import gemmi
import numpy as np
from phasebuild.geometry import template
from phasebuild.trace import write_backbone
from run_hybrid import model_fragments


def test_external_model_is_split_at_a_broken_peptide(tmp_path):
    density = SimpleNamespace(cell=gemmi.UnitCell(40,40,40,90,90,90),
                              sg=gemmi.find_spacegroup_by_name('P 1'))
    xyz = template(-63.,-42.,length=4) + [10.,10.,10.]
    xyz[-1] += [10.,0.,0.]
    path = tmp_path / 'external.cif'
    write_backbone(path, [xyz], density)
    pieces = model_fragments(path, density)
    assert len(pieces) == 1
    assert len(pieces[0]) == 3
    np.testing.assert_allclose(pieces[0],xyz[:3],atol=1e-6)


def test_buccaneer_n_ca_c_only_input_is_not_silently_discarded(tmp_path):
    density = SimpleNamespace(cell=gemmi.UnitCell(40,40,40,90,90,90),
                              sg=gemmi.find_spacegroup_by_name('P 1'),
                              sample=lambda points: np.ones(points.shape[:-1]))
    xyz = template(-63.,-42.,length=4) + [10.,10.,10.]
    path = tmp_path / 'without_oxygen.cif'
    write_backbone(path, [xyz], density)
    structure = gemmi.read_structure(str(path))
    for residue in structure[0][0]:
        del residue[3]
    structure.make_mmcif_document().write_file(str(path))
    pieces = model_fragments(path, density)
    assert len(pieces) == 1 and len(pieces[0]) == 4
    assert np.isfinite(pieces[0]).all()
    np.testing.assert_allclose(pieces[0][:,:3],xyz[:,:3],atol=1e-6)
    np.testing.assert_allclose(pieces[0][:-1,3],xyz[:-1,3],atol=1e-6)
