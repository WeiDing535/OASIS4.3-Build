from types import SimpleNamespace
import gemmi
import numpy as np
from phasebuild.geometry import template
from phasebuild.completion import complete_internal_gaps
from phasebuild.validation import geometry_report


def density(score=3.):
    s = gemmi.Structure()
    s.cell = gemmi.UnitCell(100,100,100,90,90,90)
    s.spacegroup_hm = 'I 21 3'
    s.setup_cell_images()
    def distance(a,b):
        return s.cell.find_nearest_image(gemmi.Position(*a),gemmi.Position(*b)).dist()
    return SimpleNamespace(cell=s.cell,score=lambda xyz:np.full(xyz.shape[:-2],score),
                           sample=lambda xyz:np.ones(xyz.shape[:-1])*2.,distance=distance)


def test_recover_gap_with_one_symmetry_registration_for_entire_fragment():
    d = density()
    parent = template(-63.,-42.,length=8) + [30.,34.,41.]
    op = list(gemmi.find_spacegroup_by_name('I 21 3').operations())[7]
    transformed = np.array([list(d.cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(
        list(d.cell.fractionalize(gemmi.Position(*p))))))) for p in parent[5:].reshape(-1,3)]).reshape(3,4,3)
    transformed += [100.,0.,0.]
    chains,report = complete_internal_gaps([parent[:3],transformed],[parent],d)
    assert len(chains)==1 and len(chains[0])==8
    assert report['added_residues']==2
    np.testing.assert_allclose(chains[0],parent,atol=1e-6)
    assert not geometry_report(chains,d)['severe_geometry_issues']


def test_weak_density_or_unsupported_gap_is_not_filled():
    parent = template(-63.,-42.,length=8)+[30.,34.,41.]
    chains,report = complete_internal_gaps([parent[:3],parent[5:]],[parent],density(score=1.))
    assert len(chains)==2 and report['added_residues']==0
    chains,report = complete_internal_gaps([parent[:3],parent[5:]],[parent],density(),max_gap=1)
    assert len(chains)==2 and report['added_residues']==0
