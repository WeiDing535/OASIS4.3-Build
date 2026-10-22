from types import SimpleNamespace
import gemmi
from phasebuild.atom_validation import atom_geometry_report
from phasebuild.geometry import template
from phasebuild.trace import write_backbone


def test_unknown_backbone_bonds_are_not_clashes_but_overlap_is(tmp_path):
    density=SimpleNamespace(cell=gemmi.UnitCell(40,40,40,90,90,90),sg=gemmi.find_spacegroup_by_name('P 1'))
    xyz=template(-63.,-42.,length=2)+[20.,20.,20.]
    path=tmp_path/'backbone.cif'
    write_backbone(path,[xyz],density)
    structure=gemmi.read_structure(str(path))
    report=atom_geometry_report(structure,{'entries':{}})
    assert report['heavy_atoms']==8
    assert not report['severe_nonbonded_contacts_below_1p8A']
    atom=gemmi.Atom();atom.name='CB';atom.element=gemmi.Element('C')
    atom.pos=gemmi.Position(*(xyz[0,1]+[.5,0.,0.]))
    structure[0][0][0].add_atom(atom)
    report=atom_geometry_report(structure,{'entries':{}})
    assert report['severe_contacts_involving_sidechains']>0


def test_audit_includes_crystal_mates(tmp_path):
    density=SimpleNamespace(cell=gemmi.UnitCell(40,40,40,90,90,90),sg=gemmi.find_spacegroup_by_name('P -1'))
    xyz=template(-63.,-42.,length=1)
    xyz+=-xyz[0,1]+[.4,0.,0.]
    path=tmp_path/'mate.cif'
    write_backbone(path,[xyz],density)
    report=atom_geometry_report(gemmi.read_structure(str(path)),{'entries':{}})
    assert any(c['symmetry_image']>0 for c in report['severe_nonbonded_contacts_below_1p8A'])
