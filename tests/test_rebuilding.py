from types import SimpleNamespace
import gemmi
import numpy as np
from phasebuild.geometry import template
from phasebuild.trace import write_backbone
from phasebuild.rebuilding import prune_weak_fragments,transfer_properties


def model(tmp_path,group='P 1'):
    d=SimpleNamespace(cell=gemmi.UnitCell(40,40,40,90,90,90),sg=gemmi.find_spacegroup_by_name(group))
    path=tmp_path/'input.cif'
    write_backbone(path,[template(-63.,-42.,length=12)+[20.,20.,20.],
                         template(-63.,-42.,length=3)+[8.,8.,8.]],d)
    s=gemmi.read_structure(str(path))
    for r in s[0][0]:
        r.name='ALA'
        for a in r:a.b_iso=30.
    for r in s[0][1]:
        for a in r:a.b_iso=60.
    return s


def test_pruning_requires_short_unknown_weak_density_and_high_relative_b(tmp_path):
    d=SimpleNamespace(score=lambda xyz:np.full(xyz.shape[:-2],2.5))
    s=model(tmp_path)
    report=prune_weak_fragments(s,d)
    assert report['removed_residues']==3 and len(s[0])==1
    s=model(tmp_path)
    for r in s[0][1]:r.name='ALA'
    assert prune_weak_fragments(s,d)['removed_residues']==0
    s=model(tmp_path)
    for r in s[0][1]:
        for a in r:a.b_iso=35.
    assert prune_weak_fragments(s,d)['removed_residues']==0
    s=model(tmp_path)
    strong=SimpleNamespace(score=lambda xyz:np.full(xyz.shape[:-2],4.))
    assert prune_weak_fragments(s,strong)['removed_residues']==0


def test_property_transfer_uses_whole_residue_crystal_image(tmp_path):
    old=model(tmp_path,'P 21 21 21')
    del old[0][1]
    while len(old[0][0])>1:del old[0][0][-1]
    r=old[0][0][0];r.name='CYS'
    for name,offset in [('CB',[0.,0.,1.5]),('SG',[0.,0.,3.3])]:
        a=gemmi.Atom();a.name=name;a.element=gemmi.Element('S' if name=='SG' else 'C')
        a.pos=gemmi.Position(*(np.array(list(r.find_atom('CA','*').pos))+offset));a.b_iso=44.
        r.add_atom(a)
    new=old.clone();op=list(old.find_spacegroup().operations())[1]
    expected={}
    for a in new[0][0][0]:
        p=old.cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(list(old.cell.fractionalize(a.pos)))))
        a.pos=gemmi.Position(p.x+40.,p.y,p.z-40.);expected[a.name]=np.array(list(a.pos));a.b_iso=80.
    del new[0][0][0][-1];del new[0][0][0][-1]
    report=transfer_properties(old,new)
    assert report['matched_previous_residues']==1
    for a in new[0][0][0]:
        np.testing.assert_allclose(list(a.pos),expected[a.name],atol=1e-6)
        assert a.b_iso==(44. if a.name in ['CB','SG'] else 30.)
