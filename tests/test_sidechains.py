from types import SimpleNamespace
import gemmi
import numpy as np
from phasebuild.chemistry import rotatable_bonds,ring_groups
from phasebuild.sidechains import rotate_chis,retained_atoms,neighbor_atoms,build_sidechains
from phasebuild.profile import frames
from phasebuild.geometry import template
from phasebuild.trace import write_backbone


def aromatic_component():
    edges=[('N','CA'),('CA','C'),('C','O'),('CA','CB'),('CB','CG'),('CG','CD1'),
           ('CD1','CE1'),('CE1','CZ'),('CZ','CE2'),('CE2','CD2'),('CD2','CG')]
    names=['N','CA','C','O','CB','CG','CD1','CD2','CE1','CE2','CZ']
    return {'name':'PHE','elements':{n:'C' for n in names},'sidechain_names':names[4:],
            'bonds':[[a,b,'SINGLE',1.5] for a,b in edges],'angles':[]}


def test_rotors_preserve_ring_geometry_and_pruning_does_not_leave_broken_ring():
    component=aromatic_component()
    rotors=rotatable_bonds(component)
    assert [r[:2] for r in rotors]==[('CA','CB'),('CB','CG')]
    assert set(ring_groups(component)[0])=={'CG','CD1','CD2','CE1','CE2','CZ'}
    names=list(component['elements'])
    raw=np.full(7,3.);raw[names[4:].index('CE1')]=0.
    assert retained_atoms(raw,np.full(7,3.),names,component)==[0]


def test_chi_rotation_preserves_bonds_and_fixed_backbone():
    names=['N','CA','C','O','CB','SG']
    original=np.array([[-1.,1.,0.],[0.,0.,0.],[1.5,0.,0.],[2.,1.,0.],[-.5,-.7,-1.2],[-1.,-2.,-2.]])
    changed=rotate_chis(original,[('CA','CB',['SG'])],names,[1.2])
    np.testing.assert_allclose(changed[:5],original[:5])
    assert np.linalg.norm(changed[5]-original[5])>.5
    np.testing.assert_allclose(np.linalg.norm(changed[5]-changed[4]),np.linalg.norm(original[5]-original[4]))


def test_own_backbone_is_not_mistaken_for_an_external_clash(tmp_path):
    d=SimpleNamespace(cell=gemmi.UnitCell(50,50,50,90,90,90),sg=gemmi.find_spacegroup_by_name('P 1'))
    path=tmp_path/'one.cif'
    xyz=template(-63.,-42.,length=1)+[20.,20.,20.]
    write_backbone(path,[xyz],d)
    structure=gemmi.read_structure(str(path))
    points,_=neighbor_atoms(structure,(0,0),d,xyz[0,1])
    assert len(points)==0


def test_existing_complete_sidechain_can_seed_refitting_and_fallback_counts_are_accurate(tmp_path):
    d=SimpleNamespace(cell=gemmi.UnitCell(50,50,50,90,90,90),sg=gemmi.find_spacegroup_by_name('P 1'))
    path=tmp_path/'seed.cif'
    xyz=template(-63.,-42.,length=1)+[20.,20.,20.]
    write_backbone(path,[xyz],d)
    s=gemmi.read_structure(str(path));res=s[0][0][0];res.name='ALA'
    point=xyz[0,1]+[0.,0.,1.5]
    atom=gemmi.Atom();atom.name='CB';atom.element=gemmi.Element('C');atom.pos=gemmi.Position(*point)
    res.add_atom(atom)
    component={'name':'ALA','elements':{'N':'N','CA':'C','C':'C','O':'O','CB':'C'},
               'sidechain_names':['CB'],'bonds':[['N','CA','SINGLE',1.458],['CA','C','SINGLE',1.525],
                                               ['C','O','DOUBLE',1.231],['CA','CB','SINGLE',1.5]],'angles':[]}
    ca,axes=frames(xyz[0]);bad=np.vstack([xyz[0],ca+[0.,0.,-1.5]])
    library={'entries':{'ALA':{'component':component,'conformers':[((bad-ca)@axes).tolist()]}}}
    d.sample=lambda coordinates:3.*np.exp(-np.sum((coordinates-point)**2,axis=-1)/.2)
    report=build_sidechains(s,d,library,progress=lambda _:None)
    assert report['sidechain_atoms_built']==1 and report['residues'][0]['current_complete_conformer_seed']
    np.testing.assert_allclose(list(s[0][0][0].find_atom('CB','*').pos),point)
    np.testing.assert_allclose([list(s[0][0][0].find_atom(n,'*').pos) for n in ['N','CA','C','O']],xyz[0])
    d.sample=lambda coordinates:np.zeros(coordinates.shape[:-1])
    report=build_sidechains(s,d,library,progress=lambda _:None)
    assert report['sidechain_atoms_built']==1 and report['residues'][0]['built_atoms']==['CB']
