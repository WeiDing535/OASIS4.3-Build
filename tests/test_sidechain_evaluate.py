import gemmi
import numpy as np
from phasebuild.sidechain_evaluate import evaluate_sidechains
from phasebuild.geometry import template
from phasebuild.trace import write_backbone
from types import SimpleNamespace


def test_carboxyl_oxygen_names_are_equivalent_without_excusing_wrong_residue_type(tmp_path):
    d=SimpleNamespace(cell=gemmi.UnitCell(50,50,50,90,90,90),sg=gemmi.find_spacegroup_by_name('P 1'))
    ref,built=tmp_path/'ref.cif',tmp_path/'built.cif'
    write_backbone(ref,[template(-63.,-42.,length=1)+[20.,20.,20.]],d)
    s=gemmi.read_structure(str(ref));r=s[0][0][0];r.name='ASP'
    for name,pos in [('CB',[19.,19.,19.]),('CG',[18.,18.,18.]),('OD1',[17.,18.,18.]),('OD2',[18.,17.,18.])]:
        a=gemmi.Atom();a.name=name;a.element=gemmi.Element('O' if name.startswith('O') else 'C');a.pos=gemmi.Position(*pos)
        r.add_atom(a)
    s.make_mmcif_document().write_file(str(ref))
    a,b=r.find_atom('OD1','*'),r.find_atom('OD2','*')
    p,q=list(a.pos),list(b.pos);a.pos=gemmi.Position(*q);b.pos=gemmi.Position(*p)
    s.make_mmcif_document().write_file(str(built))
    report=evaluate_sidechains(built,ref)
    assert report['sidechain_atoms_within_cutoff']==4
    assert report['fully_correct_complete_sidechains']==1
    r.name='ASN';s.make_mmcif_document().write_file(str(built))
    report=evaluate_sidechains(built,ref)
    assert report['sidechain_atoms_within_cutoff']==0
