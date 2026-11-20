import gemmi
import pytest
from phasebuild.experiment_chemistry import apply_semet


def fixture_structure():
    st = gemmi.Structure(); model = gemmi.Model('1'); chain = gemmi.Chain('A')
    for index, name in enumerate(['MET', 'MSE', 'CYS', 'MET']):
        residue = gemmi.Residue(); residue.name=name; residue.seqid=gemmi.SeqId(index+1,' ')
        names = ['N','CA','C'] + (['SD'] if name=='MET' and index==0 else ['SE'] if name=='MSE' else ['SG'] if name=='CYS' else [])
        for i, atom_name in enumerate(names):
            atom=gemmi.Atom();atom.name=atom_name;atom.element=gemmi.Element('Se' if atom_name=='SE' else 'S' if atom_name in ['SD','SG'] else 'N' if atom_name=='N' else 'C');atom.pos=gemmi.Position(index*10+i,.25,1.75);atom.occ=.85;atom.b_iso=32.5;residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain);st.add_model(model);return st


def test_declared_semet_changes_only_met_identity_and_preserves_every_coordinate():
    original=fixture_structure();converted,report=apply_semet(original,True)
    assert [r.name for r in original[0][0]]==['MET','MSE','CYS','MET']
    assert [r.name for r in converted[0][0]]==['MSE','MSE','CYS','MSE']
    assert report['converted_met_residues']==2 and report['converted_sd_atoms']==1
    assert report['converted_residues_without_sd']==1
    assert converted[0][0][0][-1].name=='SE' and converted[0][0][0][-1].element.name=='Se'
    for before,after in zip(original[0][0],converted[0][0]):
        assert len(before)==len(after) and str(before.seqid)==str(after.seqid)
        for a,b in zip(before,after):
            assert a.pos.dist(b.pos)==0 and a.occ==b.occ and a.b_iso==b.b_iso
    assert converted[0][0][2][-1].name=='SG'


def test_conversion_is_explicit_and_idempotent():
    original=fixture_structure();off,report=apply_semet(original,False)
    assert off.make_mmcif_document().as_string()==original.make_mmcif_document().as_string()
    converted,_=apply_semet(original,True);again,report=apply_semet(converted,True)
    assert again.make_mmcif_document().as_string()==converted.make_mmcif_document().as_string()
    assert report['converted_met_residues']==0
    with pytest.raises(ValueError,match='boolean'):apply_semet(original,'Se')
