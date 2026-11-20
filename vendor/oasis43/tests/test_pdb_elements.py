import gemmi
import pytest
from oasis_py.data import read_atoms, read_structure_checked
from test_science import dataset


def test_shifted_selenium_is_repaired_from_input_evidence_without_rewriting(tmp_path):
    path=tmp_path/'sites.pdb'
    # Real export defect: Se occupies columns 78–79, standard parser sees S.
    line='HETATM  473 SE   MSE A  78      -8.121  28.502  17.912  1.00 195.22          Se\n'
    path.write_text(line+'END\n');before=path.read_bytes()
    assert gemmi.read_structure(str(path))[0][0][0][0].element.name=='S'
    with pytest.warns(UserWarning,match='S -> Se'):
        st,repairs=read_structure_checked(path)
    assert st[0][0][0][0].element.name=='Se'
    assert len(repairs)==1 and repairs[0]['element']=='Se'
    assert path.read_bytes()==before
    with pytest.warns(UserWarning):
        atoms=read_atoms(path,dataset())
    assert atoms.elements==['Se']


def test_valid_two_letter_element_and_protein_ca_are_unchanged(tmp_path):
    path=tmp_path/'valid.pdb'
    path.write_text('HETATM  473 SE   MSE A  78      -8.121  28.502  17.912  1.00195.22          Se\n'
                    'ATOM    474  CA  ALA A  79       1.000   2.000   3.000  1.00 20.00           C\nEND\n')
    st,repairs=read_structure_checked(path)
    assert repairs==[]
    assert [a.element.name for c in st[0] for r in c for a in r]==['Se','C']


def test_conflicting_atom_name_is_not_used_to_guess_chemistry(tmp_path):
    path=tmp_path/'ambiguous.pdb'
    path.write_text('HETATM  473 SG   CYS A  78      -8.121  28.502  17.912  1.00 195.22          Se\nEND\n')
    st,repairs=read_structure_checked(path)
    assert repairs==[] and st[0][0][0][0].element.name=='S'
