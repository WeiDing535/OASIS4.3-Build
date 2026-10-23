"""Independent molecular formulas, Matthews multiplicity and GUI command wiring."""
from pathlib import Path
import numpy as np
import gemmi
import pytest
from oasis_py.sequence import protein_sequence,sequence_asu
from oasis_py.desktop_job import build_sad_command
from test_desktop import job


def write_seq(tmp_path,text):
    path=tmp_path/'protein.fasta';path.write_text(text,encoding='utf-8');return path


def test_peptide_bonds_termini_and_mass_agree_with_gemmi(tmp_path):
    path=write_seq(tmp_path,'>one\nACDM\n');info=protein_sequence(path)
    assert info['composition_per_copy']=={'C':15,'H':26,'N':4,'O':7,'S':2}
    assert info['mass_per_copy_da']==pytest.approx(gemmi.calculate_sequence_weight(['ALA','CYS','ASP','MET']),abs=.005)
    assert info['residues_per_copy']==4


def test_all_twenty_residues_match_independent_sequence_mass(tmp_path):
    seq='ACDEFGHIKLMNPQRSTVWY';info=protein_sequence(write_seq(tmp_path,'>all\n'+seq))
    names=[gemmi.expand_one_letter(x,gemmi.ResidueKind.AA) for x in seq]
    # Gemmi's monomer table uses protonated HIS/LYS/ARG; our formulas are neutral.
    neutral_mass=gemmi.calculate_sequence_weight(names)-3*gemmi.Element('H').weight
    assert info['mass_per_copy_da']==pytest.approx(neutral_mass,abs=.02)


def test_multiple_records_are_chains_not_duplicates(tmp_path):
    info=protein_sequence(write_seq(tmp_path,'>a\nAA\n>b\nAA\n'))
    assert info['composition_per_copy']=={'C':12,'H':24,'N':4,'O':6}
    assert info['chains_per_copy']==2 and info['residues_per_copy']==4


def test_centering_and_full_file_copies_are_applied_once(tmp_path):
    path=write_seq(tmp_path,'>a\nACDM\n>b\nGG\n');info=protein_sequence(path)
    sg=gemmi.SpaceGroup('C 1 2 1');copies=3
    volume=info['mass_per_copy_da']*copies*len(sg.operations())*2.46
    cell=gemmi.UnitCell(volume**(1/3),volume**(1/3),volume**(1/3),90,90,90)
    result=sequence_asu(path,copies,cell,sg)
    assert result['solvent_fraction']==pytest.approx(.5)
    assert result['matthews_vm']==pytest.approx(2.46)
    assert result['asu_chains']==6 and result['asu_residues']==18
    assert result['asu_composition']['C']==info['composition_per_copy']['C']*3
    assert result['symmetry_multiplicity']==4


@pytest.mark.parametrize('text',['ACDM\n','>chain\nacdm\n','>P1;chain\ndescription\nACDM*\n','\ufeff>chain\nACDM\n'])
def test_supported_sequence_text_formats(tmp_path,text):
    assert protein_sequence(write_seq(tmp_path,text))['residues_per_copy']==4


def test_explicit_selenium_labels_do_not_add_native_sulfur(tmp_path):
    info=protein_sequence(write_seq(tmp_path,'>labelled\nM(MSE)(SEC)U\n'))
    assert info['composition_per_copy']['S']==1
    assert info['composition_per_copy']['Se']==3
    assert info['residues_per_copy']==4


@pytest.mark.parametrize('text',['','>empty\n','>a\nAC-X\n','>a\nAXB\n','>a\nAC*DE\n','>a\nAC1D\n','>a\n(SEP)\n'])
def test_unknown_and_aligned_sequences_are_not_silently_approximated(tmp_path,text):
    with pytest.raises(ValueError):protein_sequence(write_seq(tmp_path,text))


@pytest.mark.parametrize('copies',[0,-1,1.5,True])
def test_copies_must_be_positive_integer(tmp_path,copies):
    with pytest.raises(ValueError,match='copies'):sequence_asu(write_seq(tmp_path,'>a\nACDM'),copies,gemmi.UnitCell(50,50,50,90,90,90),gemmi.SpaceGroup('P 1'))


def test_impossible_packing_fails_instead_of_clipping_solvent(tmp_path):
    path=write_seq(tmp_path,'>a\n'+'W'*100)
    with pytest.raises(ValueError,match='不相容'):sequence_asu(path,100,gemmi.UnitCell(10,10,10,90,90,90),gemmi.SpaceGroup('P 1'))


def test_sequence_job_recomputes_solvent_and_has_no_legacy_asu_arguments(tmp_path):
    spec=job(tmp_path);spec.pop('composition');spec.pop('solvent')
    path=write_seq(tmp_path,'>a\n'+'ACDEFGHIKLMNPQRSTVWY'*2)
    spec.update(sequence=str(path),copies=1)
    cmd=build_sad_command(spec)
    assert cmd[cmd.index('--sequence')+1]==str(path)
    assert cmd[cmd.index('--copies')+1]=='1'
    assert '--composition' not in cmd and '--residues' not in cmd
    assert '--solvent-fraction' not in cmd # CLI derives it from the same sequence


def test_gui_job_requires_a_sequence_even_if_old_composition_is_present(tmp_path):
    spec=job(tmp_path);spec.update(sequence='',copies=1)
    with pytest.raises(ValueError,match='序列文件'):build_sad_command(spec)
