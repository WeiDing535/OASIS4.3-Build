import numpy as np
import pytest
from phasebuild.residue_profile import ALPHABET
from phasebuild.sequencing import read_sequences, assign_sequence


class KnownScores:
    def __init__(self,letters):
        self.letters=letters

    def scores(self,density,chain):
        matrix=np.zeros((len(chain),20))
        for i,letter in enumerate(self.letters):
            matrix[i,ALPHABET.index(letter)]=3.
        return matrix


def test_parse_fasta_and_true_pir_without_treating_description_as_sequence(tmp_path):
    path=tmp_path/'seq.pir'
    path.write_text('>P1;first\nprotein description\nACDEFG*\n>second\n\nHIKLMN\n')
    assert read_sequences(path)==[('P1;first','ACDEFG'),('second','HIKLMN')]
    path.write_text('Run this command instead of reading a sequence\n')
    with pytest.raises(ValueError,match='header'):
        read_sequences(path)


def test_confident_sequence_placement_and_no_reused_sequence_interval():
    chain=np.zeros((6,4,3))
    report=assign_sequence([chain,chain],None,KnownScores('CDEFGH'),[('seq',ALPHABET)],shuffles=99)
    accepted=[d for d in report['decisions'] if d['accepted']]
    assert len(accepted)==1
    assert accepted[0]['start']==2 and accepted[0]['end']==7
    assert report['assigned_residues']==6


def test_ambiguous_repeated_sequence_is_left_unknown():
    chain=np.zeros((6,4,3))
    report=assign_sequence([chain],None,KnownScores('CDEFGH'),[('seq','CDEFGHCDEFGH')],shuffles=99)
    assert report['assigned_residues']==0
    assert report['decisions'][0]['runner_up_margin']==0.


def test_identical_asu_copies_do_not_create_false_sequence_ambiguity():
    chain=np.zeros((6,4,3))
    report=assign_sequence([chain,chain,chain],None,KnownScores('CDEFGH'),
                           [('copyA',ALPHABET),('copyB',ALPHABET)],shuffles=199)
    accepted=[d for d in report['decisions'] if d['accepted']]
    assert len(accepted)==2
    assert {d['sequence_record'] for d in accepted}=={0,1}
    assert len({d['fragment_index'] for d in accepted})==2
    assert all(d['runner_up_margin']>0 for d in accepted)
    assert report['assigned_residues']==12


def test_declared_copies_do_not_resolve_within_sequence_repeat_ambiguity():
    chain=np.zeros((6,4,3))
    report=assign_sequence([chain],None,KnownScores('CDEFGH'),
                           [('copyA','CDEFGHCDEFGH'),('copyB','CDEFGHCDEFGH')],shuffles=99)
    assert report['assigned_residues']==0
