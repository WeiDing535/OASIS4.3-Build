import numpy as np
from phasebuild.local_sequencing import (local_maxima, local_permutation_maxima,
                                         assign_local_sequence)
from phasebuild.residue_profile import ALPHABET


def brute_force_maximum(scores, sequences, minimum_length, penalty):
    scores = np.column_stack([scores, np.zeros(len(scores))])
    best = -np.inf
    for sequence in sequences:
        codes = [ALPHABET.index(c) if c != 'X' else 20 for c in sequence]
        for length in range(minimum_length, min(len(scores), len(codes))+1):
            for begin in range(len(scores)-length+1):
                for start in range(len(codes)-length+1):
                    value = sum(scores[begin+i, codes[start+i]]-penalty for i in range(length))
                    best = max(best, value)
    return best


def test_local_dynamic_search_matches_exhaustive_subsegments_and_sequence_windows():
    scores = np.random.RandomState(8).normal(size=(7, 20))
    sequences = ['ACDEXFGH', 'WYAXCD']
    catalogs = [np.array([[ALPHABET.index(c) if c != 'X' else 20 for c in s]]) for s in sequences]
    actual, candidates = local_maxima(scores, catalogs, 3, .4, True)
    expected = brute_force_maximum(scores, sequences, 3, .4)
    np.testing.assert_allclose(actual, [expected], atol=1e-12)
    assert max(c[0] for c in candidates) == actual[0]
    for value, ri, delta, begin, end in candidates:
        assert end-begin >= 3 and begin+delta >= 0 and end+delta <= len(sequences[ri])


def test_local_null_repeats_full_length_and_offset_search_with_identical_rng_order():
    scores = np.random.RandomState(3).normal(size=(6, 20))
    sequences = ['ACDEXFG', 'YACDE']
    random = np.random.RandomState(11)
    expected = []
    for _ in range(23):
        shuffled = [''.join(random.permutation(list(sequence))) for sequence in sequences]
        expected.append(brute_force_maximum(scores, shuffled, 3, .5))
    actual = local_permutation_maxima(scores, sequences, 23, np.random.RandomState(11), 3, .5, 4)
    np.testing.assert_allclose(actual, expected, rtol=0, atol=1e-12)


class CoreTarget:
    def scores(self, density, chain):
        scores = np.full((len(chain), 20), -3.)
        for i, letter in enumerate('CDEFGHIK'):
            scores[5+i, ALPHABET.index(letter)] = 3.
        return scores


def test_local_assignment_recovers_strong_core_without_assigning_noisy_ends():
    chain = np.arange(18*4*3).reshape(18, 4, 3).astype(float)
    split, report = assign_local_sequence([chain], None, CoreTarget(), [('one', ALPHABET)], shuffles=99)
    accepted = [d for d in report['decisions'] if d['accepted']]
    assert len(accepted) == 1
    anchor = accepted[0]
    assert (anchor['source_begin'], anchor['source_end'], anchor['start'], anchor['end']) == (5, 13, 2, 9)
    assert report['assigned_residues'] == 8
    assert [len(c) for c in split] == [5, 8, 5]
    np.testing.assert_array_equal(np.concatenate(split), chain)
    assert anchor['shuffle_tail_fraction'] <= report['familywise_tail_threshold']


def test_declared_identical_copies_supply_capacity_without_competing_sequence_hypotheses():
    chain = np.zeros((18, 4, 3))
    _, report = assign_local_sequence([chain, chain], None, CoreTarget(),
                                      [('copy_a', ALPHABET), ('copy_b', ALPHABET)], shuffles=99)
    accepted = [d for d in report['decisions'] if d['accepted']]
    assert len(accepted) == 2
    assert {d['sequence_record'] for d in accepted} == {0, 1}
    assert report['assigned_residues'] == 16
    assert report['familywise_tail_threshold'] == .025


def test_no_evidence_abstains_and_preserves_the_entire_trace():
    class NoEvidence:
        def scores(self, density, chain):
            return np.zeros((len(chain), 20))
    chain = np.zeros((12, 4, 3))
    split, report = assign_local_sequence([chain], None, NoEvidence(), [('one', ALPHABET)], shuffles=99)
    assert report['assigned_residues'] == 0
    assert len(split) == 1 and len(split[0]) == len(chain)
    assert report['chain_searches'][0]['best_tail_fraction'] == 1.
