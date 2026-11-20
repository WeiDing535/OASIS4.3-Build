import numpy as np
from phasebuild.sequence_null import permutation_maxima,required_shuffles
from phasebuild.sequencing import windows,placement_scores

def test_vectorized_null_matches_composition_shuffle_window_search_exactly():
    random=np.random.RandomState(4);scores=random.normal(size=(6,20));sequences=['ACDEFGHIKLMN','QRSTVWYACDEFG','ACD']
    reference_random=np.random.RandomState(7);expected=[]
    for _ in range(39):
        decoys=[(str(i),''.join(reference_random.permutation(list(s)))) for i,s in enumerate(sequences)]
        expected.append(float(placement_scores(scores,windows(decoys,6)).max()))
    actual=permutation_maxima(scores,sequences,39,np.random.RandomState(7),batch_size=7)
    np.testing.assert_array_equal(actual,expected)


def test_tail_resolution_can_reach_multiple_fragment_threshold_without_relaxing_it():
    assert required_shuffles(499,2)==499
    assert required_shuffles(499,30)==599
    assert 1/(required_shuffles(499,80)+1)<=.05/80


def test_fast_sequence_placement_preserves_legacy_decisions_for_attainable_thresholds():
    from phasebuild.fast_sequencing import assign_sequence as fast
    from phasebuild.sequencing import assign_sequence as legacy
    from phasebuild.residue_profile import ALPHABET
    class Target:
        def scores(self,density,chain):
            result=np.zeros((len(chain),20))
            for i,x in enumerate('CDEFGH'):result[i,ALPHABET.index(x)]=3
            return result
    chains=[np.zeros((6,4,3))]*2;records=[('a',ALPHABET),('b',ALPHABET)]
    a=legacy(chains,None,Target(),records,shuffles=99);b=fast(chains,None,Target(),records,shuffles=99)
    assert a['decisions']==b['decisions'] and a['assigned_residues']==b['assigned_residues']
