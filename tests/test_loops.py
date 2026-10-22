from types import SimpleNamespace
import numpy as np
from phasebuild.geometry import template
from phasebuild.loops import loop_candidates,sequence_hypothesis
from phasebuild.polish import rotate


def fake_density(score=3.2):
    return SimpleNamespace(score=lambda xyz:np.full(xyz.shape[:-2],score),
                           sample=lambda xyz:np.full(xyz.shape[:-1],3.))


def test_de_novo_loop_closes_without_a_parent_trace_and_weak_map_is_rejected():
    ideal=template(-63.,-42.,length=9)+[20.,20.,20.]
    left,right=ideal[:3],ideal[4:]
    candidates,_=loop_candidates(left,right,1,fake_density(),None,'A',starts=3,right_letter='A')
    assert candidates
    best=max(candidates,key=lambda c:c['quality'])
    assert len(best['xyz'])==len(ideal)
    assert best['closure_rmsd']<.25
    np.testing.assert_allclose(best['xyz'][:3,:3],left[:,:3],atol=1e-6)
    np.testing.assert_allclose(best['xyz'][-len(right):],right,atol=1e-6)
    candidates,_=loop_candidates(left,right,1,fake_density(0.),None,'A',starts=3,right_letter='A')
    assert not candidates


def test_bounded_flexibility_keeps_outer_anchor_fixed():
    ideal=template(-63.,-42.,length=9)+[20.,20.,20.]
    left,right=ideal[:3],ideal[4:].copy()
    right[0,0]=rotate(right[0,0:1],right[0,1],right[0,2]-right[0,1],[45.])[0,0]
    candidates,_=loop_candidates(left,right,1,fake_density(),None,'A',starts=3,
                                 flexible_right=2,flexible_sequence='AA',right_letter='A')
    assert candidates
    best=max(candidates,key=lambda c:c['quality'])
    assert best['maximum_anchor_atom_displacement']<=1.4
    np.testing.assert_allclose(best['xyz'][-3:],right[2:],atol=1e-6)


def test_sequence_order_does_not_force_a_join_between_incompatible_records():
    left={'record':0,'start':0,'end':6};right={'record':1,'start':8,'end':15}
    assert sequence_hypothesis(6,7,left,right,2,[('a','A'*20),('b','A'*20)]) is None
    assert sequence_hypothesis(6,7,left,{'record':0,'start':9,'end':16},2,[('a','A'*20)]) is None
