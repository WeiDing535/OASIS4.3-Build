from types import SimpleNamespace
import numpy as np
from scipy.special import logit
from oasis_py.relations import Relations
from phasebuild.pplus_stability import relation_groups,stability_shrink
from phasebuild.pplus_trust import bounded_branch_mass,bernoulli_kl

def test_shrink_uses_agreement_and_rejects_sparse_or_opposed_updates():
 b=np.full(4,.5);p,a=stability_shrink(b,np.full(4,.6),[.6,.6,.6,.55],[.6,.4,.6,.6],[4,4,3,4],[4]*4)
 np.testing.assert_allclose(p,[.6,.5,.5,.5+.1*2/3]);np.testing.assert_allclose(a,[1,0,0,2/3])

def test_convex_shrink_respects_original_probability_support_and_KL():
 rng=np.random.default_rng(4);b=rng.uniform(.001,.999,200);lm=np.log(np.column_stack([b,1-b]));full=bounded_branch_mass(lm,rng.normal(0,10,200),.0005);p,a=stability_shrink(b,full,full,b+(full-b)*.5,np.full(200,8),np.full(200,8))
 assert np.max(bernoulli_kl(p,lm))<=.0005+1e-12;assert np.all(abs(p-b)<=abs(full-b))

def test_equivalent_triplets_share_group_and_HKL_row_order_is_irrelevant():
 h=np.array([[1,1,0],[1,0,0],[0,1,0]])
 graph=Relations(np.array([0,0,1]),np.array([1,2,0]),np.array([2,1,2]),np.array([1,1,1]),np.array([1,1,-1]),np.zeros(3),np.ones(3),3,{})
 groups=relation_groups(SimpleNamespace(hkl=h),graph);assert len(set(groups))==1
 order=np.array([2,0,1]);inverse=np.argsort(order);other=Relations(inverse[graph.target],inverse[graph.left],inverse[graph.right],graph.left_sign,graph.right_sign,graph.shift,graph.kappa,3,{})
 assert relation_groups(SimpleNamespace(hkl=h[order]),other)==groups

def test_phase_shift_global_sign_equivalence_is_grouped():
 h=np.array([[1,1,0],[1,0,0],[0,1,0]]);graph=Relations(np.array([0,1]),np.array([1,0]),np.array([2,2]),np.array([1,1]),np.array([1,-1]),np.array([np.pi/3,-np.pi/3]),np.ones(2),3,{})
 assert len(set(relation_groups(SimpleNamespace(hkl=h),graph)))==1
