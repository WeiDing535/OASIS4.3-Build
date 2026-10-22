import numpy as np
import pytest
from phasebuild.weak_trace import growth_evidence,supported_intervals

def test_bridges_only_short_troughs_and_preserves_positive_evidence():
 assert supported_intervals([3.,3.,2.1,3.,3.])==[(0,5)]
 assert supported_intervals([3.,3.,2.1,2.1,3.,3.])==[(0,6)]
 assert supported_intervals([3.,3.,1.,3.,3.])==[]
 assert supported_intervals([3.,3.,2.,2.,2.,3.,3.])==[]

def test_negative_growth_evidence_has_bounded_budget():
 gain,weak,ok=growth_evidence(1.5,previous_weak_steps=0);assert gain<0 and weak==1 and ok
 _,weak,ok=growth_evidence(1.5,previous_weak_steps=weak);assert weak==2 and ok
 _,weak,ok=growth_evidence(1.5,previous_weak_steps=weak);assert weak==3 and not ok
 _,weak,ok=growth_evidence(2.5,previous_weak_steps=2);assert weak==0 and ok


def test_recovers_ideal_peptide_across_two_weak_residues_without_coordinate_distortion():
 """A synthetic density trough blocks legacy growth but permits a correct bridge."""
 from phasebuild.geometry import template
 from phasebuild.weak_trace import grow_weak
 from phasebuild.trace import grow
 peptide=template(-63.,-42.,8)
 class SyntheticField:
  profile=True
  def sample(self,xyz):return np.ones(np.asarray(xyz).shape[:-1])*2.
  def distance(self,a,b):return np.linalg.norm(a-b)
  def score(self,xyz):
   xyz=np.asarray(xyz)
   distance=((xyz[...,None,:3,:]-peptide[:,:3])**2).sum(axis=(-1,-2))/3
   peaks=np.array([3.,3.,3.,1.5,1.5,4.,4.,4.])
   return peaks[distance.argmin(axis=-1)]-3.*distance.min(axis=-1)
 field=SyntheticField();seed=peptide[:3].copy()
 hard=grow(field,seed,-42.,8,4,2.)
 recovered=grow_weak(field,seed,8,4)
 assert len(hard)==3 and len(recovered)==8
 np.testing.assert_allclose(recovered[:,:3],peptide[:,:3],atol=1e-5)
