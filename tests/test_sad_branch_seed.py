import numpy as np
import gemmi
import pytest
from oasis_py.data import Reflections
from oasis_py.phasing import sad_hl
from oasis_py.sad import posterior_moments
from phasebuild.sad_branch_seed import branch_weighted_seed


def data(n,sg='P 1'):
    return Reflections([[i+1,0,1] for i in range(n)],np.ones(n),np.ones(n),gemmi.UnitCell(30,30,30,90,90,90),gemmi.SpaceGroup(sg))


def test_balanced_branches_recover_symmetric_observation_moment():
    d=data(3);ref=np.array([.1,.7,-1.2]);hl=sad_hl(ref,[2,3,4],[1,2,-2],[1,1,1]);before=hl.copy()
    seed=branch_weighted_seed(d,hl,ref,[.5]*3)
    phase,fom=posterior_moments(d,hl)
    np.testing.assert_allclose(seed.fom*np.exp(1j*seed.phase),fom*np.exp(1j*phase),atol=2e-7)
    np.testing.assert_array_equal(hl,before);assert seed.hl is None


def test_branch_weight_changes_sign_and_matches_independent_dense_integration():
    d=data(2);ref=np.array([.3,.3]);hl=sad_hl(ref,[2,2],[1,1],[.8,.8]);p=np.array([.85,.15])
    seed=branch_weighted_seed(d,hl,ref,p,tolerance=1e-7)
    theta=(np.arange(65536)+.5)*2*np.pi/65536+.3
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)]);lp=hl[0]@basis;w=np.exp(lp-lp.max());positive=np.sin(theta-.3)>0
    zp=np.sum(w[positive]*np.exp(1j*theta[positive]))/sum(w[positive]);zm=np.sum(w[~positive]*np.exp(1j*theta[~positive]))/sum(w[~positive])
    np.testing.assert_allclose(seed.fom*np.exp(1j*seed.phase),p*zp+(1-p)*zm,atol=2e-7)
    assert np.imag(np.exp(-.3j)*seed.fom*np.exp(1j*seed.phase))[0]>0


def test_centric_branch_probability_stays_on_allowed_states():
    d=data(2,'P 21');seed=branch_weighted_seed(d,np.zeros((2,4)),np.full(2,np.pi/2),[.8,.2])
    np.testing.assert_allclose(seed.fom,[.6,.6]);np.testing.assert_allclose(np.cos(seed.phase),[-1,1])


def test_invalid_probability_or_nonconverged_grid_is_explicit():
    d=data(1)
    with pytest.raises(ValueError):branch_weighted_seed(d,np.zeros((1,4)),[0],[1.1])
    with pytest.raises(RuntimeError,match='did not converge'):
        branch_weighted_seed(d,np.array([[1.,0.,-2.,0.]]),[0],[.9],tolerance=1e-14,max_samples=1440)
