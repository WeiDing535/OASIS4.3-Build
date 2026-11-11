"""Numerical controls and separation of map information from likelihood."""
import numpy as np
import pytest
from test_sad_pipeline import inputs
from oasis_py.data import read_mtz
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from oasis_py import native_density as legacy
from phasebuild.dual_channel_dm import density_modify_dual_channel


def setup(inputs):
    data,_=read_mtz(inputs['sad'],'sad',free_label='FreeR_flag',free_value=7)
    hl=np.zeros((len(data.f),4));hl[:,0]=.6;hl[:,2]=1.1
    phase,fom=posterior_moments(data,hl)
    return data,PhaseResult(phase,fom,None,hl,{'method':'synthetic_observation_HL'})


@pytest.mark.parametrize('histogram',['observation','map'])
def test_identical_seed_and_prior_reproduce_preserved_native_engine(inputs,histogram):
    data,prior=setup(inputs)
    expected=legacy.density_modify_native(data,prior,.45,inputs['reference'],cycles=2,seed=17)
    result=density_modify_dual_channel(data,prior,prior,.45,inputs['reference'],cycles=2,seed=17,histogram=histogram)
    # Compare the circular moment, not a non-unique wrapped angle.
    # Equivalent FFT/BLAS evaluation can differ at a few ulps.
    np.testing.assert_allclose(result.fom*np.exp(1j*result.phase),expected.fom*np.exp(1j*expected.phase),atol=1e-12,rtol=0)
    np.testing.assert_allclose(result.hl,expected.hl,atol=1e-12,rtol=1e-12)


def test_map_seed_is_used_without_replacing_or_accumulating_observation_hl(inputs,monkeypatch):
    data,prior=setup(inputs);original=prior.hl.copy();seen=[];first_fc=[]
    def field(data,fc,hl):
        seen.append(hl.copy());assert not hl.flags.writeable
        if len(seen)%3==1:first_fc.append(fc.copy())
        return np.full(len(fc),.2+.1j),[]
    monkeypatch.setattr(legacy,'fit_shell_likelihood',field)
    for shift in [0.,1.3]:
        seed=PhaseResult(prior.phase+shift,prior.fom,None,np.full_like(prior.hl,999.),{})
        result=density_modify_dual_channel(data,prior,seed,.45,inputs['reference'],cycles=3,seed=17)
        expected=original.copy();expected[:,0]+=.2;expected[:,1]+=.1
        np.testing.assert_array_equal(result.hl,expected)
    assert not np.allclose(first_fc[0],first_fc[1])
    assert all(np.array_equal(h,original) for h in seen)
    np.testing.assert_array_equal(prior.hl,original)


def test_contradictory_fixed_likelihood_is_rejected_instead_of_reconstructed(inputs):
    data,prior=setup(inputs);seed=PhaseResult(prior.phase+.5,prior.fom,None,prior.hl,{})
    with pytest.raises(ValueError,match='disagree'):
        density_modify_dual_channel(data,seed,prior,.45,inputs['reference'])
    with pytest.raises(ValueError,match='settings'):
        density_modify_dual_channel(data,prior,seed,.45,inputs['reference'],histogram='automatic')
    seed.fom[0]=np.nan
    with pytest.raises(ValueError,match='map phase/FOM'):
        density_modify_dual_channel(data,prior,seed,.45,inputs['reference'])
