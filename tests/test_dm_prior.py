"""Guard against comparing different phase priors as a DM-engine effect."""
import gemmi
import numpy as np
import pytest
from oasis_py.data import Reflections,write_mtz
from oasis_py.phasing import PhaseResult
from oasis_py.sad import posterior_moments
from phasebuild.dm_prior import prepare_dm_prior
from phasebuild.sad_input import read_phase_prior
from phasebuild.pplus_experimental import bimodal_prior


def observations():
    return Reflections([[3,1,1],[1,1,1],[2,1,1]],[30,10,20],[1,1,1],
        gemmi.UnitCell(30,30,30,90,90,90),gemmi.SpaceGroup('P 1'))


def test_conflicting_pplus_hl_and_phase_require_explicit_handoff(tmp_path):
    data=observations();hl=np.zeros((3,4));hl[:,2]=2
    initial=PhaseResult(np.array([.3,.8,1.1]),np.array([.4,.6,.7]),np.full(3,.8),hl,{'method':'pplus'})
    raw=hl.copy();obs=prepare_dm_prior(data,initial,representation='observation_hl')
    posterior=prepare_dm_prior(data,initial,representation='phase_posterior')
    assert obs.fom.max()<1e-12
    np.testing.assert_allclose(posterior.fom*np.exp(1j*posterior.phase),initial.fom*np.exp(1j*initial.phase),atol=1e-6)
    assert 'von_Mises' in posterior.metadata['effective_distribution']
    np.testing.assert_array_equal(initial.hl,raw)
    path=tmp_path/'shared.mtz';write_mtz(path,data,posterior.phase,posterior.fom,posterior.hl)
    stored=read_phase_prior(path,data)
    phase,fom=posterior_moments(data,stored.hl)
    np.testing.assert_allclose(stored.fom*np.exp(1j*stored.phase),fom*np.exp(1j*phase),atol=1e-6)


def test_consistent_bimodal_posterior_is_not_collapsed():
    data=observations();initial=bimodal_prior(data,np.zeros(3),np.full(3,np.pi/2),np.full(3,.5),np.full(3,.9))
    result=prepare_dm_prior(data,initial,representation='phase_posterior')
    assert result.fom.max()<1e-7
    assert np.all(result.hl[:,2]<-1)
    assert result.metadata['effective_distribution']=='original_consistent_HL'


def test_missing_or_ambiguous_prior_choice_fails():
    data=observations();initial=PhaseResult(np.zeros(3),np.zeros(3),None,None,{})
    with pytest.raises(TypeError):prepare_dm_prior(data,initial)
    with pytest.raises(ValueError,match='Observation HL required'):prepare_dm_prior(data,initial,representation='observation_hl')
    with pytest.raises(ValueError,match='Choose'):prepare_dm_prior(data,initial,representation='auto')
