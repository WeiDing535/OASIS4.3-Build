"""Observation/seed separation and frozen-state contracts for opt-in P+."""
import numpy as np
import gemmi
import pytest
from oasis_py.phasing import sad_hl
from oasis_py.data import structure_factors,write_mtz
from phasebuild import pplus_calibrated as module
from test_pplus_observed_scale import synthetic,COMPOSITION,FPP,OPTIONS


def prepare(monkeypatch,with_centric=False,placeholder=1.):
    data,heavy,single,_,raw=synthetic(with_centric=with_centric,placeholder=placeholder)
    monkeypatch.setattr(module,'fit_phase_marginal_scale',
        lambda *a,**kw:(.23,{'test_known_fixed_scale':.23}))
    context=module.prepare_observation_pplus(data,COMPOSITION,heavy,FPP,
        centric_singletons=single,**OPTIONS)
    return data,heavy,single,raw,context


def test_reported_sigmas_define_both_prior_and_conditional_seed(monkeypatch):
    data,heavy,single,raw,context=prepare(monkeypatch)
    expected=sad_hl(np.angle(raw),.23*abs(raw),data.dano,data.sigdano)
    np.testing.assert_array_equal(context.observations,expected)
    prior=module.observation_pplus_prior(data,context)
    seed=module.phase_observation_pplus(data,context,cycles=1)
    np.testing.assert_allclose(seed.fom*np.exp(1j*seed.phase),
        prior.fom*np.exp(1j*prior.phase),atol=2e-14,rtol=0)
    assert context.metadata['joint_bijvoet_likelihood'] is False
    assert context.metadata['fully_calibrated_SAD_model'] is False


def test_seed_cannot_be_mistaken_for_experimental_HL(monkeypatch,tmp_path):
    data,_,_,_,context=prepare(monkeypatch)
    before=context.observations.tobytes()
    prior=module.observation_pplus_prior(data,context)
    seed=module.phase_observation_pplus(data,context)
    assert seed.hl is None and seed.probability is not None
    assert prior.probability is None
    prior.hl[0,0]+=1
    assert context.observations.tobytes()==before
    path=tmp_path/'seed.mtz'
    write_mtz(path,data,seed.phase,seed.fom,seed.hl,seed.probability)
    labels=gemmi.read_mtz_file(str(path)).column_labels()
    assert 'PPLUS' in labels and not any(k in labels for k in ['HLA','HLB','HLC','HLD'])


def test_singleton_placeholder_cannot_change_observation_or_seed(monkeypatch):
    a=prepare(monkeypatch,True,1.)
    b=prepare(monkeypatch,True,1e12)
    outputs=[module.phase_observation_pplus(x[0],x[-1]) for x in [a,b]]
    np.testing.assert_array_equal(a[-1].observations[a[2]],0.)
    # Independent direct-summation calls can differ by one floating-point ULP.
    np.testing.assert_allclose(a[-1].observations,b[-1].observations,atol=1e-13,rtol=1e-14)
    np.testing.assert_allclose(outputs[0].fom*np.exp(1j*outputs[0].phase),
        outputs[1].fom*np.exp(1j*outputs[1].phase),atol=1e-11,rtol=0)


def test_changed_data_graph_or_metadata_rejects_reuse(monkeypatch):
    data,_,_,_,context=prepare(monkeypatch)
    for a in [context.observations,context.sim,context.graph.kappa]:
        assert not a.flags.writeable
    data.dano[0]+=1
    with pytest.raises(ValueError,match='Frozen'):module.phase_observation_pplus(data,context)
    data.dano[0]-=1
    context.metadata['phase_samples']=1440
    with pytest.raises(ValueError,match='Frozen'):module.observation_pplus_prior(data,context)


def test_invalid_singleton_class_rejected_before_scale_fit(monkeypatch):
    data,heavy,single,_,_=synthetic()
    def forbidden(*a,**kw):raise AssertionError('Reached fit')
    monkeypatch.setattr(module,'fit_phase_marginal_scale',forbidden)
    single[0]=True
    with pytest.raises(ValueError,match='singleton'):
        module.prepare_observation_pplus(data,COMPOSITION,heavy,FPP,
            centric_singletons=single,**OPTIONS)
