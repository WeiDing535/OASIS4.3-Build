from pathlib import Path
import json
import numpy as np
import pytest
from oasis_py.clipper_likelihood import hl_target


def fixture():
    saved=json.loads((Path(__file__).parent/'fixtures/clipper_mlhl.json').read_text())
    d={k:np.asarray(v) for k,v in saved['inputs'].items()}
    args=[d['fo'],d['sig'],d['fc_real']+1j*d['fc_imag'],d['hl'],d['s'],d['w'],
          d['epsilon'],d['centric'],d['allowed']]
    return saved,args


@pytest.mark.parametrize('samples',[24,72,720])
def test_matches_unmodified_cpp_clipper_oracle(samples):
    saved,args=fixture();v=hl_target(*args,samples=samples)
    actual=np.column_stack([v.value,v.ds,v.dw,v.dss,v.dww])
    expected=np.asarray(saved['reference'][str(samples)])
    assert np.max(abs(actual-expected)/(1+abs(expected)))<1e-9
    posterior=np.asarray(saved['posterior'][str(samples)])
    np.testing.assert_allclose(v.fom*np.exp(1j*v.phase),posterior[:,1]*np.exp(1j*posterior[:,0]),atol=1e-11)
    np.testing.assert_allclose(v.hl,posterior[:,2:],rtol=1e-13,atol=1e-12)


def test_scale_and_variance_derivatives_match_independent_finite_differences():
    _,args=fixture();v=hl_target(*args,samples=72)
    for pos,grad,curvature in [(4,v.ds,v.dss),(5,v.dw,v.dww)]:
        step=1e-5*np.maximum(args[pos],.01)
        plus=list(args);minus=list(args)
        plus[pos]=args[pos]+step;minus[pos]=args[pos]-step
        vp=hl_target(*plus,samples=72);vm=hl_target(*minus,samples=72)
        numerical=(vp.value-vm.value)/(2*step)
        np.testing.assert_allclose(grad,numerical,rtol=2e-5,atol=2e-6)
        gplus=vp.ds if pos==4 else vp.dw;gminus=vm.ds if pos==4 else vm.dw
        np.testing.assert_allclose(curvature,(gplus-gminus)/(2*step),rtol=2e-5,atol=2e-6)


def test_centric_restrictions_and_stable_large_second_harmonics():
    _,args=fixture();args[3]=args[3].copy();args[3][:,2:]=[10000,-10000]
    v=hl_target(*args,samples=72)
    assert all(np.all(np.isfinite(x)) for x in [v.value,v.ds,v.dw,v.phase,v.fom])
    np.testing.assert_allclose(np.sin(v.phase[args[7]]-args[8][args[7]]),0,atol=1e-12)
    assert np.all((v.fom>=0)&(v.fom<=1))


def test_zero_calculated_amplitude_preserves_prior_moments():
    from oasis_py.sad import posterior_moments
    from test_science import dataset
    d=dataset('P 21 21 21')
    from oasis_py.clipper_likelihood import reflection_classes
    eps,cen,allowed=reflection_classes(d);n=len(d.f)
    prior=np.random.default_rng(7).normal(size=(n,4))
    v=hl_target(d.f,d.sigf,np.zeros(n,complex),prior,1.,1.,eps,cen,allowed,samples=720)
    phase,fom=posterior_moments(d,prior)
    np.testing.assert_allclose(v.fom*np.exp(1j*v.phase),fom*np.exp(1j*phase),atol=1e-13)


@pytest.mark.parametrize('samples',[0,25,73])
def test_rejects_angular_grids_incompatible_with_centric_symmetry(samples):
    _,args=fixture()
    with pytest.raises(ValueError,match='angular grid'):hl_target(*args,samples=samples)


def test_rejects_nonfinite_inputs_and_zero_sigmas():
    _,args=fixture();args[1]=np.zeros_like(args[1])
    with pytest.raises(ValueError,match='Invalid likelihood'):hl_target(*args)
