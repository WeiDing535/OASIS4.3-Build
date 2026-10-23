import numpy as np
from oasis_py.clipper_likelihood import reflection_classes,hl_target
from oasis_py.hl_constraints import canonicalize_centric_hl
from oasis_py.sad import posterior_moments
from test_science import dataset


def test_canonical_centric_hl_preserves_prior_and_combined_posteriors():
    data=dataset('P 21 21 21',extent=4);n=len(data.f)
    rng=np.random.default_rng(20261017);hl=rng.normal(size=(n,4))*30
    eps,cen,theta=reflection_classes(data)
    # A large forbidden component models the exact observed failure condition.
    hl[cen,0]+=-np.sin(theta[cen])*1400
    hl[cen,1]+=np.cos(theta[cen])*1400
    canonical,meta=canonicalize_centric_hl(data,hl)
    np.testing.assert_array_equal(canonical[~cen],hl[~cen])
    assert meta['phase_likelihood_ratios_preserved'] and not meta['target_reference_used']
    for extra in [np.zeros((n,4)),rng.normal(size=(n,4))*2]:
        p,f=posterior_moments(data,hl+extra)
        q,g=posterior_moments(data,canonical+extra)
        np.testing.assert_allclose(f*np.exp(1j*p),g*np.exp(1j*q),atol=1e-11)


def test_canonicalization_is_idempotent_and_does_not_mutate_inputs():
    data=dataset('P 21 21 21');hl=np.random.default_rng(6).normal(size=(len(data.f),4))
    before=hl.copy();a,_=canonicalize_centric_hl(data,hl);b,_=canonicalize_centric_hl(data,a)
    np.testing.assert_array_equal(before,hl);np.testing.assert_allclose(a,b,atol=1e-14)


def test_large_forbidden_component_is_harmless_to_stable_kernel():
    data=dataset('P 21 21 21');eps,cen,allowed=reflection_classes(data);n=len(data.f)
    hl=np.zeros((n,4));hl[cen,0]=-1500*np.sin(allowed[cen]);hl[cen,1]=1500*np.cos(allowed[cen])
    a,_=canonicalize_centric_hl(data,hl)
    x=hl_target(data.f,data.sigf,np.ones(n,complex),hl,1,1,eps,cen,allowed,samples=72)
    y=hl_target(data.f,data.sigf,np.ones(n,complex),a,1,1,eps,cen,allowed,samples=72)
    np.testing.assert_allclose(x.fom*np.exp(1j*x.phase),y.fom*np.exp(1j*y.phase),atol=1e-11)
