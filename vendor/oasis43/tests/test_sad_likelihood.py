import numpy as np
import gemmi,pytest
from oasis_py.data import Reflections
from oasis_py.sad import fit_anomalous_scale,posterior_moments,phase_sad_likelihood
from oasis_py.phasing import sad_hl


def simple(n=5000,sg='P 1'):
    cell=gemmi.UnitCell(30,32,34,90,90,90)
    return Reflections(np.column_stack([np.arange(1,n+1),np.ones(n),np.ones(n)]),np.full(n,100.),np.ones(n),cell,gemmi.SpaceGroup(sg),dano=np.zeros(n),sigdano=np.ones(n)*.2)


def test_anomalous_scale_recovers_known_scale_without_protein_phases():
    rng=np.random.default_rng(876);d=simple();a=rng.uniform(.5,2,len(d.f));scale=2.3
    d.dano=2*scale*a*np.cos(rng.uniform(-np.pi,np.pi,len(d.f)))+rng.normal(0,.2,len(d.f))
    fitted,meta=fit_anomalous_scale(d,1j*a)
    assert abs(fitted/scale-1)<.025 and meta['scale_fit_uses_reference_phases'] is False
    d.dano*=13;d.sigdano*=13
    rescaled,_=fit_anomalous_scale(d,1j*a)
    assert abs(rescaled/fitted-13)<1e-5


def test_sad_posterior_retains_both_modes_instead_of_clipping_noisy_observations():
    d=simple(1);hl=sad_hl(np.array([.3]),np.array([2.]),np.array([0.]),np.array([.2]))
    phase,fom=posterior_moments(d,hl)
    assert fom[0]<1e-12
    # A noisy observation outside the nominal doublet remains integrable.
    hl=sad_hl(np.array([.3]),np.array([2.]),np.array([9.]),np.array([3.]))
    phase,fom=posterior_moments(d,hl)
    assert np.isfinite(phase[0]) and 0<fom[0]<1


def test_centrosymmetric_reflection_moments_stay_on_allowed_symmetry_phases():
    cell=gemmi.UnitCell(30,32,34,90,90,90);sg=gemmi.SpaceGroup('P 21 21 21')
    hkl=gemmi.make_miller_array(cell,sg,4.,unique=True)
    d=Reflections(hkl,np.ones(len(hkl))*100,np.ones(len(hkl)),cell,sg)
    rng=np.random.default_rng(54);hl=rng.normal(size=(len(hkl),4))
    phase,fom=posterior_moments(d,hl)
    for i in np.flatnonzero(d.centric):
        h=d.hkl[i].tolist();op=next(op for op in sg.operations() if op.apply_to_hkl(h)==[-v for v in h])
        assert abs(np.sin(phase[i]+op.phase_shift(h)/2))<1e-10
    assert np.all((fom>=0)&(fom<=1))


def test_posterior_rejects_nonfinite_coefficients():
    d=simple(1)
    with pytest.raises(ValueError,match='Invalid SAD'):posterior_moments(d,[[np.nan,0,0,0]])


def test_likelihood_output_hl_and_first_moments_agree():
    from test_science import dataset
    from oasis_py.data import Atoms,structure_factors
    d=dataset('P 21 21 21',extent=4);rng=np.random.default_rng(7)
    full=Atoms(rng.random((22,3)),['C']*20+['Se']*2,np.ones(22),np.full(22,15.))
    heavy=full.subset([20,21]);p=structure_factors(d,full);q=structure_factors(d,heavy,{'Se':3.8})
    d.f=np.maximum(abs(p),.01);d.sigf=np.maximum(d.f*.02,.1);d.dano=abs(p+q)-abs(np.conj(p)-np.conj(q));d.sigdano=np.ones(len(d.f))*.5
    result=phase_sad_likelihood(d,{'C':20,'Se':2},heavy,{'Se':3.8},cycles=3,max_neighbors=120,max_per_target=32)
    phase,fom=posterior_moments(d,result.hl)
    np.testing.assert_allclose(np.exp(1j*phase),np.exp(1j*result.phase),atol=1e-13)
    np.testing.assert_allclose(fom,result.fom,atol=1e-13)
    assert result.probability is None


def test_scale_search_can_find_signal_when_unweighted_excess_variance_is_negative():
    rng=np.random.default_rng(901);d=simple(12000)
    a=np.r_[np.full(3000,2.),np.full(9000,.05)]
    d.sigdano=np.r_[np.full(3000,.2),np.full(9000,20.)]
    # Most noisy observations have conservative sigmas; global moment is
    # negative, but the strong-site likelihood identifies a positive scale.
    d.dano=2*1.7*a*np.cos(rng.uniform(-np.pi,np.pi,len(a)))+rng.normal(size=len(a))*np.r_[np.full(3000,.2),np.full(9000,2.)]
    assert np.sum(d.dano**2-d.sigdano**2)<0
    alpha,meta=fit_anomalous_scale(d,1j*a)
    assert alpha==pytest.approx(1.7,rel=.04)
    assert meta['loglikelihood_gain_vs_null']>100


def test_zero_anomalous_difference_can_select_exact_null_scale():
    d=simple(100);alpha,meta=fit_anomalous_scale(d,np.ones(100,complex))
    assert alpha==0 and meta['loglikelihood_gain_vs_null']==0
