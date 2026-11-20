"""Independent circular-moment checks and unchanged legacy-control contracts."""
from types import SimpleNamespace
import numpy as np
import gemmi
import pytest
from scipy.special import i0e,i1e,ive,logsumexp
from phasebuild.pplus_experimental import (rank_cosine,kernel_moments,
    maximum_entropy_hl,bimodal_prior,phase_pplus_trial,mtz_safe_bimodal_prior)
from test_sad_pipeline import inputs


def test_rank_fit_is_uniform_in_cosine_preserves_ties_and_centric_angles():
    d=np.array([0.,np.pi/3,np.pi/2,np.pi,.2])
    result=rank_cosine(d,[False]*4+[True])
    np.testing.assert_allclose(np.cos(result[:4]),[.75,.25,-.25,-.75],atol=1e-14)
    assert result[-1]==d[-1]
    tied=rank_cosine(np.array([0.,0.,np.pi]),np.zeros(3,bool))
    assert tied[0]==tied[1]
    np.testing.assert_array_equal(d,[0.,np.pi/3,np.pi/2,np.pi,.2])


def test_kernel_moments_match_independent_high_resolution_quadrature():
    kappas=np.array([0.,.3,2.,10.,100.])
    first=i1e(kappas)/i0e(kappas)
    a,b=kernel_moments(first)
    theta=np.arange(32768)*2*np.pi/32768
    lp=kappas[:,None]*np.cos(theta)
    weights=np.exp(lp-logsumexp(lp,axis=1)[:,None])
    np.testing.assert_allclose(a,weights@np.cos(theta),atol=1e-12)
    np.testing.assert_allclose(b,weights@np.cos(2*theta),atol=1e-12)


def test_two_harmonic_projection_retains_ambiguous_antipodal_mass():
    data=SimpleNamespace(f=np.ones(1),centric=np.zeros(1,bool))
    prior=bimodal_prior(data,[0.],[np.pi/2],[.5],[.9])
    theta=np.arange(1440)*2*np.pi/1440
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    lp=prior.hl[0]@basis;w=np.exp(lp-logsumexp(lp))
    assert abs(w@np.exp(1j*theta))<1e-7
    assert (w@np.cos(2*theta))<-.6
    assert w[360]>w[0] and w[1080]>w[720]
    assert np.argmax(w[:720])==360 and np.argmax(w[720:])+720==1080
    assert prior.metadata['maximum_first_moment_difference']<1e-7


def test_projection_matches_both_moments_over_sharp_and_asymmetric_mixtures():
    rng=np.random.default_rng(7);n=18
    ref=rng.uniform(-np.pi,np.pi,n);d=rng.uniform(0,np.pi,n)
    p=rng.uniform(0,1,n);g=np.tile([0.,.1,.5,.85,.97,.9999],3)
    data=SimpleNamespace(f=np.ones(n),centric=np.zeros(n,bool))
    prior=bimodal_prior(data,ref,d,p,g)
    first,second=kernel_moments(g)
    expected1=first*np.exp(1j*ref)*(np.cos(d)+1j*(2*p-1)*np.sin(d))
    expected2=second*np.exp(2j*ref)*(np.cos(2*d)+1j*(2*p-1)*np.sin(2*d))
    theta=np.arange(720)*2*np.pi/720
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    lp=prior.hl@basis;w=np.exp(lp-logsumexp(lp,axis=1)[:,None])
    np.testing.assert_allclose(w@np.exp(1j*theta),expected1,atol=1e-7)
    np.testing.assert_allclose(w@np.exp(2j*theta),expected2,atol=1e-7)


def test_centric_projection_remains_on_allowed_two_state_phases():
    from oasis_py.data import Reflections
    data=Reflections([[1,0,1],[1,1,1]],[10,10],[1,1],
        gemmi.UnitCell(30,30,30,90,90,90),gemmi.SpaceGroup('P 21'))
    assert data.centric.tolist()==[True,False]
    result=bimodal_prior(data,[np.pi/2,0],[np.pi/2,np.pi/3],[.7,.5],[.8,.6])
    assert abs(np.sin(result.phase[0]))<1e-12
    assert result.fom[0]==pytest.approx(.32,abs=1e-10)


def test_legacy_control_and_original_observation_prior_are_unchanged(inputs):
    from oasis_py.data import read_mtz,read_atoms
    from oasis_py.sequence import sequence_asu
    from oasis_py.phasing import phase_native
    data,_=read_mtz(inputs['sad'],'sad',free_label='FreeR_flag',free_value=7)
    atoms=read_atoms(inputs['heavy'],data)
    composition=sequence_asu(inputs['sequence'],2,data.cell,data.spacegroup)['asu_composition']
    options=dict(cycles=2,max_neighbors=40,max_per_target=12)
    baseline=phase_native(data,'sad',composition,heavy=atoms,fpp={'S':.56},**options)
    control=phase_pplus_trial(data,composition,atoms,{'S':.56},**options)
    for name in ['phase','fom','probability','hl']:
        np.testing.assert_allclose(getattr(control.raw,name),getattr(baseline,name),rtol=0,atol=1e-14)
    ranked=phase_pplus_trial(data,composition,atoms,{'S':.56},cosine_rank=True,**options)
    np.testing.assert_array_equal(control.raw.hl,ranked.raw.hl)
    assert not np.allclose(control.raw.fom,ranked.raw.fom)


def test_invalid_or_infeasible_moment_constraints_raise():
    with pytest.raises(ValueError):kernel_moments([np.nan])
    with pytest.raises(ValueError):rank_cosine(np.array([-1.]),[False])
    with pytest.raises(RuntimeError,match='did not converge'):
        maximum_entropy_hl([[2.,0.,0.,0.]],max_iterations=5)


def test_close_sharp_doublets_converge_and_survive_actual_float32_mtz(tmp_path):
    from oasis_py.data import Reflections,write_mtz
    from oasis_py.sad import posterior_moments
    # Previously stalled at the fixed coefficient-step cap, then exceeded
    # float32 moment tolerance even when the unquantized solve converged.
    ref=np.array([-.78767774448362,-.4441303836241395,-.3780789263732425,.6261796931789809])
    d=np.array([.052002278118614076,.056888289803994026,.09750478487491918,.13964582691575075])
    p=np.array([.7,.3,.3,.5]);g=np.full(4,.9999)
    data=Reflections([[1,1,1],[2,1,1],[3,1,1],[4,1,1]],np.ones(4),np.ones(4),
        gemmi.UnitCell(30,30,30,90,90,90),gemmi.SpaceGroup('P 1'))
    original=bimodal_prior(data,ref,d,p,g)
    safe=mtz_safe_bimodal_prior(data,original)
    path=tmp_path/'safe.mtz';write_mtz(path,data,safe.phase,safe.fom,safe.hl)
    m=gemmi.read_mtz_file(str(path));a=np.array(m);labels=m.column_labels()
    stored=a[:,[labels.index(x) for x in ['HLA','HLB','HLC','HLD']]]
    phase,fom=posterior_moments(data,stored.astype(float))
    np.testing.assert_allclose(fom*np.exp(1j*phase),original.fom*np.exp(1j*original.phase),atol=1e-5)
    z=a[:,labels.index('FOM')]*np.exp(1j*np.deg2rad(a[:,labels.index('PHIB')]))
    np.testing.assert_allclose(fom*np.exp(1j*phase),z,atol=1e-6)
    assert safe.metadata['mtz_hl_quantization']['rows_with_neighbor_rounding']>0
    np.testing.assert_array_equal(original.probability,p)


def test_reader_keeps_centric_singletons_and_excludes_missing_acentric_pairs(tmp_path):
    from phasebuild.sad_input import read_sad_pairs
    m=gemmi.Mtz(with_base=True);m.spacegroup=gemmi.SpaceGroup('P 21');m.set_cell_for_all(gemmi.UnitCell(30,30,30,90,90,90));m.add_dataset('test')
    for name,kind in [('F(+)','G'),('SIGF(+)','L'),('F(-)','G'),('SIGF(-)','L'),('FreeR_flag','I')]:m.add_column(name,kind)
    m.set_data(np.array([[1,0,1,10,1,np.nan,np.nan,1],[1,1,1,12,1,10,1,1],
        [2,1,1,12,1,np.nan,np.nan,1],[2,0,1,9,1,np.nan,np.nan,0]],np.float32))
    path=tmp_path/'centric.mtz';m.write_to_file(str(path));before=path.read_bytes()
    data,single,meta=read_sad_pairs(path)
    assert len(data.f)==2 and single.sum()==1 and np.all(data.centric[single])
    assert meta['noncentric_missing_pairs_excluded']==1
    assert data.f[single][0]==10 and data.sigf[single][0]==1
    assert path.read_bytes()==before


def test_singleton_anomalous_likelihood_and_sigma_placeholder_are_not_used(inputs):
    from oasis_py.data import Reflections,Atoms
    data=Reflections([[1,0,1],[1,1,1],[2,1,1]],[10,12,14],[1,1,1],
        gemmi.UnitCell(30,30,30,90,90,90),gemmi.SpaceGroup('P 21'),dano=[0,1,2],sigdano=[1,1,1])
    atoms=Atoms([[.13,.17,.19]],['Se'],[1],[10]);single=data.centric
    args=dict(cycles=2,max_neighbors=20,max_per_target=8,centric_singletons=single)
    a=phase_pplus_trial(data,{'C':200,'N':50,'O':60,'Se':1},atoms,{'Se':6.},**args)
    data.sigdano[single]=1e12
    b=phase_pplus_trial(data,{'C':200,'N':50,'O':60,'Se':1},atoms,{'Se':6.},**args)
    np.testing.assert_array_equal(a.raw.hl[single],0)
    np.testing.assert_array_equal(a.raw.phase,b.raw.phase)
    np.testing.assert_array_equal(a.raw.fom,b.raw.fom)
    assert a.raw.metadata['calibrated_error_mean']==pytest.approx(.5)


def test_shared_normalization_is_reused_and_rejects_changed_input(inputs,monkeypatch):
    import phasebuild.pplus_experimental as module
    from oasis_py.data import read_mtz,read_atoms
    from oasis_py.sequence import sequence_asu
    data,_=read_mtz(inputs['sad'],'sad',free_label='FreeR_flag',free_value=7)
    atoms=read_atoms(inputs['heavy'],data)
    comp=sequence_asu(inputs['sequence'],2,data.cell,data.spacegroup)['asu_composition']
    context=module.prepare_pplus_normalization(data,comp)
    calibration=module.prepare_pplus_calibration(data,comp,atoms,{'S':.56},
        normalization_context=context,max_neighbors=40,max_per_target=12)
    def forbidden(*args,**kwargs):raise AssertionError('Unexpected repeated fit')
    monkeypatch.setattr(module,'normalize',forbidden)
    monkeypatch.setattr(module,'structure_factors',forbidden)
    options=dict(cycles=2,max_neighbors=40,max_per_target=12,calibration_context=calibration)
    a=module.phase_pplus_trial(data,comp,atoms,{'S':.56},**options)
    b=module.phase_pplus_trial(data,comp,atoms,{'S':.56},cosine_rank=True,**options)
    np.testing.assert_array_equal(a.raw.hl,b.raw.hl)
    assert a.raw.metadata['normalization_context_sha256']==b.raw.metadata['normalization_context_sha256']
    data.f[0]*=1.01
    with pytest.raises(ValueError,match='does not match'):
        module.phase_pplus_trial(data,comp,atoms,{'S':.56},**options)


def test_persisted_prior_is_aligned_by_hkl_not_mtz_row_order(tmp_path):
    from oasis_py.data import Reflections,write_mtz
    from phasebuild.sad_input import read_phase_prior
    data=Reflections([[3,1,1],[1,1,1],[2,1,1]],[30,10,20],[1,1,1],
        gemmi.UnitCell(30,30,30,90,90,90),gemmi.SpaceGroup('P 1'))
    phase=np.array([.3,.1,.2]);fom=np.array([.6,.2,.4]);hl=np.arange(12).reshape(3,4)/10
    p=np.array([.7,.8,.9]);path=tmp_path/'prior.mtz'
    write_mtz(path,data,phase,fom,hl,p)
    # Writer sorts MTZ rows; the observations deliberately retain another order.
    assert not np.array_equal(gemmi.read_mtz_file(str(path)).make_miller_array(),data.hkl)
    result=read_phase_prior(path,data)
    np.testing.assert_allclose(result.phase,phase,atol=1e-7)
    np.testing.assert_allclose(result.fom,fom,atol=1e-7)
    np.testing.assert_allclose(result.hl,hl,atol=1e-7)
    np.testing.assert_allclose(result.probability,p,atol=1e-7)
    data.hkl[0]=[4,1,1]
    with pytest.raises(ValueError,match='HKL sets differ'):read_phase_prior(path,data)
