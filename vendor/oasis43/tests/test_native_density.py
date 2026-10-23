import gemmi
import numpy as np
import pytest
from oasis_py.native_density import metric_filter,smooth,protein_mask,histogram_transform,apply_transform,fit_shell_likelihood
from oasis_py.sad import posterior_moments
from oasis_py.substructure import coefficient_map
from test_science import dataset


@pytest.mark.parametrize('sg',['P 1','P 21 21 21','P 31 2 1'])
def test_fft_roundtrip_preserves_complex_structure_factors(sg):
    from oasis_py.data import Atoms,structure_factors
    d=dataset(sg,extent=5)
    atoms=Atoms(np.array([[.17,.26,.38],[.41,.52,.73]]),['C','O'],np.ones(2),np.ones(2)*10)
    sf=structure_factors(d,atoms)
    grid=coefficient_map(d,sf)
    actual=gemmi.transform_map_to_f_phi(grid).get_value_by_hkl(d.hkl)
    np.testing.assert_allclose(actual,sf,atol=2e-5,rtol=2e-5)


def test_skew_cell_filter_uses_cartesian_metric_and_preserves_constant():
    cell=gemmi.UnitCell(24,24,28,90,90,120);shape=(24,24,28)
    kernel=metric_filter(cell,shape,2.)
    for h in [[1,0,0],[1,1,0],[2,-1,1]]:
        expected=np.exp(-2*np.pi**2*4*cell.calculate_1_d2(h))
        assert kernel[tuple(h)]==pytest.approx(expected)
    np.testing.assert_allclose(smooth(np.ones(shape)*3,kernel),3,atol=1e-12)


def test_mask_solvent_volume_and_screw_symmetry():
    d=dataset('P 21 21 21',extent=5)
    rng=np.random.default_rng(13);grid=coefficient_map(d,rng.normal(size=len(d.f))+1j*rng.normal(size=len(d.f)))
    mask=protein_mask(grid,.63,metric_filter(d.cell,grid.shape,4.))
    assert mask.mean()==pytest.approx(.37,abs=1e-7)
    grid.array[:]=mask;before=np.array(grid.array);grid.symmetrize_avg()
    np.testing.assert_allclose(before,grid.array,atol=2e-6)


def test_histogram_mapping_is_monotonic_and_flattens_solvent():
    rho=np.linspace(-3,5,1000).reshape(10,10,10)
    mask=np.zeros_like(rho);mask[5:]=1
    transformed,mapping=histogram_transform(rho,mask,np.linspace(-1,2,513)**3)
    assert np.std(transformed[mask==0])<1e-12
    assert np.all(np.diff(transformed[mask==1])>=0)
    np.testing.assert_allclose(apply_transform(rho,mask,mapping),transformed)


def test_shell_mlhl_is_scale_invariant_and_does_not_mutate_prior():
    d=dataset('P 1',extent=5);rng=np.random.default_rng(76)
    d.f=rng.uniform(3,10,len(d.f));d.sigf=d.f*.03
    phase=rng.uniform(0,2*np.pi,len(d.f));fc=.7*d.f*np.exp(1j*phase)
    prior=np.column_stack([np.cos(phase),np.sin(phase),np.zeros((len(d.f),2))]);original=prior.copy()
    a,_=fit_shell_likelihood(d,fc,prior)
    d.f*=31;d.sigf*=31
    b,_=fit_shell_likelihood(d,fc*31,prior)
    np.testing.assert_allclose(a,b,atol=1e-5)
    np.testing.assert_array_equal(prior,original)
    assert np.all(np.isfinite(a)) and np.max(np.abs(a))<=20.00001


def test_legacy_pplus_bridge_preserves_centroid_but_does_not_use_observation_hl():
    from oasis_py.native_density import initial_prior
    from oasis_py.phasing import PhaseResult
    d=dataset('P 21 21 21',extent=4)
    rng=np.random.default_rng(63);hl=rng.normal(size=(len(d.f),4))
    phase,fom=posterior_moments(d,hl);fom=np.minimum(fom,.8)
    initial=PhaseResult(phase,fom,None,np.zeros_like(hl),{'method':'pplus'})
    prior,kind,_=initial_prior(d,initial)
    p,f=posterior_moments(d,prior)
    np.testing.assert_allclose(f*np.exp(1j*p),fom*np.exp(1j*phase),atol=1e-9)
    assert 'branch_shape_not_preserved' in kind


def test_native_pipeline_never_executes_external_program(tmp_path,monkeypatch):
    import subprocess
    from oasis_py.native_density import density_modify_native
    from oasis_py.phasing import PhaseResult
    d=dataset('P 21 21 21',extent=4);rng=np.random.default_rng(9)
    d.f=rng.uniform(30,90,len(d.f));hl=rng.normal(size=(len(d.f),4))*.3
    phase,fom=posterior_moments(d,hl);initial=PhaseResult(phase,fom,None,hl,{'method':'sad_likelihood'})
    ref=tmp_path/'reference.mtz'
    m=gemmi.Mtz(with_base=True);m.cell=d.cell;m.spacegroup=d.spacegroup;m.add_dataset('reference')
    m.add_column('FC_BEST.F_phi.F','F');m.add_column('FC_BEST.F_phi.phi','P')
    m.set_data(np.column_stack([d.hkl,d.f,np.rad2deg(phase)]).astype(np.float32));m.write_to_file(str(ref))
    def fail(*a,**k):raise AssertionError('External program execution forbidden')
    monkeypatch.setattr(subprocess,'run',fail)
    result=density_modify_native(d,initial,.5,ref,cycles=2)
    p,f=posterior_moments(d,result.hl)
    np.testing.assert_allclose(f*np.exp(1j*p),result.fom*np.exp(1j*result.phase),atol=1e-12)
    assert np.isfinite(result.hl).all() and result.metadata['prior_repeatedly_accumulated'] is False
    np.testing.assert_array_equal(hl,initial.hl)


def test_grid_resource_fallback_is_symmetry_compatible_and_oversampled():
    from oasis_py.native_density import working_grid_size
    class Asu:
        def get_size_for_hkl(self,sample_rate):return {3.:[180,320,240],2.5:[144,256,192],2.1:[120,216,160]}[sample_rate]
    shape,rate=working_grid_size(Asu(),8_000_000)
    assert shape==[144,256,192] and rate==2.5
    with pytest.raises(ValueError,match='memory limit'):working_grid_size(Asu(),100)
