import numpy as np
from oasis_py.phasing import PhaseResult
from oasis_py.ensemble import mix_phases


def result(phase,fom):return PhaseResult(np.array(phase,dtype=float),np.array(fom,dtype=float),None,None,{'method':'test'})


def test_duplicate_posterior_is_not_counted_twice():
    a=result([.2,1.3],[.4,.8]);b=mix_phases(a,a)
    np.testing.assert_allclose(b.fom,a.fom,atol=1e-15)
    np.testing.assert_allclose(np.exp(1j*b.phase),np.exp(1j*a.phase),atol=1e-15)
    assert b.hl is None and b.probability is None


def test_disagreeing_equal_components_reduce_confidence():
    a=result([0],[.8]);b=result([np.pi],[.8]);m=mix_phases(a,b)
    assert m.fom[0]<1e-15


def test_mixture_matches_exact_complex_first_moment():
    a=result([.2,1.3],[.4,.8]);b=result([2.7,-.6],[.5,.9]);m=mix_phases(a,b,.3)
    expected=.3*a.fom*np.exp(1j*a.phase)+.7*b.fom*np.exp(1j*b.phase)
    np.testing.assert_allclose(m.fom*np.exp(1j*m.phase),expected,atol=1e-15)


def test_native_density_modifies_each_hypothesis_before_mixing(monkeypatch):
    from oasis_py import phasing,sad,native_density
    from oasis_py.ensemble import phase_sad_ensemble
    a=result([0],[.5]);a.metadata={'method':'pplus'}
    b=result([1],[.5]);b.metadata={'method':'sad_likelihood','anomalous_scale_fit_version':'wide_grid_null_v2'}
    monkeypatch.setattr(phasing,'phase_native',lambda *args,**kwargs:a)
    monkeypatch.setattr(sad,'phase_sad_likelihood',lambda *args,**kwargs:b)
    seen=[]
    def modify(data,initial,**kwargs):
        seen.append(initial.metadata['method'])
        return result([float(initial.phase[0]**2+1)],[.7])
    monkeypatch.setattr(native_density,'density_modify_native',modify)
    mixed=phase_sad_ensemble(None,None,None,None,density_options={})
    expected=.35*(np.exp(1j*1)+np.exp(1j*2))
    np.testing.assert_allclose(mixed.fom*np.exp(1j*mixed.phase),[expected],atol=1e-15)
    assert seen==['pplus','sad_likelihood']
    assert mixed.metadata['ensemble_stage']=='after_separate_native_dm'


def test_mtz_mixture_obeys_circular_wrap_and_omits_invented_HL(tmp_path):
    from test_science import dataset
    from oasis_py.data import write_mtz
    from oasis_py.ensemble import mix_mtz_files
    import gemmi
    d=dataset('P 1');n=len(d.f)
    write_mtz(tmp_path/'a.mtz',d,np.full(n,np.deg2rad(179)),np.full(n,.8))
    write_mtz(tmp_path/'b.mtz',d,np.full(n,np.deg2rad(-179)),np.full(n,.8))
    mix_mtz_files(tmp_path/'a.mtz',tmp_path/'b.mtz',tmp_path/'mixed.mtz')
    m=gemmi.read_mtz_file(str(tmp_path/'mixed.mtz'))
    np.testing.assert_allclose(m.column_with_label('PHIB').array,180,atol=1e-5)
    assert m.column_with_label('HLA') is None and m.column_with_label('PPLUS') is None


def test_mtz_mixture_refuses_different_observations(tmp_path):
    from test_science import dataset
    from oasis_py.data import write_mtz
    from oasis_py.ensemble import mix_mtz_files
    import pytest
    d=dataset('P 1');n=len(d.f)
    write_mtz(tmp_path/'a.mtz',d,np.zeros(n),np.full(n,.8))
    d.f*=2;write_mtz(tmp_path/'b.mtz',d,np.zeros(n),np.full(n,.8))
    with pytest.raises(ValueError,match='same observed amplitudes'):mix_mtz_files(tmp_path/'a.mtz',tmp_path/'b.mtz',tmp_path/'mixed.mtz')
