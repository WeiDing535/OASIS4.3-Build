import numpy as np
import gemmi,pytest
from test_science import dataset
from oasis_py.anisotropy import fit_anisotropy,corrected_reflections


def test_known_anisotropic_wilson_tensor_recovered_and_symmetry_enforced():
    d=dataset('P 21 21 21',extent=20);rng=np.random.default_rng(231)
    q=d.hkl@np.array(d.cell.frac.mat);truth=np.diag([18.,-5.,-13.])
    mu=10000*np.exp(-.5*np.einsum('ni,ij,nj->n',q,truth+np.eye(3)*22,q))
    d.f=np.sqrt(d.epsilon*mu*rng.exponential(size=len(q)));d.sigf=d.f/20
    factor,meta=fit_anisotropy(d)
    assert meta['status']=='fitted' and meta['independent_parameters']==2
    tensor=np.array(meta['tensor_cartesian_B'])
    np.testing.assert_allclose(tensor,truth,atol=3.5)
    assert abs(np.trace(tensor))<1e-10
    assert np.allclose(tensor,np.diag(np.diag(tensor)))
    assert factor.min()>=1/3-1e-10 and factor.max()<=3+1e-10


def test_cubic_symmetry_cannot_have_directional_B_tensor():
    d=dataset('P 2 3',extent=5);d.cell=gemmi.UnitCell(24,24,24,90,90,90);factor,meta=fit_anisotropy(d)
    assert meta['status']=='isotropic_by_crystal_symmetry'
    np.testing.assert_array_equal(factor,1)


def test_all_amplitudes_and_sigmas_transform_consistently_without_mutation():
    d=dataset();d.dano=np.ones(len(d.f))*2;d.sigdano=np.ones(len(d.f))*.3
    factors=np.linspace(.5,2,len(d.f));other=corrected_reflections(d,factors)
    np.testing.assert_allclose(other.f/other.sigf,d.f/d.sigf)
    np.testing.assert_allclose(other.dano/other.sigdano,d.dano/d.sigdano)
    np.testing.assert_array_equal(d.f,100)
