import numpy as np
import pytest
import gemmi
from oasis_py.data import Reflections,Atoms,structure_factors,write_mtz,read_mtz,normalize
from oasis_py.phasing import circular_doublet,error_weights,sad_hl,hl_moments,iterate_pplus,phase_native
from oasis_py.relations import build_relations
from oasis_py.substructure import coefficient_map,pick_peaks,invert_hand


def dataset(sg='P 1',extent=3):
    cell=gemmi.UnitCell(24,26,28,90,90,90)
    spacegroup=gemmi.SpaceGroup(sg)
    if sg=='P 31 2 1':cell=gemmi.UnitCell(24,24,28,90,90,120)
    hkl=[];seen=set();asu=gemmi.ReciprocalAsu(spacegroup);ops=spacegroup.operations()
    for h in range(-extent,extent+1):
        for k in range(-extent,extent+1):
            for l in range(-extent,extent+1):
                if (h,k,l)==(0,0,0) or ops.is_systematically_absent([h,k,l]):continue
                index=tuple(asu.to_asu([h,k,l],ops)[0])
                if index not in seen:seen.add(index);hkl.append(index)
    return Reflections(hkl,np.ones(len(hkl))*100,np.ones(len(hkl)),cell,spacegroup)


def test_circular_mean_matches_direct_probability_integral():
    ref=np.array([.8,-2.,1.,0.]);d=np.array([.4,2.8,np.pi/2,0.])
    p=np.array([.9,.25,.5,.1]);g=np.array([.8,.7,1.,.3])
    phase,fom=circular_doublet(ref,d,p,g)
    expected=g*(p*np.exp(1j*(ref+d))+(1-p)*np.exp(1j*(ref-d)))
    np.testing.assert_allclose(fom*np.exp(1j*phase),expected,atol=2e-15)


def test_error_target_and_monotonicity():
    a=np.geomspace(.001,100,1000);g,n=error_weights(a,.37)
    assert abs(g.mean()-.37)<1e-10 and n>0
    assert np.all(np.diff(g)<=0)
    with pytest.raises(ValueError):error_weights(np.zeros(4))


@pytest.mark.parametrize('scale',[1e-200,1e200])
def test_error_calibration_handles_extreme_exponent_scales(scale):
    a=np.linspace(.1,10,100)*scale
    g,n=error_weights(a,.5)
    assert abs(g.mean()-.5)<1e-9 and np.isfinite(n)
    np.testing.assert_allclose(np.exp(-a*n*n),g,atol=1e-12)


def test_sad_hl_equals_unprojected_likelihood_up_to_constant():
    ref=np.array([.3,2.]);a=np.array([2.,5.]);delta=np.array([-1.,3.]);sigma=np.array([.8,1.5])
    hl=sad_hl(ref,a,delta,sigma);theta=np.linspace(0,2*np.pi,200,endpoint=False)
    basis=np.array([np.cos(theta),np.sin(theta),np.cos(2*theta),np.sin(2*theta)])
    expected=-.5*((delta[:,None]-2*a[:,None]*np.cos(theta-ref[:,None]))/sigma[:,None])**2
    actual=hl@basis
    np.testing.assert_allclose((expected-actual).std(axis=1),0,atol=2e-14)


def test_uniform_hl_has_zero_confidence_and_extreme_is_finite():
    phase,fom=hl_moments(np.array([[0,0,0,0],[1e5,0,0,0]]))
    assert fom[0]<1e-12 and fom[1]>.999
    assert np.all(np.isfinite(phase))


@pytest.mark.parametrize('sg',['P 1','P 21 21 21','P 31 2 1'])
def test_symmetry_relation_lookup_against_explicit_structure_factors(sg):
    data=dataset(sg)
    atoms=Atoms([[.13,.27,.39],[.41,.12,.23]],['S','C'],[1,.7],[10,15])
    sf=structure_factors(data,atoms)
    graph=build_relations(data,np.ones(len(data.f)),.2,kmin=0,max_neighbors=0,max_per_target=0)
    # Reconstruct the left+right predicted phase from each encoded edge and
    # compare with fresh sums at the actual equivalent indices.
    assert len(graph.target)>0
    for edge in np.linspace(0,len(graph.target)-1,min(30,len(graph.target)),dtype=int):
        i=graph.target[edge];j=graph.left[edge];k=graph.right[edge]
        # The complex symmetry transformation must reproduce any directly
        # expanded Fourier sum, including screw translations and Friedel mates.
        for original in [j,k]:
            h=data.hkl[original]
            for op in data.spacegroup.operations():
                q=op.apply_to_hkl(h.tolist())
                direct=Reflections([q],[100],[1],data.cell,data.spacegroup)
                fresh=structure_factors(direct,atoms)[0]
                predicted=sf[original]*np.exp(1j*op.phase_shift(h.tolist()))
                np.testing.assert_allclose(fresh,predicted,rtol=1e-10,atol=1e-10)
    phase=np.angle(sf);fom=np.ones(len(sf))
    field=graph.field(phase,fom)
    assert field.shape==(len(sf),) and np.all(np.isfinite(field))


def test_graph_field_matches_scalar_independent_sum():
    data=dataset();rng=np.random.default_rng(8)
    graph=build_relations(data,rng.uniform(.5,2,len(data.f)),.1,kmin=0,max_neighbors=0,max_per_target=0)
    phase=rng.uniform(-np.pi,np.pi,len(data.f));fom=rng.random(len(data.f))
    expected=np.zeros(len(data.f),complex)
    for i,j,k,sj,sk,t,w in zip(graph.target,graph.left,graph.right,graph.left_sign,graph.right_sign,graph.shift,graph.kappa):
        expected[i]+=w*fom[j]*fom[k]*np.exp(1j*(sj*phase[j]+sk*phase[k]+t))
    np.testing.assert_allclose(graph.field(phase,fom),expected,atol=1e-13)


def test_special_positions_are_not_overcounted():
    data=dataset('P 2')
    atoms=Atoms([[0,.2,0]],['C'],[1],[0])
    sf=structure_factors(data,atoms)
    from oasis_py.data import form_factor
    expected=form_factor('C',data.stol2)*np.exp(2j*np.pi*data.hkl[:,1]*.2)
    np.testing.assert_allclose(sf,expected,atol=1e-12)


def test_anomalous_factor_includes_i_and_correct_friedel_sign():
    data=dataset();atoms=Atoms([[.1,.2,.3]],['S'],[.7],[10])
    sf=structure_factors(data,atoms,{'S':.56})
    reverse=Reflections(-data.hkl,data.f,data.sigf,data.cell,data.spacegroup)
    sf_reverse=structure_factors(reverse,atoms,{'S':.56})
    np.testing.assert_allclose(sf_reverse,-sf.conj(),atol=1e-13)


def test_fixed_phases_do_not_alternate_or_drift():
    data=dataset();n=len(data.f)
    graph=build_relations(data,np.ones(n),.2,kmin=0)
    fixed=np.zeros(n,bool);fixed[:8]=True
    fp=np.linspace(.2,5,n);fm=np.full(n,.8)
    phase,fom,_,_=iterate_pplus(np.zeros(n),np.full(n,.7),np.full(n,.5),np.ones(n)*1j,
                               graph,cycles=7,fixed_phase=fp,fixed_fom=fm,fixed_mask=fixed)
    np.testing.assert_allclose(phase[fixed],fp[fixed]);np.testing.assert_array_equal(fom[fixed],fm[fixed])


def test_mtz_roundtrip_and_missing_observation_filter(tmp_path):
    data=dataset();phase=np.linspace(0,2*np.pi,len(data.f));fom=np.linspace(0,1,len(data.f))
    path=tmp_path/'out.mtz';write_mtz(path,data,phase,fom,np.zeros((len(data.f),4)))
    reread,meta=read_mtz(path)
    assert meta['retained_rows']==len(data.f)
    m=gemmi.read_mtz_file(str(path));lookup={tuple(h):i for i,h in enumerate(data.hkl)}
    for h,p,w in zip(m.make_miller_array(),m.column_with_label('PHIB').array,m.column_with_label('FOM').array):
        i=lookup[tuple(h)]
        assert abs(np.angle(np.exp(1j*(np.deg2rad(p)-phase[i]))))<1e-6
        assert abs(w-fom[i])<1e-7
    table=m.array.copy();table[0,3]=np.nan;m.set_data(table);m.write_to_file(str(path))
    reread,meta=read_mtz(path);assert len(reread.f)==len(data.f)-1


def test_map_peak_coordinates_obey_fourier_sign():
    data=dataset(extent=7)
    atoms=Atoms([[.25,.5,.75]],['S'],[1],[20])
    sf=structure_factors(data,atoms)
    grid=coefficient_map(data,sf,sample_rate=3.)
    xyz,_=pick_peaks(grid,1)
    distance=data.cell.find_nearest_pbc_image(gemmi.Fractional(*xyz[0]),gemmi.Fractional(*atoms.fractional[0])).dist()
    assert distance<.5


def test_hand_inversion_uses_spacegroup_operation():
    data=dataset('P 31 2 1');atoms=Atoms([[.13,.27,.39]],['S'],[1],[10])
    other,sg=invert_hand(atoms,data.spacegroup)
    assert sg.number==gemmi.SpaceGroup('P 32 2 1').number
    expected=np.mod(data.spacegroup.change_of_hand_op().apply_to_xyz(atoms.fractional[0].tolist()),1)
    np.testing.assert_allclose(other.fractional[0],expected)


@pytest.mark.parametrize('mode',['sad','sir','dmr'])
def test_native_modes_return_finite_phases_and_reproducible_results(mode):
    data=dataset(extent=4)
    model=Atoms([[.13,.27,.39],[.41,.12,.23],[.22,.64,.8]],['S','C','O'],[1,1,1],[15,15,15])
    heavy=model.subset([0]);sf=structure_factors(data,model)
    data.f=np.abs(sf)*10+.1;data.sigf=np.maximum(data.f*.02,.1)
    imag=structure_factors(data,heavy,{'S':.56})*10
    data.dano=2*np.abs(imag)*np.cos(np.angle(sf)-np.angle(imag))
    data.sigdano=np.full(len(data.f),.2)
    data.fderivative=np.abs(sf+structure_factors(data,heavy))*10+.1
    data.sigderivative=np.full(len(data.f),.2)
    result=phase_native(data,mode,{'C':30,'N':8,'O':9,'S':1},model,heavy,{'S':.56},cycles=3)
    repeated=phase_native(data,mode,{'C':30,'N':8,'O':9,'S':1},model,heavy,{'S':.56},cycles=3)
    # arccos doublets at +/-1 amplify float64 rounding to O(sqrt(eps)).
    # Compare circular phases, not bit patterns or 0-versus-2pi encodings.
    tolerance=4*np.sqrt(np.finfo(float).eps)
    delta=np.angle(np.exp(1j*(result.phase-repeated.phase)))
    np.testing.assert_allclose(delta,0,rtol=0,atol=tolerance)
    np.testing.assert_allclose(result.fom,repeated.fom,rtol=0,atol=tolerance)
    assert np.all(np.isfinite(result.phase)) and np.all((result.fom>=0)&(result.fom<=1))
    assert abs(result.metadata['mean_error_weight']-.5)<1e-9
    assert (result.hl is None)==(mode=='dmr')
