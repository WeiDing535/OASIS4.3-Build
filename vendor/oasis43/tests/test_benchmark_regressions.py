"""Independent crystallographic checks and regressions discovered on real inputs."""
import numpy as np
import gemmi
import pytest
from oasis_py.data import Reflections,Atoms,structure_factors
from oasis_py.backends import prepare_phaser_pairs

@pytest.mark.parametrize('sg',['P 1','P 21','C 2','P 21 21 21','P 41','P 31','R 3:H','P 31 2 1','I 21 3'])
def test_numpy_structure_factors_against_independent_gemmi_calculator(sg):
    group=gemmi.SpaceGroup(sg);system=group.crystal_system_str()
    if system in ['trigonal','hexagonal']:cell=gemmi.UnitCell(24,24,30,90,90,120)
    elif system=='cubic':cell=gemmi.UnitCell(28,28,28,90,90,90)
    elif system=='tetragonal':cell=gemmi.UnitCell(24,24,30,90,90,90)
    elif system=='monoclinic':cell=gemmi.UnitCell(24,26,30,90,103,90)
    else:cell=gemmi.UnitCell(24,26,30,90,90,90)
    hkl=gemmi.make_miller_array(cell,group,4.,unique=True)
    data=Reflections(hkl,np.ones(len(hkl)),np.ones(len(hkl)),cell,group)
    rng=np.random.default_rng(42)
    atoms=Atoms(rng.random((6,3)),['C','N','O','S','Se','Pt'],[1,.7,.5,1,.8,.3],[10,12,14,16,18,20])
    st=gemmi.Structure();st.cell=cell;st.spacegroup_hm=group.hm
    model=gemmi.Model('1');chain=gemmi.Chain('A');res=gemmi.Residue();res.name='UNK';res.seqid=gemmi.SeqId(1,' ')
    for x,e,o,b in zip(atoms.fractional,atoms.elements,atoms.occupancy,atoms.b_iso):
        atom=gemmi.Atom();atom.name=e;atom.element=gemmi.Element(e);atom.pos=cell.orthogonalize(gemmi.Fractional(*x));atom.occ=o;atom.b_iso=b;res.add_atom(atom)
    chain.add_residue(res);model.add_chain(chain);st.add_model(model);st.setup_cell_images()
    calculator=gemmi.StructureFactorCalculatorX(st.cell)
    independent=np.array([calculator.calculate_sf_from_model(st[0],h.tolist()) for h in hkl])
    np.testing.assert_allclose(structure_factors(data,atoms),independent,rtol=3e-6,atol=3e-5)


def pair_mtz(path,values,types=('G','L','G','L')):
    mtz=gemmi.Mtz(with_base=True);mtz.spacegroup=gemmi.SpaceGroup('P 1');mtz.set_cell_for_all(gemmi.UnitCell(20,20,20,90,90,90));mtz.add_dataset('test')
    for label,kind in zip(['F+','SIGF+','F-','SIGF-'],types):mtz.add_column(label,kind)
    mtz.add_column('FreeR_flag','I')
    rows=np.array([[i+1,0,0,*v,i%2] for i,v in enumerate(values)],dtype=np.float32)
    mtz.set_data(rows);mtz.write_to_file(str(path))
    return ['F+','SIGF+','F-','SIGF-']


def test_phaser_filters_negative_missing_and_zero_sigma_pairs_without_losing_free_flags(tmp_path):
    original=tmp_path/'input.mtz';output=tmp_path/'filtered.mtz'
    labels=pair_mtz(original,[[10,1,11,1],[-1,1,2,1],[3,0,4,1],[3,1,np.nan,1],[12,2,13,2]])
    prepared,stats=prepare_phaser_pairs(original,labels,output)
    assert stats['retained_rows']==2 and stats['removed_rows']==3
    result=gemmi.read_mtz_file(prepared)
    np.testing.assert_array_equal(result.make_miller_array()[:,0],[1,5])
    np.testing.assert_array_equal(result.column_with_label('FreeR_flag').array,[0,0])
    assert gemmi.read_mtz_file(str(original)).nreflections==5


def test_phaser_reuses_valid_input(tmp_path):
    path=tmp_path/'input.mtz';labels=pair_mtz(path,[[10,1,11,1]])
    prepared,stats=prepare_phaser_pairs(path,labels,tmp_path/'unused.mtz')
    assert prepared==str(path.resolve()) and stats['removed_rows']==0
    assert not (tmp_path/'unused.mtz').exists()


@pytest.mark.parametrize('values,types,error',[
    ([[-1,1,2,1]],('G','L','G','L'),'No complete'),
    ([[10,1,11,1]],('J','Q','J','Q'),'Invalid amplitude')])
def test_phaser_rejects_nonphysical_or_intensity_inputs(tmp_path,values,types,error):
    path=tmp_path/'input.mtz';labels=pair_mtz(path,values,types)
    with pytest.raises(ValueError,match=error):prepare_phaser_pairs(path,labels,tmp_path/'out.mtz')
