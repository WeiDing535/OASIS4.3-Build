import numpy as np
from phasebuild.seed_search import score_oriented_profiles
from phasebuild.data import Density
from phasebuild.profile import DensityProfile
from phasebuild.geometry import template

def test_chunked_context_search_matches_direct_three_residue_profile_scoring(tmp_path):
    import gemmi
    m=gemmi.Mtz(with_base=True);m.cell=gemmi.UnitCell(24,24,24,90,90,90);m.spacegroup=gemmi.SpaceGroup('P 1');m.add_dataset('test');m.add_column('FWT','F');m.add_column('PHWT','P')
    h=gemmi.make_miller_array(m.cell,m.spacegroup,3.);random=np.random.RandomState(3);m.set_data(np.column_stack([h,random.uniform(1,100,len(h)),random.uniform(0,360,len(h))]).astype('float32'));path=tmp_path/'map.mtz';m.write_to_file(str(path));d=Density(path)
    offsets=random.normal(size=(17,3));d.profile=DensityProfile(offsets,random.normal(size=17),np.ones(17),{})
    oriented=np.array([template(-63,-42),template(-120,130)]);centers=random.uniform(0,24,(5,3));xyz=centers[:,None,None,None,:]+oriented[None]
    values=d.profile.score(d,xyz);expected=values.mean(axis=-1)-.15*values.std(axis=-1)
    np.testing.assert_allclose(score_oriented_profiles(d,oriented,centers,block=2),expected,rtol=1e-12,atol=1e-12)
