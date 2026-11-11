"""Independent Gemmi SF oracle for heavy-only anisotropic ADP support."""
from pathlib import Path
import gemmi,numpy as np,pytest
from oasis_py.data import Reflections,read_atoms,structure_factors
from phasebuild.heavy_scattering import read_heavy,anomalous_factors,AnisotropicHeavy

def model(tmp_path,sg,anisotropic=True):
 cell=gemmi.UnitCell(20,20 if sg!='P 1' else 23,30,90,90,90)
 data=Reflections(np.array([[1,2,3],[2,0,1],[1,1,2],[3,1,0]]),np.full(4,10.),np.ones(4),cell,gemmi.SpaceGroup(sg))
 st=gemmi.Structure();st.cell=cell;st.spacegroup_hm=sg;m=gemmi.Model(1);chain=gemmi.Chain('A')
 for i,xyz in enumerate([[.13,.27,.39],[.31,.19,.23]]):
  res=gemmi.Residue();res.name='SE';res.seqid=gemmi.SeqId(i+1,' ');a=gemmi.Atom();a.name='SE';a.element=gemmi.Element('Se');a.pos=cell.orthogonalize(gemmi.Fractional(*xyz));a.occ=.7+i*.2;a.b_iso=12+i*3
  if anisotropic and i==0:a.aniso=gemmi.SMat33f(.12,.23,.34,.014,.025,-.031)
  res.add_atom(a);chain.add_residue(res)
 m.add_chain(chain);st.add_model(m);path=tmp_path/'heavy.pdb';st.write_pdb(str(path));return data,path

@pytest.mark.parametrize('sg',['P 1','P 4','I 4'])
def test_anisotropic_factors_against_independent_gemmi_normal_sf(tmp_path,sg):
 data,path=model(tmp_path,sg);heavy=read_heavy(path,data);assert isinstance(heavy,AnisotropicHeavy)
 small=gemmi.SmallStructure();small.cell=data.cell;small.spacegroup_hm=sg;small.spacegroup=gemmi.SpaceGroup(sg)
 for chain in heavy.structure[0]:
  for residue in chain:
   for a in residue:
    s=gemmi.SmallStructure.Site();s.label=str(len(small.sites));s.type_symbol='Se';s.element=a.element;s.fract=data.cell.fractionalize(a.pos);s.occ=a.occ;s.u_iso=a.b_iso/(8*np.pi**2)
    if a.aniso.nonzero():s.aniso=gemmi.SMat33d(a.aniso.u11,a.aniso.u22,a.aniso.u33,a.aniso.u12,a.aniso.u13,a.aniso.u23)
    small.add_site(s)
 small.setup_cell_images();calculator=gemmi.StructureFactorCalculatorX(small.cell)
 expected=np.array([calculator.calculate_sf_from_small_structure(small,list(h))/(gemmi.Element('Se').it92.calculate_sf(float(s2)))*6j for h,s2 in zip(data.hkl,data.stol2)])
 actual=anomalous_factors(data,heavy,{'Se':6});np.testing.assert_allclose(actual,expected,rtol=3e-7,atol=3e-7)

def test_isotropic_route_bitwise_unchanged(tmp_path):
 data,path=model(tmp_path,'I 4',False);heavy=read_heavy(path,data);assert not isinstance(heavy,AnisotropicHeavy)
 np.testing.assert_array_equal(anomalous_factors(data,heavy,{'Se':6}),structure_factors(data,read_atoms(path,data),{'Se':6}))

def test_nonpositive_anisotropic_tensor_refused(tmp_path):
 data,path=model(tmp_path,'P 1');st=gemmi.read_structure(str(path));st[0][0][0][0].aniso=gemmi.SMat33f(-.1,.2,.3,0,0,0);st.write_pdb(str(path))
 with pytest.raises(ValueError,match='tensor'):read_heavy(path,data)

@pytest.mark.parametrize('sg,cell_parameters',[
 ('P 1',(20,23,30,90,105,90)),
 ('P 21',(20,23,30,90,105,90)),
 ('P 31',(20,20,30,90,90,120)),
])
def test_cartesian_tensor_and_symmetry_in_nonorthogonal_cells(tmp_path,sg,cell_parameters):
 data,path=model(tmp_path,sg)
 data.cell=gemmi.UnitCell(*cell_parameters)
 st=gemmi.read_structure(str(path));st.cell=data.cell;st.spacegroup_hm=sg;st.write_pdb(str(path))
 heavy=read_heavy(path,data);heavy.structure.setup_cell_images();calculator=gemmi.StructureFactorCalculatorX(heavy.structure.cell)
 expected=np.array([calculator.calculate_sf_from_model(heavy.structure[0],list(h))/gemmi.Element('Se').it92.calculate_sf(float(s2))*6j for h,s2 in zip(data.hkl,data.stol2)])
 np.testing.assert_allclose(anomalous_factors(data,heavy,{'Se':6}),expected,rtol=1e-6,atol=1e-6)
