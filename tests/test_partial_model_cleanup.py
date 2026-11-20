import gemmi,numpy as np,pytest
from phasebuild.partial_model_cleanup import cleanup_partial_single_copy,work_density
from test_sad_pipeline import mtz

SEQ=''.join(np.random.default_rng(125).choice(list('ACDEFGHIKLMNPQRSTVWY'),80))
class Density:
 def interpolate_value(self,pos):return 2.-.1*abs(pos.y)
def structure(parts):
 st=gemmi.Structure();st.cell=gemmi.UnitCell(400,60,60,90,90,90);st.spacegroup_hm='P 1';model=gemmi.Model(1)
 for ci,(start,text,shift) in enumerate(parts):
  ch=gemmi.Chain(chr(65+ci))
  for j,letter in enumerate(text):
   r=gemmi.Residue();r.name='UNK' if letter=='X' else gemmi.expand_one_letter(letter,gemmi.ResidueKind.AA);r.seqid=gemmi.SeqId(j+1,' ')
   for name,dx in [('N',-1.2),('CA',0),('C',1.2),('O',1.5)]:
    a=gemmi.Atom();a.name=name;a.element=gemmi.Element('N' if name=='N' else 'O' if name=='O' else 'C');a.pos=gemmi.Position(3.8*(start+j)+dx,shift,0);r.add_atom(a)
   ch.add_residue(r)
  model.add_chain(ch)
 st.add_model(model);return st

def test_unique_partial_anchors_remove_duplicates_junk_without_changing_atoms():
 st=structure([(0,SEQ[:50],0),(40,SEQ[40:],0),(0,'XXXXXXXX',20),(10,SEQ[10:16],20)]);before=st.make_mmcif_document().as_string();out,r=cleanup_partial_single_copy(st,SEQ,Density(),asu_copies=1)
 assert r['output_protein_residues']==80 and r['removed_protein_residues']==24
 assert r['overlapping_sequence_positions']==10 and len(out[0])==1
 assert [x.seqid.num for x in out[0][0]]==list(range(1,81))
 for residue in out[0][0]:assert residue.find_atom('CA','*').pos.x==pytest.approx(3.8*(residue.seqid.num-1))
 assert before==st.make_mmcif_document().as_string() and not r['selection_uses_Rfree']

def test_nonphysical_join_is_left_as_a_chain_break():
 out,r=cleanup_partial_single_copy(structure([(0,SEQ[:50],0),(50,SEQ[50:],20)]),SEQ,Density(),asu_copies=1)
 assert len(out[0])==2 and r['output_protein_residues']==80
 assert not r['forced_nonphysical_joins']

@pytest.mark.parametrize('copies',[0,2,True,1.0])
def test_explicit_copy_gate(copies):
 with pytest.raises(ValueError,match='single-copy'):cleanup_partial_single_copy(structure([(0,SEQ,0)]),SEQ,Density(),asu_copies=copies)

def test_ambiguous_sequence_and_insufficient_coverage_refused():
 with pytest.raises(ValueError,match='anchor'):cleanup_partial_single_copy(structure([(0,'A'*40,0)]),'A'*80,Density(),asu_copies=1)
 with pytest.raises(ValueError,match='coverage'):cleanup_partial_single_copy(structure([(0,SEQ[:45],0)]),SEQ,Density(),asu_copies=1)

def test_real_density_excludes_free_coefficients_and_unobserved_hkl(tmp_path):
 h=np.array([[1,0,0],[0,1,0],[0,0,1]],dtype=float);native=tmp_path/'native.mtz';mtz(native,[('F','F'),('SIGF','Q'),('FreeR_flag','I')],np.column_stack([h,[10,20,30],[1,1,1],[1,7,1]]))
 maps=[]
 for i,value in enumerate([0,1e10]):
  p=tmp_path/('map%d.mtz'%i);mtz(p,[('FWT','F'),('PHWT','P')],np.column_stack([h,[10,value,30],[20,40,60]]));g,r=work_density(p,native,free_value=7);assert r['free_coefficients_removed']==1 and r['measured_work_coefficients']==2;maps.append(np.array(g,copy=True))
 np.testing.assert_array_equal(*maps)

def test_short_protein_split_across_chains_does_not_need_a_40_residue_fragment():
 seq=SEQ[:51]
 st=structure([(0,seq[:21],0),(22,seq[22:],20)])
 cleaned,report=cleanup_partial_single_copy(st,seq,Density(),asu_copies=1,min_long_anchor=None)
 assert report['thresholds']['min_long_anchor']==26
 assert report['output_protein_residues']==50 and report['coverage']==pytest.approx(50/51)
 assert len(cleaned[0])==2 and not report['forced_nonphysical_joins']
 with pytest.raises(ValueError,match='anchor'):
  cleanup_partial_single_copy(st,seq,Density(),asu_copies=1,min_long_anchor=40)

def test_short_anchor_adaptation_keeps_coverage_and_ambiguity_gates():
 seq=SEQ[:51]
 with pytest.raises(ValueError,match='coverage'):
  cleanup_partial_single_copy(structure([(0,seq[:29],0)]),seq,Density(),asu_copies=1,min_long_anchor=None)
 with pytest.raises(ValueError,match='anchor'):
  cleanup_partial_single_copy(structure([(0,'A'*30,0)]),'A'*51,Density(),asu_copies=1,min_long_anchor=None)
