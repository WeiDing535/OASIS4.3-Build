import gemmi
import numpy as np
from phasebuild.structure_pruning import prune_structure


class Map:
    cell=gemmi.UnitCell(40,40,40,90,90,90)
    sg=gemmi.find_spacegroup_by_name('P 1')


def structure():
    st=gemmi.Structure();st.cell=Map.cell;st.spacegroup_hm='P 1'
    model=gemmi.Model('1');chain=gemmi.Chain('A')
    for i in range(12):
        residue=gemmi.Residue();residue.name='MSE';residue.seqid=gemmi.SeqId(i+11,' ')
        for j,name in enumerate(['N','CA','C','O','SE']):
            atom=gemmi.Atom();atom.name=name;atom.element=gemmi.Element('Se' if name=='SE' else name[0])
            atom.pos=gemmi.Position(i,j,.125);atom.b_iso=17.25+i;atom.occ=.7
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain);st.add_model(model);return st


def test_prune_preserves_whole_residue_atoms_and_splits_removed_intervals():
    class Classifier:
        metadata={'calibration':{'thresholds':{'0.05':{'logit':0.}}}}
        def score(self,density,xyz):return np.array([3.,3.,3.,3.,-1.,-1.,3.,3.,3.,3.,3.,3.])
    st=structure();result,report=prune_structure(st,Map(),Classifier())
    assert len(st[0][0])==12
    assert [len(c) for c in result[0]]==[4,6]
    assert result[0][0].name=='A' and result[0][1].name!='A'
    for output,indices in zip(result[0],[range(4),range(6,12)]):
        for residue,i in zip(output,indices):
            old=st[0][0][i]
            assert residue.name=='MSE' and str(residue.seqid)==str(old.seqid)
            assert [(a.name,a.element.name,a.b_iso,a.occ,list(a.pos)) for a in residue]==[
                    (a.name,a.element.name,a.b_iso,a.occ,list(a.pos)) for a in old]
    assert report['output_residues']==10


def test_single_internal_weak_residue_keeps_original_chain_and_seqids():
    class Classifier:
        metadata={'calibration':{'thresholds':{'0.05':{'logit':0.}}}}
        def score(self,density,xyz):return np.array([3.,3.,3.,3.,-1.,3.,3.,3.,3.,3.,3.,3.])
    st=structure();result,report=prune_structure(st,Map(),Classifier())
    assert len(result[0])==1 and len(result[0][0])==12
    assert result[0][0].name=='A' and str(result[0][0][4].seqid)=='15'
