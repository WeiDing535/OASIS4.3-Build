import gemmi
import numpy as np
import pytest
from phasebuild.geometry import template
from phasebuild.fragment_join import append_torsions,overlap_graph,assemble_fragments
from phasebuild.validation import geometry_report


class ConstantMap:
    def __init__(self,spacegroup='P 1'):
        self.cell=gemmi.UnitCell(60,60,60,90,90,90)
        self.sg=gemmi.find_spacegroup_by_name(spacegroup)
    def score(self,xyz):
        return np.full(len(xyz),3.)
    def distance(self,a,b):
        return self.cell.find_nearest_image(gemmi.Position(*a),gemmi.Position(*b)).dist()


def test_overlap_extension_retains_ideal_geometry_and_original_curve():
    chain=template(-63.,-42.,9)+[20.,20.,20.]
    result=append_torsions(chain[:5],chain[3:],2)
    assert np.allclose(result[:,:3],chain[:,:3],atol=1e-8)
    assert geometry_report([result])['severe_geometry_issues']==[]


def test_overlap_uses_order_and_atom_names_not_only_ca_distance():
    chain=template(-120.,130.,6)+[20.,20.,20.]
    fragments=[chain[:3],chain[1:4],chain[2:5],chain[3:]]
    graph=overlap_graph(fragments,ConstantMap())
    assert any(edge[0]==1 and edge[1]==2 for edge in graph[0])
    reversed_fragments=[chain[:3],chain[1:4][::-1]]
    assert overlap_graph(reversed_fragments,ConstantMap())==[[],[]]
    atom_swapped=chain[1:4].copy();atom_swapped[:,[0,2]]=atom_swapped[:,[2,0]]
    assert overlap_graph([chain[:3],atom_swapped],ConstantMap())==[[],[]]


def test_assembly_joins_tripeptides_and_does_not_double_count_evidence():
    chain=template(-120.,130.,7)+[20.,20.,20.]
    fragments=[chain[i:i+3] for i in range(5)]
    merged,report=assemble_fragments(ConstantMap(),fragments,capacity=7)
    assert len(merged)==1 and len(merged[0])==7
    assert np.allclose(merged[0][:,:3],chain[:,:3])
    assert report['map_objective']==pytest.approx(7*.5)
    assert geometry_report(merged)['severe_geometry_issues']==[]


def test_join_handles_lattice_translation_and_crystallographic_image():
    density=ConstantMap('P 21 21 21')
    chain=template(-120.,130.,6)+[20.,20.,20.]
    op=list(density.sg.operations())[1]
    other=np.array([list(density.cell.orthogonalize(gemmi.Fractional(*op.apply_to_xyz(
        list(density.cell.fractionalize(gemmi.Position(*p)))))))
        for p in chain[1:4].reshape(-1,3)]).reshape(3,4,3)
    other+=[60.,-60.,0.]
    graph=overlap_graph([chain[:3],other],density)
    assert any(edge[0]==1 and edge[1]==2 for edge in graph[0])


def test_long_overlapping_traces_are_not_limited_to_short_seed_overlaps():
    chain=template(-120.,130.,15)+[30.,30.,30.]
    graph=overlap_graph([chain[:12],chain[1:13]],ConstantMap())
    assert any(edge[0]==1 and edge[1]==11 for edge in graph[0])


def test_assembly_respects_capacity_and_rejects_invalid_overlap():
    chain=template(-120.,130.,7)+[20.,20.,20.]
    result,_=assemble_fragments(ConstantMap(),[chain[i:i+3] for i in range(5)],capacity=4)
    assert sum(map(len,result))<=4
    with pytest.raises(ValueError):
        append_torsions(chain[:3],chain[1:4],1)


def test_map_objective_includes_previous_terminal_oxygen_change():
    class OxygenSensitiveMap(ConstantMap):
        def score(self,xyz):
            return 3.+.02*np.asarray(xyz)[:,3,2]
    density=OxygenSensitiveMap()
    chain=template(-120.,130.,7)+[20.,20.,20.]
    fragments=[chain[i:i+3].copy() for i in range(5)]
    # Terminal O has no following peptide restraint in an independently fitted
    # seed; a genuine join must replace it and account for the score difference.
    from phasebuild.geometry import oxygen
    for fragment in fragments:
        fragment[-1,3]=oxygen(*fragment[-1,:3],-40.)
    merged,report=assemble_fragments(density,fragments,capacity=7)
    assert len(merged)==1 and len(merged[0])==7
    expected=sum(np.sum(density.score(c)-2.5) for c in merged)
    assert report['map_objective']==pytest.approx(expected,abs=1e-10)
