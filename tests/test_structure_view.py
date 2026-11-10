"""The 3D page's geometry: traces, rotation, fitting and projection."""
import numpy as np
import pytest

from phasebuild import structure_view


def test_chain_traces_keep_drawable_chains_only():
    preview=[{'chain':'A','points':[[0,0,0,'1','MET'],[3.8,0,0,'2','ALA'],[7.6,0,0,'3','GLY']]},
             {'chain':'B','points':[[0,5,0,'1','SER']]},
             {'chain':'C','points':[]}]
    traces=structure_view.chain_traces(preview)
    assert [chain for chain,_ in traces]==['A']
    assert traces[0][1].shape==(3,3)
    assert traces[0][1][1].tolist()==[3.8,0.0,0.0]
    assert structure_view.chain_traces(None)==[]


def test_chain_traces_skip_rows_that_are_not_coordinates():
    preview=[{'chain':'A','points':[[0,0,0],[1,1,1],['x','y','z'],[2,2,2]]}]
    assert structure_view.chain_traces(preview)[0][1].tolist()==[[0,0,0],[1,1,1],[2,2,2]]


def test_rotation_is_orthonormal_and_keeps_the_axis():
    matrix=structure_view.rotation(0.7,-0.4)
    assert np.allclose(matrix@matrix.T,np.eye(3))
    assert np.isclose(np.linalg.det(matrix),1.0)
    assert np.allclose(structure_view.rotation(0,0),np.eye(3))


def test_fit_puts_every_point_inside_the_canvas():
    points=np.array([[0,0,0],[20,0,0],[0,15,0],[0,0,12]],float)
    center,scale=structure_view.fit(points,400,300)
    x,y,_=structure_view.project(points,np.eye(3),center,scale,400,300)
    assert (x>=26-1e-6).all() and (x<=400-26+1e-6).all()
    assert (y>=26-1e-6).all() and (y<=300-26+1e-6).all()


def test_projection_centre_maps_to_the_middle_and_flips_y():
    center,scale=structure_view.fit(np.array([[0,0,0],[1,0,0]],float),200,100)
    x,y,_=structure_view.project(np.array([[center[0],center[1],center[2]]]),np.eye(3),center,scale,200,100)
    assert np.allclose([x[0],y[0]],[100,50])
    x,y,_=structure_view.project(np.array([[center[0],center[1]+2,center[2]]]),np.eye(3),center,scale,200,100)
    assert y[0]<50


def test_empty_input_is_harmless():
    center,scale=structure_view.fit(np.zeros((0,3)),200,200)
    assert scale>0 and center.tolist()==[0,0,0]


def test_read_atoms_reads_a_small_model(tmp_path):
    model=tmp_path/'tiny.pdb'
    model.write_text(
        'ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00 20.00           N\n'
        'ATOM      2  CA  MET A   1       1.458   0.000   0.000  1.00 20.00           C\n'
        'ATOM      3  CA  GLY A   2       2.916   0.000   0.000  1.00 20.00           C\n'
        'END\n')
    chains=structure_view.read_atoms(model)
    assert [chain for chain,_,_ in chains]==['A']
    _,coordinates,elements=chains[0]
    assert coordinates.shape==(3,3) and elements.tolist()==['N','C','C']
    with pytest.raises(Exception):
        structure_view.read_atoms(tmp_path/'missing.pdb')
