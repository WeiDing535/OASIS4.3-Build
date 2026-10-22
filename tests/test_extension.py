from types import SimpleNamespace
import gemmi
import numpy as np
import pytest
from phasebuild.extension import free_positions,terminal_candidates
from phasebuild.geometry import template


def synthetic_density(xyz):
    points=xyz.reshape(-1,3)
    cell=gemmi.UnitCell(80,80,80,90,90,90)
    def sample(p):
        distances=np.sum((np.asarray(p)[...,None,:]-points)**2,axis=-1)
        return 3.5*np.exp(-distances.min(axis=-1)/.18)
    return SimpleNamespace(sample=sample,score=lambda p:sample(p).mean(axis=-1),
                           cell=cell,sg=gemmi.find_spacegroup_by_name('P 1'),
                           distance=lambda a,b:cell.find_nearest_image(gemmi.Position(*a),gemmi.Position(*b)).dist())


@pytest.mark.parametrize('direction,phi,psi,letter',[('C',-63.,-42.,'A'),('N',-63.,-42.,'A'),('C',60.,40.,'G')])
def test_terminal_growth_finds_an_ideal_chain_and_keeps_interior_anchor_fixed(direction,phi,psi,letter):
    ideal=template(phi,psi,length=8)+[35.,35.,35.]
    density=synthetic_density(ideal)
    if direction=='C':
        anchor=ideal[:4];positions=[4,5]
    else:
        anchor=ideal[4:];positions=[3,2]
    candidates,_=terminal_candidates(anchor,positions,letter*8,direction,density,None,
                                     beam_width=2,maximum_refinements=4)
    assert candidates
    best=candidates[0]
    assert best['added']==2
    if direction=='C':
        np.testing.assert_allclose(best['xyz'][:4,:3],anchor[:,:3],atol=1e-6)
        assert np.max(np.linalg.norm(best['xyz'][4:,:3]-ideal[4:6,:3],axis=-1))<.3
    else:
        np.testing.assert_allclose(best['xyz'][-4:],anchor,atol=1e-6)
        assert np.max(np.linalg.norm(best['xyz'][:2,:3]-ideal[2:4,:3],axis=-1))<.3


def test_sequence_occupancy_and_bounds_limit_extensions_without_known_chain_breaks():
    records=[('catalog','A'*20)]
    a={'record':0,'start':0,'end':7};b={'record':0,'start':9,'end':16}
    assert free_positions(a,[a,b],records,'C')==[7,8]
    assert free_positions(b,[a,b],records,'N')==[8,7]
    assert free_positions(a,[a,b],records,'N')==[]
    assert free_positions(b,[a,b],records,'C')==[16,17,18,19]


def test_no_growth_is_forced_through_a_weak_map():
    anchor=template(-63.,-42.,length=4)+[35.,35.,35.]
    d=synthetic_density(anchor)
    d.sample=lambda p:np.zeros(p.shape[:-1]);d.score=lambda p:np.zeros(p.shape[:-2])
    candidates,_=terminal_candidates(anchor,[4,5],'A'*6,'C',d,None,maximum_refinements=2)
    assert not candidates
