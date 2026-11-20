import gemmi
import numpy as np
from oasis_py.data import Reflections
import sys
from pathlib import Path
from support.experimental_sad_errors import missing_substructure_variance


def observations(extra):
    rng=np.random.default_rng(48)
    cell=gemmi.UnitCell(80,90,100,90,90,90);sg=gemmi.SpaceGroup('P 1')
    h=gemmi.make_miller_array(cell,sg,5.,unique=True);n=len(h)
    q=rng.normal(size=n)+1j*rng.normal(size=n)
    difference=rng.normal(size=n)*np.sqrt(1+2*.7**2*abs(q)**2+extra)
    data=Reflections(h,np.ones(n)*100,np.ones(n),cell,sg,
                     dano=difference,sigdano=np.ones(n))
    return data,q


def test_unmodeled_variance_is_detected_and_scale_is_recovered():
    data,q=observations(6.)
    sf,sigma,meta=missing_substructure_variance(data,q)
    assert meta['accepted'] and not meta['target_reference_used']
    assert abs(np.median(abs(sf/q))-.7)<.1
    assert abs(np.median(sigma**2-1)-6)<1
    assert np.all(sigma>=data.sigdano)


def test_holdout_rejects_unneeded_shell_complexity_and_is_order_invariant():
    data,q=observations(0.)
    sf,sigma,meta=missing_substructure_variance(data,q)
    assert not meta['accepted']
    np.testing.assert_array_equal(sigma,data.sigdano)
    order=np.random.default_rng(77).permutation(len(q))
    other,other_sigma,other_meta=missing_substructure_variance(data.subset(order),q[order])
    assert other_meta['accepted']==meta['accepted']
    np.testing.assert_allclose(other,sf[order],rtol=1e-7)
    np.testing.assert_array_equal(other_sigma,sigma[order])
