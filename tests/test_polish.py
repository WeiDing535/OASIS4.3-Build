from types import SimpleNamespace
import numpy as np
from phasebuild.geometry import template
from phasebuild.polish import fit_terminals, rotate
from phasebuild.validation import geometry_report


def test_recover_free_terminal_n_without_moving_ca_or_internal_peptide_links():
    expected=template(-63.,-42.,length=3)+[10.,10.,10.]
    damaged=expected.copy()
    damaged[0,0]=rotate(damaged[0,0:1],damaged[0,1],damaged[0,2]-damaged[0,1],[90.])[0,0]
    def score(xyz):
        return -np.sum((xyz[...,None,:,:]-expected)**2,axis=(-1,-2)).min(axis=-1)
    density=SimpleNamespace(score=score,distance=lambda a,b:float(np.linalg.norm(a-b)))
    fitted,report=fit_terminals([damaged],density)
    assert report['changes']
    np.testing.assert_allclose(fitted[0][:,1],damaged[:,1],atol=1e-8)
    np.testing.assert_allclose(fitted[0],expected,atol=1e-6)
    assert not geometry_report(fitted)['severe_geometry_issues']
