import gemmi
import numpy as np
import pytest
from phasebuild.geometry import template
from phasebuild.discriminant_assembly import positive_evidence,assemble_discriminant_fragments
from phasebuild.validation import geometry_report


class ConstantMap:
    def __init__(self):
        self.cell=gemmi.UnitCell(80,80,80,90,90,90)
        self.sg=gemmi.find_spacegroup_by_name('P 1')
    def score(self,xyz):
        return np.full(len(xyz),3.)
    def distance(self,a,b):
        return self.cell.find_nearest_image(gemmi.Position(*a),gemmi.Position(*b)).dist()


class ConstantClassifier:
    metadata={'calibration':{'thresholds':{'0.05':{'logit':0.}}}}
    def score(self,density,xyz):
        return np.full(len(xyz),3.)


def test_positive_evidence_counts_only_fixed_threshold_excess():
    assert positive_evidence([-10.,-.1,1.,3.],.5)==pytest.approx(3.)


def test_classifier_paths_keep_peptide_geometry_and_count_each_residue_once():
    chain=template(-120.,130.,7)+[30.,30.,30.]
    chains,report=assemble_discriminant_fragments(ConstantMap(),
                 [chain[i:i+3] for i in range(5)],ConstantClassifier(),capacity=7)
    assert len(chains)==1 and len(chains[0])==7
    assert np.allclose(chains[0][:,:3],chain[:,:3])
    assert report['classifier_objective']==pytest.approx(21.)
    assert geometry_report(chains)['severe_geometry_issues']==[]


def test_classifier_ranking_can_recover_alternative_hidden_by_higher_map_score():
    chain=template(-120.,130.,3)+[30.,30.,30.]
    wrong=chain+[0.,0.,.4]
    def nearest(xyz):
        return np.sqrt(np.mean((xyz[:,None,:3]-chain[None,:,:3])**2,axis=(-1,-2))).min(axis=1)
    class AdversarialMap(ConstantMap):
        def score(self,xyz):
            return 3.+2.*nearest(np.asarray(xyz))
    class GeometryClassifier(ConstantClassifier):
        def score(self,density,xyz):
            return 3.-5.*nearest(np.asarray(xyz))
    density=AdversarialMap()
    assert density.score(wrong).sum()>density.score(chain).sum()
    output,report=assemble_discriminant_fragments(density,[wrong,chain],
                                                 GeometryClassifier(),capacity=3)
    assert len(output)==1 and np.allclose(output[0],chain)
    assert report['classifier_supported_pieces']==2


def test_unrelated_calibration_entry_is_fixed_and_unsupported_candidates_drop():
    chain=template(-120.,130.,3)+[30.,30.,30.]
    class NegativeClassifier(ConstantClassifier):
        def score(self,density,xyz):
            return np.full(len(xyz),-1.)
    output,report=assemble_discriminant_fragments(ConstantMap(),[chain],
                                                 NegativeClassifier(),capacity=3)
    assert output==[] and report['classifier_supported_pieces']==0
    with pytest.raises(ValueError):
        assemble_discriminant_fragments(ConstantMap(),[chain],ConstantClassifier(),
                                         capacity=3,calibration_decoy_rate='0.1')
