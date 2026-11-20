"""Refinement evidence must not change with phase-improvement iterations."""
import numpy as np
import pytest
from test_sad_pipeline import inputs
from oasis_py.data import read_mtz, read_atoms
from oasis_py.sequence import sequence_asu
from oasis_py.sad import phase_sad_likelihood
from oasis_py.phasing import phase_native
from phasebuild.observation_prior import observation_prior


def arguments(inputs):
    data, _ = read_mtz(inputs['sad'], 'sad', free_label='FreeR_flag', free_value=7)
    heavy = read_atoms(inputs['heavy'], data)
    sequence = sequence_asu(inputs['sequence'], 2, data.cell, data.spacegroup)
    return data, heavy, sequence['asu_composition'], {'S': .56}


def test_refinement_likelihood_is_invariant_to_triplet_recycling(inputs):
    data, heavy, composition, fpp = arguments(inputs)
    first = phase_sad_likelihood(data, composition, heavy, fpp, cycles=1, max_neighbors=80)
    recycled = phase_sad_likelihood(data, composition, heavy, fpp, cycles=3,
                                    triplet_weight=.7, max_neighbors=80)
    assert not np.array_equal(first.hl, recycled.hl)
    saved = recycled.hl.copy()
    a = observation_prior(data, first, heavy, fpp, 'likelihood')
    b = observation_prior(data, recycled, heavy, fpp, 'likelihood')
    np.testing.assert_array_equal(a.hl, b.hl)
    np.testing.assert_array_equal(a.phase, b.phase)
    np.testing.assert_array_equal(a.fom, b.fom)
    np.testing.assert_array_equal(recycled.hl, saved)
    assert not np.array_equal(b.hl, recycled.hl)
    assert b.metadata['structural_triplet_field_included'] is False
    assert b.probability is None


def test_pplus_observation_hl_is_retained_without_branch_posterior(inputs):
    data, heavy, composition, fpp = arguments(inputs)
    result = phase_native(data, 'sad', composition, heavy=heavy, fpp=fpp,
                          cycles=2, max_neighbors=80)
    prior = observation_prior(data, result, heavy, fpp, 'pplus')
    np.testing.assert_array_equal(prior.hl, result.hl)
    assert not np.shares_memory(prior.hl, result.hl)
    assert prior.probability is None
    assert np.isfinite(prior.phase).all() and np.all((prior.fom >= 0) & (prior.fom <= 1))


def test_likelihood_prior_rejects_unidentified_posterior(inputs):
    data, heavy, composition, fpp = arguments(inputs)
    result = phase_native(data, 'sad', composition, heavy=heavy, fpp=fpp,
                          cycles=1, max_neighbors=80)
    with pytest.raises(ValueError, match='matching SAD likelihood metadata'):
        observation_prior(data, result, heavy, fpp, 'likelihood')


@pytest.mark.parametrize('collision',['existing','same_as_phase'])
def test_phase_export_refuses_prior_output_collision_before_writing(tmp_path,collision):
    import run_oasis_stage
    target=tmp_path/'phase.mtz'
    prior=tmp_path/'prior.mtz' if collision=='existing' else target
    if collision=='existing':prior.write_bytes(b'preserve existing prior')
    with pytest.raises(FileExistsError if collision=='existing' else ValueError):
        run_oasis_stage.main(['phase','--sad','unused','--sequence','unused','--heavy','unused',
            '--fpp','S=0.56','--output',str(target),'--observation-output',str(prior)])
    assert not target.exists()
    if collision=='existing':assert prior.read_bytes()==b'preserve existing prior'
