import copy
import json
from pathlib import Path
import pytest
from phasebuild.data import fingerprint
from phasebuild.hybrid_recycle import reuse_first_build


@pytest.fixture
def cached(tmp_path):
    old = tmp_path / 'old'
    source = old / 'round_01_build'
    source.mkdir(parents=True)
    command = ['/same/runtime/cbuccaneer', '-mtzin', '/same/map.mtz', '-pdbin', '/same/seed.cif',
               '-pdbout', str(source / 'model.pdb'), '-xmlout', str(source / 'buccaneer.xml'),
               '-cycles', '3', '-jobs', '1', '-model-filter', '-no-correlation-mode']
    (source / 'command.json').write_text(json.dumps(command))
    (source / 'model.pdb').write_text('existing construction PDB')
    (source / 'model.cif').write_text('existing construction CIF')
    (source / 'buccaneer.log').write_text('cbuccaneer: Normal termination')
    (source / 'buccaneer.xml').write_text('<Buccaneer><Cycles>' + ''.join(
        '<Cycle><Number>%d</Number></Cycle>' % i for i in range(1, 4)) + '</Cycles></Buccaneer>')
    report = {'parameters': {'rounds': 5, 'refine_cycles': 10},
              'inputs': {'model': {'path': '/same/seed.cif', 'sha256': 'seed'}},
              'binary': {'sha256': 'binary'}, 'training_library': [{'sha256': 'training'}]}
    (old / 'hybrid_report.json').write_text(json.dumps(report))
    destination = tmp_path / 'new/round_01_build'
    expected = list(command)
    expected[expected.index('-pdbout') + 1] = str(destination / 'model.pdb')
    expected[expected.index('-xmlout') + 1] = str(destination / 'buccaneer.xml')
    return source, destination, expected, report


def test_verified_recovery_retains_original_files_and_command(cached):
    source, destination, expected, report = cached
    result = reuse_first_build(source, destination, expected, report)
    assert result['construction_executed_again'] is False
    assert (destination / 'command.json').read_bytes() == (source / 'command.json').read_bytes()
    assert result['copied_model']['sha256'] == fingerprint(source / 'model.cif')['sha256']
    assert all(row['original']['sha256'] == row['copied']['sha256'] for row in result['files'])


@pytest.mark.parametrize('changed', ['input', 'binary', 'training', 'budget', 'cycles'])
def test_cached_first_construction_is_rejected_on_provenance_or_protocol_change(cached, changed):
    source, destination, command, report = cached
    current = copy.deepcopy(report)
    if changed == 'input': current['inputs']['model']['sha256'] = 'other'
    elif changed == 'binary': current['binary']['sha256'] = 'other'
    elif changed == 'training': current['training_library'][0]['sha256'] = 'other'
    elif changed == 'budget': current['parameters']['refine_cycles'] = 20
    else: command[command.index('-cycles') + 1] = '2'
    with pytest.raises(ValueError):
        reuse_first_build(source, destination, command, current)
    assert not destination.exists()


def test_incomplete_construction_cannot_be_reused(cached):
    source, destination, command, report = cached
    (source / 'buccaneer.log').write_text('failed before normal termination')
    with pytest.raises(ValueError, match='did not complete'):
        reuse_first_build(source, destination, command, report)
    assert not destination.exists()
