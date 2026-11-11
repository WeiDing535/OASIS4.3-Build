"""Probability phases must not inherit REFMAC difference-map phase rotations."""
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import gemmi
import pytest

from phasebuild import hybrid_recycle


def test_work_round_guard_uses_work_only():
    stages=[{'statistics':{'r_work':.40,'r_free':.90}},
            {'statistics':{'r_work':.34,'r_free':.95}},
            {'statistics':{'r_work':.43,'r_free':.01}}]
    assert hybrid_recycle.select_round_by_work_fit(stages)==1
    stages[-1]['statistics']['r_work']=.38
    assert hybrid_recycle.select_round_by_work_fit(stages)==1
    stages[-1]['statistics']['r_work']=.345
    assert hybrid_recycle.select_round_by_work_fit(stages)==2
from phasebuild.data import audit
from phasebuild.geometry import template
from phasebuild.trace import write_backbone
from test_phase_refinement import make_mtz


MAP_COLUMNS = [('F', 'F'), ('SIGF', 'Q'), ('FWT', 'F'),
               ('PHWT', 'P'), ('FOM', 'W')]
MAP_ROWS = [[1, 0, 0, 10, 1, 5, 30, .5],
            [2, 0, 0, 20, 1, 10, 60, .7]]


@pytest.mark.parametrize('phase_source', ['combined', 'calculated', 'both', 'legacy'])
def test_recycling_keeps_centroid_probability_and_map_channels_distinct(
        tmp_path, monkeypatch, phase_source):
    root = tmp_path / 'runtime'
    (root / 'bin').mkdir(parents=True)
    for name in ['cbuccaneer', 'refmac5']:
        binary = root / 'bin' / name
        binary.write_text('mock engine')
        binary.chmod(0o755)
    library = root / 'lib/data/reference_structures/reference-1tqw'
    library.parent.mkdir(parents=True)
    for suffix in ['.mtz', '.pdb']:
        library.with_suffix(suffix).write_text('independent mock reference')
    maps = tmp_path / 'initial.mtz'
    source = make_mtz(maps, MAP_COLUMNS, MAP_ROWS)
    sequence = tmp_path / 'sequence.fasta'
    sequence.write_text('>A\nAA\n')
    model = tmp_path / 'builder.cif'
    write_backbone(model, [template(-63., -42., length=2) + [20, 20, 20]],
                   SimpleNamespace(cell=source.cell, sg=source.spacegroup))
    gemmi.read_structure(str(model)).write_pdb(str(tmp_path / 'builder.pdb'))
    commands = []

    class FakeProcess:
        def __init__(self, command, **kwargs):
            commands.append(command)
            shutil.copyfile(tmp_path / 'builder.pdb',
                            command[command.index('-pdbout') + 1])
            self.pid = 999999

        def wait(self, timeout=None):
            return 0

        def poll(self):
            return 0

    def refine(native, built, output, cycles, ccp4=None, **kwargs):
        out = Path(output)
        out.mkdir()
        shutil.copyfile(built, out / 'model.cif')
        extra = {'combined': [('PHCOMB', 'P')],
                 'calculated': [('PHIC', 'P')],
                 'both': [('PHIC', 'P'), ('PHCOMB', 'P')],
                 'legacy': []}[phase_source]
        columns = MAP_COLUMNS + extra
        # A signed map phase can flip by 180 degrees from the FOM phase.
        # Distinct PHIC/PHCOMB values also exercise the explicit priority.
        rows = [row + [row[6] + (120. if label == 'PHCOMB' and phase_source == 'both'
                                else 180.) for label, _ in extra]
                for row in MAP_ROWS]
        make_mtz(out / 'refined.mtz', columns, rows)
        return dict(statistics=dict(final=dict(cycle=cycles, r_work=.3, r_free=.4)),
                    free_flags_preserved=True)

    monkeypatch.setattr(hybrid_recycle.subprocess, 'Popen', FakeProcess)
    monkeypatch.setattr(hybrid_recycle, 'refmac_refine', refine)
    report = hybrid_recycle.recycle(maps, maps, None, sequence, tmp_path / 'out',
                                   rounds=2, ccp4=root, progress=lambda _: None)
    assert len(commands) == 2
    assert commands[0][commands[0].index('-colin-phifom') + 1] == '/*/*/[PHWT,FOM]'
    expected = {'combined': 'PHCOMB', 'calculated': 'PHIC',
                'both': 'PHCOMB', 'legacy': 'PHWT'}[phase_source]
    assert commands[1][commands[1].index('-colin-phifom') + 1] == '/*/*/[%s,FOM]' % expected
    assert all(command[command.index('-colin-fc') + 1] == '/*/*/[FWT,PHWT]'
               for command in commands)
    handoffs = report['construction_phase_handoffs']
    assert [record['probability_phase'] for record in handoffs] == ['PHWT', expected]
    assert handoffs[1]['map_phase'] == 'PHWT'
    assert ('legacy' in handoffs[1]['probability_contract']) is (phase_source == 'legacy')
    saved = json.loads((tmp_path / 'out/hybrid_report.json').read_text())
    assert saved['construction_phase_handoffs'] == handoffs


@pytest.mark.parametrize('extra', [
    [('PHCOMB', 'P'), ('PHCOMB', 'P')],
    [('PHCOMB', 'F')],
    [('PHIC', 'P'), ('PHIC', 'P')],
    [('PHIC', 'F')],
    [('PHCOMB', 'P'), ('PHIC', 'F')],
    [('PHCOMB', 'F'), ('PHIC', 'P')],
    [('PHWT', 'P')],
    [('FOM', 'W')],
])
def test_malformed_or_ambiguous_probability_columns_are_rejected(tmp_path, extra):
    mtz = tmp_path / 'ambiguous.mtz'
    rows = [row + [45.] * len(extra) for row in MAP_ROWS]
    make_mtz(mtz, MAP_COLUMNS + extra, rows)
    with pytest.raises(ValueError, match='ambiguous|Ambiguous|incorrectly typed'):
        hybrid_recycle._construction_phase_columns(audit(mtz))


def test_missing_map_channel_is_not_hidden_by_valid_centroid(tmp_path):
    mtz = tmp_path / 'missing_map_phase.mtz'
    columns = [column for column in MAP_COLUMNS if column[0] != 'PHWT'] + [('PHCOMB', 'P')]
    rows = [row[:6] + row[7:] + [210.] for row in MAP_ROWS]
    make_mtz(mtz, columns, rows)
    with pytest.raises(ValueError, match='required map column PHWT'):
        hybrid_recycle._construction_phase_columns(audit(mtz))
