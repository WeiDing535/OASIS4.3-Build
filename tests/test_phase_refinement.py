import gemmi
import numpy as np
import pytest
from phasebuild.phase_refinement import merge_experimental_phases


def make_mtz(path, columns, rows):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name('P 1')
    mtz.set_cell_for_all(gemmi.UnitCell(40., 40., 40., 90., 90., 90.))
    mtz.add_dataset('test')
    for name, kind in columns:
        mtz.add_column(name, kind)
    mtz.set_data(np.array(rows, dtype=np.float32))
    mtz.update_reso()
    mtz.write_to_file(str(path))
    return mtz


def test_phase_join_preserves_native_order_flags_and_missing_high_resolution(tmp_path):
    native = tmp_path / 'native.mtz'
    # Deliberately non-sorted HKLs and a reordered free column.
    original = make_mtz(native, [('FreeR_flag', 'I'), ('F', 'F'), ('SIGF', 'Q')],
                       [[8, 0, 0, 7, 30, 2], [1, 0, 0, 0, 10, 1], [2, 0, 0, 3, 20, 1]])
    phase = tmp_path / 'phase.mtz'
    make_mtz(phase, [('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A')],
             [[1, 0, 0, 2, 3, -1, 4], [2, 0, 0, 1, np.nan, 0, 0], [9, 0, 0, 4, 2, 0, 0]])
    output = tmp_path / 'merged.mtz'
    result = merge_experimental_phases(native, phase, output)
    data = np.array(gemmi.read_mtz_file(str(output)), copy=False)
    np.testing.assert_array_equal(data[:, :6], np.array(original, copy=False))
    assert np.isnan(data[[0, 2], 6:]).all()
    np.testing.assert_array_equal(data[1, 6:], [2, 3, -1, 4])
    assert result['matched_prior_rows'] == 1 and result['missing_prior_rows'] == 2
    assert result['source_rows_not_in_native'] == 1
    assert result['matched_free_prior_rows'] == 1


def test_phase_join_rejects_duplicates_and_invalid_fom(tmp_path):
    native = tmp_path / 'native.mtz'
    make_mtz(native, [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')],
             [[1, 0, 0, 10, 1, 0], [2, 0, 0, 20, 1, 3]])
    phase = tmp_path / 'phase.mtz'
    columns = [('PHIB', 'P'), ('FOM', 'W')]
    make_mtz(phase, columns, [[1, 0, 0, 120, 1.5]])
    with pytest.raises(ValueError, match='FOM outside'):
        merge_experimental_phases(native, phase, tmp_path / 'one.mtz', 'phifom')
    make_mtz(phase, columns, [[1, 0, 0, 120, .5], [1, 0, 0, 110, .6]])
    with pytest.raises(ValueError, match='Duplicate experimental'):
        merge_experimental_phases(native, phase, tmp_path / 'two.mtz', 'phifom')


def test_long_path_refinement_stages_identical_input_and_uses_short_relative_arguments(tmp_path, monkeypatch):
    from pathlib import Path
    from types import SimpleNamespace
    import phasebuild.phase_refinement as adapter
    from phasebuild.geometry import template
    from phasebuild.trace import write_backbone
    long_parent = tmp_path / ('a' * 90) / ('b' * 90) / ('c' * 90)
    long_parent.mkdir(parents=True)
    native = long_parent / ('native_' + 'd' * 80 + '.mtz')
    original = make_mtz(native, [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I'),
        ('SAD_HLA', 'A'), ('SAD_HLB', 'A'), ('SAD_HLC', 'A'), ('SAD_HLD', 'A')],
        [[1, 0, 0, 10, 1, 0, 2, 3, 1, -1], [2, 0, 0, 20, 1, 1, np.nan, np.nan, np.nan, np.nan]])
    assert len(str(native)) > 255
    model = long_parent / 'model.cif'
    write_backbone(model, [template(-63., -42., length=2) + [20., 20., 20.]],
                   SimpleNamespace(cell=original.cell, sg=original.spacegroup))
    ccp4 = tmp_path / 'ccp4'
    (ccp4 / 'bin').mkdir(parents=True)
    (ccp4 / 'bin/refmac5').write_bytes(b'fake executable')
    (ccp4 / 'bin/refmac5').chmod(0o755)
    def fake_run(command, **kwargs):
        expected = ['HKLIN', 'input.mtz', 'XYZIN', 'input.pdb', 'HKLOUT', 'refined.mtz',
                    'XYZOUT', 'refined.pdb', 'XMLOUT', 'refmac.xml']
        assert command[1:] == expected
        assert all(len(value) < 32 and not Path(value).is_absolute() for value in command[2::2])
        out = Path(kwargs['cwd'])
        assert (out / 'input.mtz').read_bytes() == native.read_bytes()
        gemmi.read_structure(str(out / 'input.pdb')).write_pdb(str(out / 'refined.pdb'))
        make_mtz(out / 'refined.mtz', [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')],
                 [[1, 0, 0, 10, 1, 0], [2, 0, 0, 20, 1, 1]])
        (out / 'refmac.xml').write_text('''<REFMAC><err_level>0</err_level><Overall_stats>
<n_reflections_all>2</n_reflections_all><n_reflections_free>1</n_reflections_free>
<resolution_low>40</resolution_low><resolution_high>20</resolution_high>
<stats_vs_cycle><new_cycle><cycle>0</cycle><r_factor>0.437</r_factor><r_free>0.433</r_free>
</new_cycle></stats_vs_cycle></Overall_stats></REFMAC>''')
        kwargs['stdout'].write('Residual : Phase modified Maximum Likelihood for Fs\n' +
            'Ncyc    Rfact    Rfree\n0 0.4375 0.4333 0.58 3000 150 0.0148 1.27 1.56 0.89 0.1\n$TEXT:')
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(adapter.subprocess, 'run', fake_run)
    result = adapter.refmac_with_experimental_prior(native, model, long_parent / 'out', cycles=0, ccp4=ccp4)
    assert result['input']['sha256'] == result['staged_input']['sha256']
    assert result['free_flags_preserved']
    assert result['phase_modified_likelihood_active']
