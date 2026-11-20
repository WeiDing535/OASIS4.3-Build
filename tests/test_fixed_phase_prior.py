from pathlib import Path
from types import SimpleNamespace
import json
import gemmi
import numpy as np
import pytest
from phasebuild.fixed_phase_prior import make_fixed_prior_refiner


def test_factory_rejects_wrong_hash_and_later_source_changes(tmp_path):
    source = tmp_path / 'source.mtz'
    source.write_bytes(b'original phases')
    with pytest.raises(ValueError, match='hash'):
        make_fixed_prior_refiner(source, expected_sha256='0' * 64)
    refine = make_fixed_prior_refiner(source)
    source.write_bytes(b'changed phases')
    with pytest.raises(RuntimeError, match='changed during recycling'):
        refine(tmp_path / 'unused_native.mtz', tmp_path / 'unused_model.cif', tmp_path / 'out')


def test_factory_reuses_audited_fixed_prior_without_changing_output_contract(tmp_path, monkeypatch):
    import phasebuild.fixed_phase_prior as adapter
    def write(path, columns, data):
        mtz = gemmi.Mtz(with_base=True)
        mtz.spacegroup = gemmi.find_spacegroup_by_name('P 1')
        mtz.set_cell_for_all(gemmi.UnitCell(40, 40, 40, 90, 90, 90))
        mtz.add_dataset('test')
        for name, kind in columns:
            mtz.add_column(name, kind)
        mtz.set_data(np.array(data, dtype=np.float32))
        mtz.update_reso()
        mtz.write_to_file(str(path))
    native, source = tmp_path / 'native.mtz', tmp_path / 'phase.mtz'
    write(native, [('F', 'F'), ('SIGF', 'Q'), ('FreeR_flag', 'I')],
          [[1, 0, 0, 10, 1, 0], [2, 0, 0, 20, 1, 1]])
    write(source, [(label, 'A') for label in ['HLA', 'HLB', 'HLC', 'HLD']],
          [[1, 0, 0, 2, 3, 1, -1]])
    inputs = []
    def fake_refine(merged, model, output, **kwargs):
        inputs.append(Path(merged))
        assert kwargs['use_prior'] and kwargs['cycles'] == 10
        out = Path(output)
        out.mkdir(exist_ok=False)
        (out / 'model.cif').write_text('mock model')
        return {'status': 'success', 'statistics': {'final': {'cycle': 10}}}
    monkeypatch.setattr(adapter, 'refmac_with_experimental_prior', fake_refine)
    refine = make_fixed_prior_refiner(source)
    for round_number in [1, 2]:
        output = tmp_path / ('round_%d_refine' % round_number)
        result = refine(native, tmp_path / 'model.cif', output)
        assert result['phase_merge_audit']['matched_prior_rows'] == 1
        assert result['fixed_original_experimental_prior']['path'] == str(source.resolve())
        assert (output / 'refinement_report.json').exists()
        assert (output / 'model.cif').exists()
    assert inputs[0] == inputs[1]
