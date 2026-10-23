"""Exact-HKL SAD-to-building handoff; original observations/flags are retained."""
import json
from pathlib import Path
import gemmi
import numpy as np
from .data import fingerprint
from .refinement import reflection_partition

FLAG = 'PIPELINE_FREE'


def columns(mtz, schema):
    labels = mtz.column_labels()
    for name, kind in schema.items():
        if labels.count(name) != 1 or mtz.column_with_label(name).type != kind:
            raise ValueError('Required unique column/type: ' + name + '/' + kind)
    return labels


def keys(mtz):
    raw = np.array(mtz, copy=False)[:, :3]
    if not np.isfinite(raw).all() or not np.equal(raw, np.rint(raw)).all():
        raise ValueError('HKLs must be finite integers')
    result = [tuple(map(int, row)) for row in raw]
    if len(set(result)) != len(result):
        raise ValueError('Duplicate HKLs are not supported')
    return result


def crystal(first, second):
    if (first.spacegroup is None or second.spacegroup is None or
            first.spacegroup.hall != second.spacegroup.hall or
            not np.allclose(first.cell.parameters, second.cell.parameters, rtol=0, atol=.01)):
        raise ValueError('Crystal cell or space-group setting differs; explicit reindexing is required')


def inspect_inputs(sad_path, native_path, free_value=0):
    native, partition = reflection_partition(native_path, free_value=free_value)
    sad = gemmi.read_mtz_file(str(sad_path));crystal(native, sad)
    nk, sk = keys(native), keys(sad)
    # The upstream reader canonicalizes its inputs. Require the same convention
    # before joining so phase signs/origin shifts are never silently guessed.
    asu = gemmi.ReciprocalAsu(native.spacegroup);ops = native.spacegroup.operations()
    if any(tuple(asu.to_asu(list(k), ops)[0]) != k for k in nk + sk):
        raise ValueError('Unified input requires merged reciprocal-ASU HKLs in the same setting')
    n = np.array(native, copy=False);nl = native.column_labels()
    lookup = {k: n[i, nl.index('FreeR_flag')] for i, k in enumerate(nk)}
    joined = np.array([lookup.get(k, np.nan) for k in sk], dtype=np.float32)
    if FLAG in sad.column_labels():
        raise ValueError('SAD input already contains reserved pipeline flag column')
    if 'FreeR_flag' in sad.column_labels():
        columns(sad, {'FreeR_flag': 'I'})
        old = np.asarray(sad.column_with_label('FreeR_flag'))
        both = np.isfinite(old) & np.isfinite(joined)
        if np.any(old[both] != joined[both]):
            raise ValueError('SAD and native original free flags disagree by HKL')
    if not np.any(np.isfinite(joined) & (joined != free_value)):
        raise ValueError('SAD/native have no shared working reflections')
    return sad, native, joined, partition


def prepare_sad(sad_path, native_path, output, free_value=0):
    sad, native, flags, partition = inspect_inputs(sad_path, native_path, free_value)
    output = Path(output)
    if output.exists() or output.with_suffix('.json').exists():raise FileExistsError(output)
    original = np.array(sad, copy=True)
    dataset = sad.add_dataset('pipeline_original_free_flags')
    sad.add_column(FLAG, 'I', dataset.id)
    sad.set_data(np.column_stack([original, flags]))
    sad.write_to_file(str(output))
    report = {'sad': fingerprint(sad_path), 'native': fingerprint(native_path),
              'output': fingerprint(output), 'free_value': free_value,
              'native_partition': partition, 'flag_column': FLAG,
              'flags_regenerated': False, 'native_input_modified': False,
              'free_exclusion': 'Work reflections only in SAD and density modification; unmatched SAD rows excluded'}
    output.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    return report


def adapt_maps(native_path, experimental_path, dm_path, output, free_value=0):
    """Use native F times final DM FOM once, and retain raw SAD HL separately."""
    output = Path(output)
    if output.exists() or output.with_suffix('.audit.json').exists():raise FileExistsError(output)
    native, partition = reflection_partition(native_path, free_value=free_value)
    original = np.array(native, copy=True);nk = keys(native)
    phase = gemmi.read_mtz_file(str(dm_path));prior = gemmi.read_mtz_file(str(experimental_path))
    for m in [phase, prior]:crystal(native, m)
    pl = columns(phase, {'PHIB':'P', 'FOM':'W'});pk = keys(phase)
    hl = columns(prior, {name:'A' for name in ['HLA','HLB','HLC','HLD']});hk = keys(prior)
    phase_data = np.array(phase, copy=False);prior_data = np.array(prior, copy=False)
    dm = {k: phase_data[i, [pl.index('PHIB'),pl.index('FOM')]] for i,k in enumerate(pk)}
    raw = {k: prior_data[i, [hl.index(x) for x in ['HLA','HLB','HLC','HLD']]] for i,k in enumerate(hk)}
    nl = native.column_labels();flags = original[:,nl.index('FreeR_flag')]
    map_values = np.full((len(nk),3), np.nan, dtype=np.float32)
    matched = 0
    for i,k in enumerate(nk):
        # Free reflections must not re-enter the initial construction map.
        if flags[i] == free_value:continue
        if k in dm:
            p, m = dm[k]
            if not np.isfinite([p,m]).all() or not 0 <= m <= 1:
                raise ValueError('Invalid final DM phase/FOM')
            f = original[i,nl.index('F')]
            if np.isfinite(f):map_values[i] = [f*m,p,m];matched += 1
    if not matched:raise ValueError('No finite working native/map pairs match by HKL')
    # A missing original prior stays missing; the DM HL is never substituted.
    if not any(k in raw and np.isfinite(raw[k]).all() for k in nk):
        raise ValueError('No original SAD HL matches native HKLs')
    for name in ['FWT','PHWT','FOM']:
        if name in nl:raise ValueError('Native must be unmodified observations without map columns: '+name)
    dataset = native.add_dataset('pipeline_work_only_dm_map')
    for name,kind in [('FWT','F'),('PHWT','P'),('FOM','W')]:native.add_column(name,kind,dataset.id)
    native.set_data(np.column_stack([original,map_values]))
    native.write_to_file(str(output))
    restored = np.array(gemmi.read_mtz_file(str(output)), copy=False)
    if not np.array_equal(restored[:,:original.shape[1]],original,equal_nan=True):
        raise RuntimeError('Original native observations or flags changed')
    report = {'native':fingerprint(native_path),'experimental_prior':fingerprint(experimental_path),
              'dm_phases':fingerprint(dm_path),'output':fingerprint(output),
              'native_columns_unchanged':True,'native_order_unchanged':True,'matched_work_map_rows':matched,
              'free_value':free_value,'native_partition':partition,
              'map_formula':'FWT = native F * final DM FOM, once; PHWT = final DM PHIB',
              'prior_policy':'Separate pre-DM SAD HLA/HLB/HLC/HLD; never substitute DM/model posterior',
              'scope':'Upstream excludes free reflections. Later external-engine maps are not certified work-only.'}
    output.with_suffix('.audit.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
