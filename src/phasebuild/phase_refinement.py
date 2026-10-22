"""Opt-in REFMAC experimental-phase refinement, with exact-HKL audited priors.

This is a separate scientific protocol from amplitude-only refinement. Missing
experimental priors remain missing. No model phases or target references enter
the join or the adapter. LABIN activates the phase-aware likelihood in REFMAC;
REFI RESI MLKF remains the documented keyword in both arms.
"""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time
import xml.etree.ElementTree as ET
import gemmi
import numpy as np
from .data import fingerprint
from .refinement import reflection_partition, verify_free_flags
from .runtime import resolve_runtime, runtime_environment


def merge_experimental_phases(native_path, phase_path, output, representation='hl', free_value=0):
    if representation not in ['hl', 'phifom']:
        raise ValueError('Phase representation must be hl or phifom')
    native, partition = reflection_partition(native_path, free_value=free_value)
    phases = gemmi.read_mtz_file(str(phase_path))
    if (native.spacegroup.hm != phases.spacegroup.hm or
            not np.allclose(native.cell.parameters, phases.cell.parameters, rtol=0., atol=.01)):
        raise ValueError('Native and experimental phase crystal settings differ')
    labels = ['HLA', 'HLB', 'HLC', 'HLD'] if representation == 'hl' else ['PHIB', 'FOM']
    types = ['A'] * 4 if representation == 'hl' else ['P', 'W']
    for label, kind in zip(labels, types):
        if phases.column_labels().count(label) != 1 or phases.column_with_label(label).type != kind:
            raise ValueError('Required unique experimental phase column/type: ' + label)
    original = np.array(native, copy=True)
    source = np.array(phases, copy=False)
    if not np.all(source[:, :3] == np.rint(source[:, :3])):
        raise ValueError('Noninteger experimental HKLs')
    lookup = {tuple(row[:3].astype(int)): row for row in source}
    if len(lookup) != len(source):
        raise ValueError('Duplicate experimental HKLs')
    keys = [tuple(row[:3].astype(int)) for row in original]
    if len(set(keys)) != len(keys):
        raise ValueError('Duplicate native HKLs')
    priors = np.full((len(original), len(labels)), np.nan, dtype=np.float32)
    columns = [phases.column_labels().index(label) for label in labels]
    matched = 0
    for i, key in enumerate(keys):
        row = lookup.get(key)
        if row is None:
            continue
        values = row[columns]
        if not np.isfinite(values).all():
            continue
        if representation == 'phifom' and not 0. <= values[1] <= 1.:
            raise ValueError('Experimental FOM outside [0,1]')
        priors[i] = values
        matched += 1
    if not matched:
        raise ValueError('No exact-HKL finite experimental phases matched native data')
    dataset = native.add_dataset('original_experimental_phase_prior')
    destination_labels = ['SAD_' + label for label in labels]
    if any(label in native.column_labels() for label in destination_labels):
        raise ValueError('Native already has requested SAD prior columns')
    for label, kind in zip(destination_labels, types):
        native.add_column(label, kind, dataset.id)
    native.set_data(np.column_stack([original, priors]))
    native.update_reso()
    output = Path(output)
    if output.exists():
        raise ValueError('Merged MTZ output already exists')
    native.write_to_file(str(output))
    restored = np.array(gemmi.read_mtz_file(str(output)), copy=False)
    if not np.array_equal(restored[:, :original.shape[1]], original, equal_nan=True):
        raise RuntimeError('Native input columns changed during phase join')
    flags = original[:, native.column_labels().index('FreeR_flag')]
    missing = ~np.isfinite(priors).all(axis=1)
    distances = native.make_d_array()
    report = {'native': fingerprint(native_path), 'phase_source': fingerprint(phase_path),
        'output': fingerprint(output), 'representation': representation,
        'source_labels': labels, 'destination_labels': destination_labels,
        'native_partition': partition, 'native_rows': len(original),
        'matched_prior_rows': matched, 'missing_prior_rows': int(missing.sum()),
        'matched_free_prior_rows': int(np.sum((flags == free_value) & ~missing)),
        'matched_work_prior_rows': int(np.sum((flags != free_value) & ~missing)),
        'source_rows_not_in_native': len(set(lookup) - set(keys)),
        'prior_resolution_high': float(distances[~missing].min()),
        'missing_prior_resolution_high': float(distances[missing].min()) if missing.any() else None,
        'native_columns_bitwise_unchanged': True, 'native_row_order_unchanged': True,
        'native_column_labels': gemmi.read_mtz_file(str(native_path)).column_labels(),
        'missing_values_policy': 'NaN; never replaced by zero, model phases, extrapolation or nearest HKL',
        'matching': 'Exact integer HKL; same space-group setting and unit cell; no reindexing or symmetry expansion'}
    output.with_suffix('.audit.json').write_text(json.dumps(report, indent=2) + '\n')
    return report


def parse_statistics(log_path, xml_path):
    root = ET.parse(str(xml_path)).getroot()
    if int(root.findtext('err_level', '999')) != 0:
        raise RuntimeError('REFMAC XML reports error')
    sections = Path(log_path).read_text().split('Ncyc    Rfact    Rfree')
    if len(sections) < 2:
        raise RuntimeError('Missing REFMAC R-factor cycle table')
    cycles = []
    for line in sections[-1].split('$TEXT:')[0].splitlines():
        fields = line.split()
        if len(fields) < 11 or not fields[0].isdigit():
            continue
        try:
            values = [float(value) for value in fields[1:11]]
        except ValueError:
            continue
        if not np.isfinite(values).all() or min(values[:2]) < 0.:
            raise RuntimeError('Invalid REFMAC numeric cycle table')
        cycles.append({'cycle': int(fields[0]), 'r_work': values[0], 'r_free': values[1],
            'bond_rmsd_angstrom': values[5], 'angle_rmsd_degrees': values[7]})
    xml_cycles = root.findall('./Overall_stats/stats_vs_cycle/new_cycle')
    if not cycles or len(cycles) != len(xml_cycles):
        raise RuntimeError('REFMAC log/XML cycle counts disagree')
    for value, xml in zip(cycles, xml_cycles):
        if value['cycle'] != int(xml.findtext('cycle')):
            raise RuntimeError('REFMAC log/XML cycle numbers disagree')
        for key, name in [('r_work', 'r_factor'), ('r_free', 'r_free')]:
            if abs(value[key] - float(xml.findtext(name))) > .0011:
                raise RuntimeError('REFMAC log/XML R factors disagree')
    total = int(root.findtext('./Overall_stats/n_reflections_all'))
    free = int(root.findtext('./Overall_stats/n_reflections_free'))
    return {'initial': cycles[0], 'final': cycles[-1], 'cycles': cycles,
            'reflections_used_total': total, 'reflections_used_free': free,
            'reflections_used_work': total - free,
            'resolution_low': float(root.findtext('./Overall_stats/resolution_low')),
            'resolution_high': float(root.findtext('./Overall_stats/resolution_high')),
            'r_factor_precision_decimals': 4}


def refmac_with_experimental_prior(merged_mtz, model_path, output, use_prior=True,
                                   cycles=10, representation='hl', timeout=3600,
                                   ccp4=None, free_value=0):
    if representation not in ['hl', 'phifom'] or not isinstance(cycles, int) or cycles < 0:
        raise ValueError('Invalid phase representation or cycle count')
    native, partition = reflection_partition(merged_mtz, free_value=free_value)
    structure = gemmi.read_structure(str(model_path))
    if (len(structure) != 1 or not structure[0] or structure.find_spacegroup() is None or
            structure.find_spacegroup().hm != native.spacegroup.hm or
            not np.allclose(structure.cell.parameters, native.cell.parameters, rtol=0., atol=.01)):
        raise ValueError('Model and native data crystal settings differ')
    out = Path(output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    ccp4 = resolve_runtime(ccp4, required=('refmac5',))
    binary = ccp4 / 'bin/refmac5'
    original_input = fingerprint(merged_mtz)
    # Legacy CCP4/Fortran file-name buffers can truncate a long absolute HKLIN.
    # Stage exactly the same bytes and pass only short names in cwd. This is an
    # execution repair; no observations, priors, columns or keywords change.
    shutil.copyfile(merged_mtz, out / 'input.mtz')
    staged_input = fingerprint(out / 'input.mtz')
    if staged_input['sha256'] != original_input['sha256']:
        raise RuntimeError('REFMAC staged input differs from original merged MTZ')
    structure.write_pdb(str(out / 'input.pdb'))
    labels = ['HLA', 'HLB', 'HLC', 'HLD'] if representation == 'hl' else ['PHIB', 'FOM']
    if use_prior:
        if any(native.column_with_label('SAD_' + label) is None for label in labels):
            raise ValueError('Missing experimental prior columns')
        phase_assignment = ' ' + ' '.join(label + '=SAD_' + label for label in labels)
    else:
        phase_assignment = ''
    low, high = partition['observed_resolution_low'] * 1.0001, partition['observed_resolution_high'] * .9999
    keywords = ('LABIN FP=F SIGFP=SIGF FREE=FreeR_flag' + phase_assignment + '\nFREE %d\n' % free_value +
        'REFI TYPE REST RESO %.6f %.6f\n' % (low, high) +
        'REFI RESI MLKF\nREFI BREF ISOT\nMAKE HYDR N HOUT N BUIL N\n' +
        'SCAL TYPE BULK\nSCAL MLSC WORK\nWEIGHT AUTO\nNCYC %d\nEND\n' % cycles)
    (out / 'keywords.txt').write_text(keywords)
    command = [str(binary), 'HKLIN', 'input.mtz', 'XYZIN', 'input.pdb',
        'HKLOUT', 'refined.mtz', 'XYZOUT', 'refined.pdb', 'XMLOUT', 'refmac.xml']
    report = {'engine': 'external-REFMAC-experimental-phase-probe',
        'source': fingerprint(__file__), 'model_input': fingerprint(model_path),
        'input': original_input, 'staged_input': staged_input,
        'input_staging': 'Exact byte copy; short relative CCP4 file arguments to avoid legacy Fortran pathname truncation',
        'working_directory': str(out), 'binary': fingerprint(binary),
        'use_experimental_prior': bool(use_prior), 'prior_representation': representation if use_prior else None,
        'keywords': keywords, 'command': command, 'requested_cycles': cycles,
        'partition': partition, 'phase_blurring': 'REFMAC default SCBLUR1 BBLUR0; no tuned confidence',
        'selection': 'Fixed final cycle; no choice by Rfree or target-reference scores',
        'limitation': 'Original upstream OASIS/DM priors already used free reflections; this is not an independent blind Rfree test.'}
    environment = runtime_environment(ccp4, scratch=out / 'scratch')
    report['scratch_directory'] = environment['CCP4_SCR']
    started = time.time()
    try:
        with (out / 'refmac.log').open('w') as log:
            completed = subprocess.run(command, input=keywords, text=True, stdout=log,
                stderr=subprocess.STDOUT, cwd=str(out), env=environment, timeout=timeout)
        report['returncode'] = completed.returncode
        if completed.returncode != 0:
            raise RuntimeError('REFMAC failed; inspect ' + str(out / 'refmac.log'))
        phase_target = bool(re.search(r'Residual\s*:\s*Phase modified Maximum Likelihood for Fs',
                                     (out / 'refmac.log').read_text()))
        if phase_target != bool(use_prior):
            raise RuntimeError('REFMAC residual type does not match requested prior use')
        stats = parse_statistics(out / 'refmac.log', out / 'refmac.xml')
        stats['usable_observations_not_in_refmac_statistics'] = partition['usable_observations'] - stats['reflections_used_total']
        if stats['final']['cycle'] != cycles:
            raise RuntimeError('REFMAC did not finish requested cycles')
        flag_audit = verify_free_flags(native, gemmi.read_mtz_file(str(out / 'refined.mtz')),
                                      'FreeR_flag', 'F', 'SIGF')
        gemmi.read_structure(str(out / 'refined.pdb')).make_mmcif_document().write_file(str(out / 'model.cif'))
        report.update(status='success', statistics=stats, phase_target_verified=True,
                      phase_modified_likelihood_active=phase_target, free_flags_preserved=True,
                      free_flag_audit=flag_audit, output_model=fingerprint(out / 'model.cif'),
                      output_mtz=fingerprint(out / 'refined.mtz'))
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['seconds'] = time.time() - started
        (out / 'refinement_report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
