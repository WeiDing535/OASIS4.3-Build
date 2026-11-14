"""Density modification with the in-tree handdm v2 engine.

The engine itself is 丁玮's ``handdm_v2.py``, vendored byte-identically under
``phasebuild/vendor/handdm``.  It is a self-contained maximum-likelihood/HL
density modifier: it never calls an external density-modification program, and
the only files it reads are the staged MTZ and a density-histogram library.

Two properties matter for the pipeline contract:

* **The crystallographic free set is never touched.**  handdm v2 needs a
  held-out reflection set for its SigmaQ/D estimate and for its early stopping.
  We give it a deterministic random subset of *work* reflections
  (``--holdout``, default 5%), so ``free_reflections_in_map`` stays 0 and
  ``rfree_used_for_selection`` stays false.  The holdout is a pure
  cross-validation device: those reflections are still in the final map.
* **The output schema is fixed by the pipeline**, so the bridge, the builder
  and the report code stay simple: F, SIGF, FWT, PHWT, PHIB, FOM and an
  HLA..HLD quartet.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

import gemmi
import numpy as np

#: Deterministic holdout fraction of WORK reflections used as handdm's
#: cross-validation set.  Never the crystallographic free set.
DEFAULT_HOLDOUT = 0.05
#: Candidate recombination weights for ``cmbwt='auto'``.
CMBWT_CANDIDATES = (0.3, 0.4, 0.5, 0.6, 0.8)
#: handdm runs at most this many cycles regardless of the requested count; its
#: own free-set early stopping normally stops well before.
MAX_CYCLES = 20
#: Cycles used per candidate during a ``cmbwt='auto'`` calibration sweep.
CAL_CYCLES = 3
CALIBRATION_SEED = 7717
DEFAULT_CMBWT = 0.4
MIN_CMBWT_SPREAD = 0.01

_IMPORT_PROBE = {}


def _vendor_dir():
    return Path(__file__).resolve().parent / 'vendor' / 'handdm'


def _probe(python):
    """True when *python* can import gemmi, numpy and scipy."""
    key = str(python)
    if key not in _IMPORT_PROBE:
        try:
            done = subprocess.run([str(python), '-c', 'import gemmi,numpy,scipy'],
                                  stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=120)
            _IMPORT_PROBE[key] = done.returncode == 0
        except (OSError, subprocess.SubprocessError):
            _IMPORT_PROBE[key] = False
    return _IMPORT_PROBE[key]


def locate(script=None, runtime=None):
    """Return ``(python, script)`` for the vendored handdm v2, else raise.

    Resolution order: an explicit script, ``PHASEBUILD_HANDDM``, the vendored
    copy.  The interpreter is the one running the pipeline when it can import
    gemmi/numpy/scipy, else ``PHASEBUILD_HANDDM_PYTHON``, else
    ``<runtime>/bin/python``.
    """
    if script:
        path = Path(script).expanduser()
        if not path.is_file():
            raise ValueError('Configured handdm script does not exist: ' + str(path))
    else:
        override = os.environ.get('PHASEBUILD_HANDDM')
        if override:
            path = Path(override).expanduser()
            if not path.is_file():
                raise ValueError('PHASEBUILD_HANDDM does not point at a file: ' + str(path))
        else:
            path = _vendor_dir() / 'handdm_v2.py'
            if not path.is_file():
                raise ValueError('The vendored handdm engine is missing: ' + str(path))
    candidates = []
    override_python = os.environ.get('PHASEBUILD_HANDDM_PYTHON')
    if override_python:
        candidates.append(override_python)
    candidates.append(sys.executable)
    if runtime:
        candidates.append(str(Path(runtime) / 'bin' / 'python'))
    for candidate in candidates:
        if candidate and Path(candidate).is_file() and _probe(candidate):
            return str(candidate), str(path)
    raise ValueError('No Python interpreter with gemmi/numpy/scipy was found for the handdm engine; '
                     'set PHASEBUILD_HANDDM_PYTHON or point --dm-executable at a handdm_v2.py')


def hist_library(explicit=None):
    """CCP4 density-histogram library used as handdm's histogram-match target."""
    for candidate in (explicit, os.environ.get('PHASEBUILD_HIST_LIB'),
                      Path(__file__).resolve().parents[2] / 'lib' / 'data' / 'hist.lib',
                      '/Applications/ccp4-9/lib/data/hist.lib'):
        if candidate and Path(candidate).is_file():
            return str(candidate)
    return None


def _hl_columns(path):
    m = gemmi.read_mtz_file(str(path))
    labels = m.column_labels()
    for prefix in ('DM', 'PIRATE'):
        quartet = [prefix + '.ABCD.' + k for k in 'ABCD']
        if all(c in labels for c in quartet):
            return quartet
    quartet = ['HL' + k for k in 'ABCD']
    if not all(c in labels for c in quartet):
        raise ValueError('No HL quartet in ' + str(path))
    return quartet


def phi_fom(hl):
    """Centroid phase (degrees) and FOM of an HLA..HLD quartet."""
    from scipy.special import i0e, i1e
    a, b = np.asarray(hl)[:, 0], np.asarray(hl)[:, 1]
    kappa = np.hypot(a, b)
    fom = np.where(kappa > 0, i1e(kappa) / np.maximum(i0e(kappa), 1e-300), 0.0)
    phase = np.degrees(np.arctan2(b, a))
    return phase, np.clip(fom, 0.0, 0.9999)


def _holdout_flags(keys, fraction, seed):
    """Deterministic 0/1 vector: 0 = holdout (CCP4 free convention), 1 = work."""
    n = len(keys)
    if fraction <= 0 or n < 40:
        return np.ones(n, np.int32)
    count = max(20, int(round(n * fraction)))
    count = min(count, n - 20)
    rng = np.random.default_rng(seed)
    order = np.argsort(rng.random(n))[:count]
    flags = np.ones(n, np.int32)
    flags[order] = 0
    return flags


def stage_input(work_mtz, hl_mtz, destination, *, holdout=DEFAULT_HOLDOUT, seed=20261023):
    """Write FP/SIGFP + HLA..HLD + DM_HOLDOUT for one handdm pass.

    Only work reflections are staged: ``work_mtz`` is the merged work-only file
    the trace-feedback stage builds, so the crystallographic free set cannot
    reach the engine even by accident.
    """
    work = gemmi.read_mtz_file(str(work_mtz))
    work.ensure_asu()
    labels = work.column_labels()
    fp_label = 'FP' if 'FP' in labels else 'F'
    sig_label = 'SIGFP' if 'SIGFP' in labels else 'SIGF'
    src = gemmi.read_mtz_file(str(hl_mtz or work_mtz))
    src.ensure_asu()
    src_labels = src.column_labels()
    quartet = _hl_columns(hl_mtz or work_mtz)
    lookup = {tuple(int(v) for v in row[:3]): row for row in np.array(src)}
    rows = []
    for row in np.array(work):
        h = tuple(int(v) for v in row[:3])
        other = lookup.get(h)
        if other is None:
            continue
        values = [other[src_labels.index(c)] for c in quartet]
        if not np.isfinite(values).all():
            continue
        rows.append([*h, row[labels.index(fp_label)], row[labels.index(sig_label)], *values])
    if not rows:
        raise ValueError('No work reflections shared with the HL prior')
    data = np.asarray(rows, np.float32)
    flags = _holdout_flags([tuple(int(v) for v in row[:3]) for row in data], holdout, seed)
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = work.spacegroup
    m.set_cell_for_all(work.cell)
    m.add_dataset('handdm_work')
    phase, fom = phi_fom(data[:, 3:7])
    m.add_column('FP', 'F')
    m.add_column('SIGFP', 'Q')
    m.add_column('PHIB', 'P')
    m.add_column('FOM', 'W')
    for k in 'ABCD':
        m.add_column('HL' + k, 'A')
    m.add_column('DM_HOLDOUT', 'I')
    m.set_data(np.column_stack([data[:, :3], data[:, 3:5], phase, fom, data[:, 5:9],
                                flags]).astype(np.float32))
    m.write_to_file(str(destination))
    return {'work_reflections': len(rows), 'holdout_reflections': int((flags == 0).sum())}


def _shell_bins(s2, n_shells):
    if len(s2) < n_shells * 5:
        return np.zeros(len(s2), int), 1
    edges = np.quantile(s2, np.linspace(0, 1, n_shells + 1)[1:-1])
    return np.searchsorted(edges, s2), n_shells


def combine_model_hl(staged, model_pdb, destination, *, phenix_bin=None, runtime=None,
                     fcalc='fc_all', x_max=12.0, n_shells=10):
    """Add a work-only partial model's Fcalc to the experimental HL quartet.

    handdm has no option for an external Fcalc, so the equivalent likelihood is
    built here: per resolution shell
    the least-squares scale ``D = <Fo Fc>/<Fc^2>`` gives the Sim weight
    ``X = 2 Fo D Fc / (1 - D^2 Fc^2/Fo^2)``.  ``X`` is in observed-amplitude
    units while the quartet handdm consumes is normalised, so each shell's
    contribution is rescaled to the local experimental ``|HL|`` and multiplied
    by the fraction of scattered intensity the partial model explains
    (``<Fc^2>/<Fo^2>``, 0 for an empty trace, ~1 for a complete model).  Without
    that step a 12-residue fragment outweighs the experimental likelihood by
    7x and the map collapses.  Observed amplitudes are untouched.
    """
    from .fcalc import attach_fc
    attach_fc(staged, model_pdb, destination, phenix_bin, runtime, fcalc)
    m = gemmi.read_mtz_file(str(destination))
    labels = m.column_labels()
    data = np.array(m, copy=True)
    fp = data[:, labels.index('FP')].astype(float)
    fc_label = 'FC_ALL' if 'FC_ALL' in labels else ('FC' if 'FC' in labels else None)
    if fc_label is None:
        raise ValueError('The Fcalc pass produced no FC/FC_ALL column')
    fc = data[:, labels.index(fc_label)].astype(float)
    phi = data[:, labels.index('PHIC_ALL' if fc_label == 'FC_ALL' and 'PHIC_ALL' in labels else 'PHIC')]
    hkl = np.column_stack([data[:, labels.index(c)] for c in ('H', 'K', 'L')]).astype(int)
    d = np.asarray(m.cell.calculate_d_array(hkl), dtype=float)
    s2 = 1.0 / np.maximum(d, 1e-6) ** 2
    shell, used = _shell_bins(s2, n_shells)
    amp_exp = np.hypot(data[:, labels.index('HLA')].astype(float),
                       data[:, labels.index('HLB')].astype(float))
    hl = np.zeros((len(data), 2))
    for index in range(used):
        mask = shell == index
        if mask.sum() < 10:
            continue
        denom = float(np.sum(fc[mask] ** 2))
        d = float(np.sum(fp[mask] * fc[mask]) / denom) if denom > 0 else 0.0
        d = min(max(d, 0.0), 1.0)
        q = 1.0 - (d * fc[mask] / np.maximum(fp[mask], 1e-6)) ** 2
        x = 2.0 * fp[mask] * d * fc[mask] / np.maximum(q, 0.05)
        x = np.clip(x, 0.0, x_max)
        ax = x * np.cos(np.radians(phi[mask]))
        bx = x * np.sin(np.radians(phi[mask]))
        # The experimental quartet is a *normalised* HL (|HL| ~ 1-2), while the
        # Sim weight above is in observed-amplitude units, so adding it raw makes
        # a 12-residue fragment outweigh six sulfur sites by an order of
        # magnitude and destroys the map.  Rescale each shell to the local
        # experimental amplitude, weighted by the fraction of the scattered
        # intensity the partial model actually explains (0 = empty, ~1 = complete).
        rms_model = float(np.sqrt(np.mean(ax ** 2 + bx ** 2))) if len(ax) else 0.0
        rms_exp = float(np.sqrt(np.mean(amp_exp[mask] ** 2)))
        explained = float(np.sum(fc[mask] ** 2) / max(np.sum(fp[mask] ** 2), 1e-9))
        explained = min(max(explained, 0.0), 1.0)
        scale = explained * rms_exp / rms_model if rms_model > 1e-9 else 0.0
        hl[mask, 0] = scale * ax
        hl[mask, 1] = scale * bx
    ia, ib = labels.index('HLA'), labels.index('HLB')
    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.set_cell_for_all(m.cell)
    out.add_dataset('handdm_work')
    keep = ['FP', 'SIGFP', 'PHIB', 'FOM', 'HLA', 'HLB', 'HLC', 'HLD', 'DM_HOLDOUT']
    kinds = {'FP': 'F', 'SIGFP': 'Q', 'PHIB': 'P', 'FOM': 'W', 'HLA': 'A', 'HLB': 'A',
             'HLC': 'A', 'HLD': 'A', 'DM_HOLDOUT': 'I'}
    for label in keep:
        out.add_column(label, kinds[label])
    base = [data[:, labels.index(label)] for label in keep]
    quartet = np.column_stack([data[:, labels.index('HL' + k)] for k in 'ABCD']) + \
        np.column_stack([hl[:, 0], hl[:, 1], np.zeros(len(hl)), np.zeros(len(hl))])
    base[2], base[3] = phi_fom(quartet)
    base[4], base[5] = quartet[:, 0], quartet[:, 1]
    out.set_data(np.column_stack([hkl.astype(float)] + base).astype(np.float32))
    out.write_to_file(str(destination))
    return {'model_pdb': str(model_pdb), 'fcalc': fcalc,
            'mean_added_hl': float(np.hypot(hl[:, 0], hl[:, 1]).mean())}


def write_schema(dm_mtz, destination):
    """Rename handdm's DM_* columns to the pipeline's 03_dm schema."""
    m = gemmi.read_mtz_file(str(dm_mtz))
    labels = m.column_labels()
    data = np.array(m, copy=True)
    quartet = ['DM_HLADM', 'DM_HLBDM', 'DM_HLCDM', 'DM_HLDDM']
    missing = [c for c in quartet + ['DM_FWT', 'DM_PHWT', 'DM_PHIDM', 'DM_FOMDM'] if c not in labels]
    if missing:
        raise ValueError('handdm output is missing columns: ' + ', '.join(missing))
    fp = 'FP' if 'FP' in labels else 'F'
    sig = 'SIGFP' if 'SIGFP' in labels else 'SIGF'
    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.set_cell_for_all(m.cell)
    out.add_dataset('handdm_dm')
    for label, kind in [('F', 'F'), ('SIGF', 'Q'), ('FWT', 'F'), ('PHWT', 'P'),
                        ('PHIB', 'P'), ('FOM', 'W'), ('HLA', 'A'), ('HLB', 'A'),
                        ('HLC', 'A'), ('HLD', 'A')]:
        out.add_column(label, kind)
    columns = [data[:, labels.index(c)] for c in ('H', 'K', 'L')]
    columns += [data[:, labels.index(fp)], data[:, labels.index(sig)],
               data[:, labels.index('DM_FWT')], data[:, labels.index('DM_PHWT')],
               data[:, labels.index('DM_PHIDM')], data[:, labels.index('DM_FOMDM')],
               *[data[:, labels.index(c)] for c in quartet]]
    out.set_data(np.column_stack(columns).astype(np.float32))
    out.write_to_file(str(destination))
    return {'mean_fom': float(np.nanmean(columns[8])),
            'columns': [c.replace('DM_', '').replace('DM', '') for c in quartet]}


def _invoke(python, script, staged, destination, *, solvent_fraction, cycles,
            cmbwt, hist_lib, auto_stop=True, quiet=True, log_path=None, timeout=3600):
    command = [python, '-u', script, '--mtzin', str(staged), '--mtzout', str(destination),
               '--fo', 'FP', '--sigf', 'SIGFP',
               '--hl', 'HLA,HLB,HLC,HLD', '--free', 'DM_HOLDOUT', '--free-marker', '0',
               '--solc', '%.10g' % solvent_fraction,
               '--max-cycles', str(int(cycles)), '--min-cycles', '2',
               '--cmbwt', '%.10g' % float(cmbwt)]
    command += ['--auto-stop'] if auto_stop else ['--no-auto-stop']
    if quiet:
        command.append('--quiet')
    if hist_lib:
        command += ['--hist-lib', hist_lib]
    else:
        command.append('--no-hist')
    stats = Path(str(destination)).with_suffix('.stats.json')
    command += ['--stats', str(stats)]
    log_path = Path(log_path or Path(str(destination)).with_suffix('.log'))
    started = time.time()
    with log_path.open('w') as log:
        done = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                              stdin=subprocess.DEVNULL, timeout=timeout)
    if done.returncode or not Path(destination).is_file():
        raise RuntimeError('handdm failed with exit %s; see %s' % (done.returncode, log_path))
    payload = {}
    if stats.is_file():
        try:
            payload = json.loads(stats.read_text())
        except (OSError, ValueError):
            payload = {}
    payload.setdefault('final', {})
    payload['final']['elapsed_sec'] = round(time.time() - started, 2)
    return payload


def _score(payload):
    """Holdout amplitude CC reached by a calibration pass (higher is better)."""
    cycles = payload.get('cycles') or []
    values = [c.get('free_cc') for c in cycles if c.get('free_cc') is not None]
    values = [v for v in values if np.isfinite(v)]
    return float(max(values)) if values else float('nan')


def calibrate_cmbwt(work_mtz, hl_mtz, workdir, solvent_fraction, *,
                    candidates=CMBWT_CANDIDATES, cycles=CAL_CYCLES, holdout=DEFAULT_HOLDOUT,
                    hist_lib=None, executable=None, runtime=None, patience=0.01, log=print):
    """Public entry point: stage the data, then sweep ``cmbwt``.

    The calibration holdout uses a different seed from the production pass, so
    the weight is not tuned on the reflections that later score the real run.
    """
    python, script = locate(executable, runtime)
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    staged = workdir / 'handdm_cal_input.mtz'
    stage_input(work_mtz, hl_mtz, staged, holdout=holdout, seed=CALIBRATION_SEED)
    return _calibrate_staged(python, script, staged, workdir, solvent_fraction=solvent_fraction,
                             hist_lib=hist_library(hist_lib),
                             candidates=candidates, cycles=cycles, patience=patience, log=log)


def _calibrate_staged(python, script, staged, workdir, *, solvent_fraction,
                      hist_lib, candidates=CMBWT_CANDIDATES, cycles=CAL_CYCLES,
                      patience=0.01, log=print):
    """Pick ``cmbwt`` by holdout cross-validation before the production pass.

    ``cmbwt`` damps the DM phase likelihood that is recombined with the original
    experimental HL.  Too large inflates FOM, too small throws the modified map
    away.  The optimum is data dependent, so it is measured rather than
    assumed: each candidate runs a short pass and is scored by the amplitude CC
    between |Fo| and the map's Fcalc on the *holdout* reflections, which are
    never used to build that map.  Ties are broken towards the smaller weight,
    because over-trusting the DM phases is the failure mode that hurts REFMAC.

    Returns ``(best, table)`` where *table* is a list of
    ``{'cmbwt', 'holdout_cc', 'cycles_run'}`` rows.
    """
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    table = []
    for weight in candidates:
        out = workdir / ('cal_%02d.mtz' % int(round(weight * 100)))
        payload = _invoke(python, script, staged, out, solvent_fraction=solvent_fraction,
                          cycles=cycles, cmbwt=weight,
                          hist_lib=hist_lib, log_path=workdir / ('cal_%02d.log' % int(round(weight * 100))))
        score = _score(payload)
        table.append({'cmbwt': float(weight), 'holdout_cc': score,
                      'cycles_run': int(payload.get('final', {}).get('n_cycles_run', 0))})
        log('cmbwt %.2f -> holdout CC %.4f' % (weight, score))
    usable = [row for row in table if np.isfinite(row['holdout_cc'])]
    if not usable:
        log('cmbwt calibration produced no finite score; keeping 0.4')
        return DEFAULT_CMBWT, table
    best = max(row['holdout_cc'] for row in usable)
    worst = min(row['holdout_cc'] for row in usable)
    # cmbwt only rescales the HL handed to the next stage, so it barely moves the
    # map that produces this score: a flat table means the criterion cannot see
    # the parameter, not that every weight is equally good.  Only override the
    # calibrated default when the sweep actually separates the candidates.
    if best - worst < MIN_CMBWT_SPREAD:
        log('cmbwt sweep is flat (%.4f..%.4f, spread %.4f < %.2f); keeping %.2f'
            % (worst, best, best - worst, MIN_CMBWT_SPREAD, DEFAULT_CMBWT))
        return DEFAULT_CMBWT, table
    for row in sorted(usable, key=lambda r: r['cmbwt']):
        if row['holdout_cc'] >= best - patience:
            return row['cmbwt'], table
    return usable[0]['cmbwt'], table


def run(work_mtz, hl_mtz, output_mtz, solvent_fraction, *, cycles=20,
        model_pdb=None, executable=None, runtime=None, metadata=None,
        scratch=None, timeout=3600, phenix_bin=None, fcalc='fc_all',
        cmbwt=0.4, holdout=DEFAULT_HOLDOUT, holdout_seed=20261023, hist_lib=None,
        progress=None):
    """Run handdm v2 on the work-only data and write the pipeline map schema.

    ``work_mtz`` must contain FP/SIGFP and FreeR_flag; free reflections are never
    staged.  ``model_pdb`` is an optional work-only partial model whose Fcalc is
    combined with the fixed experimental likelihood (the feedback loop).
    """
    if not 0 < solvent_fraction < 1:
        raise ValueError('Solvent fraction must lie between 0 and 1')
    python, script = locate(executable, runtime)
    out = Path(output_mtz)
    out.parent.mkdir(parents=True, exist_ok=True)
    staged = out.parent / 'handdm_input.mtz'
    staged_info = stage_input(work_mtz, hl_mtz, staged, holdout=holdout, seed=holdout_seed)
    model_info = None
    if model_pdb is not None:
        if not Path(model_pdb).is_file():
            raise ValueError('Feedback model is missing')
        combined = out.parent / 'handdm_with_model.mtz'
        model_info = combine_model_hl(staged, model_pdb, combined, phenix_bin=phenix_bin,
                                      runtime=runtime, fcalc=fcalc)
        staged = combined
    hist = hist_library(hist_lib)
    budget = max(2, min(int(cycles), MAX_CYCLES))
    calibration = None
    weight = cmbwt
    if isinstance(cmbwt, str) and cmbwt == 'auto':
        weight, calibration = _calibrate_staged(
            python, script, staged, out.parent / 'cmbwt_calibration',
            solvent_fraction=solvent_fraction, hist_lib=hist,
            log=progress or (lambda message: None))
    raw = out.parent / 'handdm_raw.mtz'
    payload = _invoke(python, script, staged, raw, solvent_fraction=solvent_fraction,
                      cycles=budget, cmbwt=weight,
                      hist_lib=hist, log_path=out.parent / 'handdm.log', timeout=timeout)
    info = write_schema(raw, out)
    final = payload.get('final', {})
    return {'engine': 'handdm v2 (in-tree)', 'executable': script, 'interpreter': python,
            'solvent_fraction': float(solvent_fraction),
            'cycles_requested': int(cycles), 'cycles_run': int(final.get('n_cycles_run', 0)),
            'best_cycle': int(final.get('best_cycle', 0)),
            'partial_model_used': model_pdb is not None,
            'phase_source': ('original observation HL' if model_pdb is None
                             else 'original HL combined with work-only model Fcalc'),
            'holdout_fraction': float(holdout),
            'holdout_reflections': staged_info['holdout_reflections'],
            'work_reflections': staged_info['work_reflections'],
            'free_reflections_passed': 0,
            'rfree_used_for_selection': False,
            'cmbwt': float(weight), 'cmbwt_auto': calibration is not None,
            'cmbwt_table': calibration,
            'best_free_cc': final.get('best_free_cc'),
            'elapsed_sec': final.get('elapsed_sec'),
            'sigmaq_source': final.get('sigmaq_source'),
            'hist_library': hist,
            **(model_info or {}), **(metadata or {}), **info}
