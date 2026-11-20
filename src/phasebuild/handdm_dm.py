"""Self-contained density-modification adapter (handdm), a dm/cparrot replacement.

This is an in-house engine that performs solvent flattening and histogram
matching without any CCP4 installation. It exposes the same call contract and
output schema as :mod:`phasebuild.parrot_dm`, so the pipeline can switch between
the two engines by name.

Input: a work-only MTZ (FP/SIGFP + FreeR_flag) plus the experimental HL quartet.
Output: an MTZ carrying F/SIGF/FWT/PHWT/PHIB/FOM/HLA..HLD, the schema the rest of
the pipeline consumes. The free set is never passed to the modifier.
"""
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import gemmi
import numpy as np

# Default locations for the standalone handdm deployment.
DEFAULT_SCRIPTS = (
    '/Volumes/Fast_4T/Workbuddy_Oasis/oasis-dm/src/handdm_v2.py',
)
DEFAULT_HIST_LIBS = (
    '/Volumes/Fast_4T/Workbuddy_Oasis/oasis-dm/hist.lib',
    '/Applications/ccp4-9/lib/data/hist.lib',
)


def locate(script=None, runtime=None):
    """Return a usable handdm_v2.py path, else raise with a clear message."""
    if script:
        path = Path(script).expanduser()
        if path.is_file():
            return str(path)
        raise ValueError('Configured handdm script is not a file: ' + str(path))
    if runtime:
        candidate = Path(runtime) / 'handdm_v2.py'
        if candidate.is_file():
            return str(candidate)
    for base in DEFAULT_SCRIPTS:
        if Path(base).is_file():
            return base
    raise ValueError('handdm engine was not found; set --handdm-script to handdm_v2.py')


def _hist_lib(explicit=None):
    if explicit:
        path = Path(explicit).expanduser()
        if path.is_file():
            return str(path)
        raise ValueError('Configured handdm hist.lib is not a file: ' + str(path))
    for base in DEFAULT_HIST_LIBS:
        if Path(base).is_file():
            return base
    return None


def _hl_columns(path):
    m = gemmi.read_mtz_file(str(path))
    labels = m.column_labels()
    for prefix in ('DM', 'PIRATE', 'PARROT'):
        quartet = [prefix + '.ABCD.' + k for k in 'ABCD']
        if all(c in labels for c in quartet):
            return quartet
    quartet = ['HL' + k for k in 'ABCD']
    if not all(c in labels for c in quartet):
        raise ValueError('No HL quartet in ' + str(path))
    return quartet


def _hl_to_phifom(values):
    """Von Mises centroid: HL (A,B) -> phase(deg)/FOM, matching parrot_dm."""
    from scipy.special import i0e, i1e
    values = np.asarray(values, dtype=float)
    kappa = np.hypot(values[:, 0], values[:, 1])
    phase = np.degrees(np.arctan2(values[:, 1], values[:, 0])) % 360.0
    clipped = np.clip(kappa, 0, 30)
    fom = np.where(kappa > 30, 1.0, i1e(clipped) / np.maximum(i0e(clipped), 1e-300))
    return phase, fom


def _stage_input(work_mtz, hl_mtz, destination):
    """Merge work amplitudes with the experimental HL quartet into one MTZ.

    Columns written: FOBS, SIGFOBS, PHIB, FOM, HLA..HLD, FreeR. Only work
    reflections shared with the HL prior are kept; free reflections are excluded.
    """
    work = gemmi.read_mtz_file(str(work_mtz))
    work.ensure_asu()
    labels = work.column_labels()
    fp_label = next((c for c in ('FP', 'F', 'FOBS') if c in labels), None)
    sig_label = next((c for c in ('SIGFP', 'SIGF', 'SIGFOBS') if c in labels), None)
    if fp_label is None or sig_label is None:
        raise ValueError('work MTZ lacks an amplitude/sigma pair: ' + str(work_mtz))
    free_label = 'FreeR_flag' if 'FreeR_flag' in labels else ('FreeR' if 'FreeR' in labels else None)
    hl_source = hl_mtz or work_mtz
    quartet = _hl_columns(hl_source)
    src = gemmi.read_mtz_file(str(hl_source))
    src.ensure_asu()
    src_labels = src.column_labels()
    lookup = {tuple(int(v) for v in row[:3]): row for row in np.array(src)}
    rows = []
    for row in np.array(work):
        h = tuple(int(v) for v in row[:3])
        other = lookup.get(h)
        if other is None:
            continue
        hl = [other[src_labels.index(c)] for c in quartet]
        phase, fom = _hl_to_phifom(np.asarray([hl]))
        free = row[labels.index(free_label)] if free_label else 0
        rows.append(list(h) + [row[labels.index(fp_label)], row[labels.index(sig_label)],
                               float(phase[0]), float(fom[0]), *hl, free])
    if not rows:
        raise ValueError('No work reflections shared with the HL prior')
    m = gemmi.Mtz(with_base=True)
    m.spacegroup = work.spacegroup
    m.set_cell_for_all(work.cell)
    m.add_dataset('handdm_work')
    m.add_column('FOBS', 'F')
    m.add_column('SIGFOBS', 'Q')
    m.add_column('PHIB', 'P')
    m.add_column('FOM', 'W')
    for k in 'ABCD':
        m.add_column('HL' + k, 'A')
    m.add_column('FreeR', 'I')
    m.set_data(np.asarray(rows, np.float32))
    m.write_to_file(str(destination))
    return len(rows), quartet


def _remap(raw, destination, free_value=0):
    """Rewrite a handdm output MTZ into the pipeline's F/SIGF/FWT/PHWT schema."""
    m = gemmi.read_mtz_file(str(raw))
    labels = m.column_labels()
    a = np.array(m, copy=True)

    def col(name):
        if name not in labels:
            raise ValueError('handdm output lacks column ' + name)
        return a[:, labels.index(name)]

    fo = col('FOBS')
    sig = col('SIGFOBS')
    fwt = col('DM_FWT')
    phwt = col('DM_PHWT')
    phib = col('DM_PHIDM')
    fom = col('DM_FOMDM')
    hla = col('DM_HLADM')
    hlb = col('DM_HLBDM')
    hlc = col('DM_HLCDM')
    hld = col('DM_HLDDM')
    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.set_cell_for_all(m.cell)
    out.add_dataset('handdm_dm')
    for label, kind in [('F', 'F'), ('SIGF', 'Q'), ('FWT', 'F'), ('PHWT', 'P'),
                        ('PHIB', 'P'), ('FOM', 'W'),
                        ('HLA', 'A'), ('HLB', 'A'), ('HLC', 'A'), ('HLD', 'A')]:
        out.add_column(label, kind)
    out.set_data(np.column_stack([a[:, :3], fo, sig, fwt, phwt, phib, fom,
                                  hla, hlb, hlc, hld]).astype(np.float32))
    out.write_to_file(str(destination))
    return {'mean_fom': float(np.nanmean(fom))}


def _free_rows(staged, free_value=0):
    """Count the reflections the staged input marks as free."""
    m = gemmi.read_mtz_file(str(staged))
    labels = m.column_labels()
    if 'FreeR' not in labels:
        return 0
    data = np.array(m)
    return int(np.sum(data[:, labels.index('FreeR')] == int(free_value)))


def run(work_mtz, hl_mtz, output_mtz, solvent_fraction, *, cycles=20, mask_radius=3.0,
        anisotropy=True, model_pdb=None, executable=None, runtime=None, metadata=None,
        scratch=None, timeout=3600, fcalc='fc', hist_lib=None, script=None,
        free_value=0):
    """Run handdm on work-only data and write the pipeline map schema.

    Mirrors :func:`phasebuild.parrot_dm.run`. ``model_pdb`` is accepted for
    signature compatibility but the in-house engine does not combine a partial
    model's Fcalc; the experimental HL prior is used directly and the flag is
    recorded so callers can tell the two behaviours apart.
    """
    if not 0 < solvent_fraction < 1:
        raise ValueError('Solvent fraction must lie between 0 and 1')
    script = locate(script, runtime)
    hist = _hist_lib(hist_lib)
    parent = Path(output_mtz).parent
    parent.mkdir(parents=True, exist_ok=True)
    staged = parent / 'handdm_input.mtz'
    work_reflections, quartet = _stage_input(work_mtz, hl_mtz, staged)
    raw = parent / 'handdm_raw.mtz'
    stats_path = parent / 'handdm_stats.json'
    log_path = parent / 'handdm.log'
    # Map Parrot's cycle budget onto handdm's iteration controls. The engine
    # also runs its own free-set early stop, so a generous cap is safe.
    total_cycles = max(1, int(cycles))
    # Prefer the interpreter running this pipeline: it already carries the
    # pinned gemmi/numpy/scipy stack that handdm needs.
    interpreter = os.environ.get('PHASEBUILD_PYTHON') or sys.executable or 'python3'
    command = [interpreter, str(script),
               '--mtzin', str(staged), '--mtzout', str(raw),
               '--fo', 'FOBS', '--sigf', 'SIGFOBS', '--phib', 'PHIB', '--fom', 'FOM',
               '--hl', 'HLA,HLB,HLC,HLD', '--free', 'FreeR',
               '--free-marker', str(int(free_value)),
               '--solc', '%.10g' % solvent_fraction,
               '--cycles', str(total_cycles), '--max-cycles', str(total_cycles),
               '--stats', str(stats_path), '--tag', parent.name]
    if hist:
        command += ['--hist-lib', hist]
    # The pipeline hands the engine work-only data, so no free reflections
    # survive. The engine's free-set early stop then has no criterion to
    # evaluate, and an older build would freeze its least converged iterate.
    # Turn the early stop off explicitly in that case.
    free_rows = _free_rows(staged, free_value)
    early_stop_note = None
    if free_rows == 0:
        command += ['--no-auto-stop']
        early_stop_note = ('no free reflections in the work-only input; the '
                            'free-set early stop was disabled and the full cycle '
                            'count was run')
    env = dict(os.environ)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               VECLIB_MAXIMUM_THREADS='1')
    if scratch:
        Path(scratch).mkdir(parents=True, exist_ok=True)
    with log_path.open('w') as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                   stdin=subprocess.DEVNULL, env=env, timeout=timeout)
    if completed.returncode or not raw.is_file():
        raise RuntimeError('handdm failed with exit %s; see %s'
                           % (completed.returncode, log_path))
    info = _remap(raw, output_mtz, free_value=free_value)
    stats = {}
    if stats_path.is_file():
        try:
            stats = json.loads(stats_path.read_text())
        except (OSError, ValueError):
            stats = {}
    return {'engine': 'handdm (in-house)', 'script': script, 'hist_lib': hist,
            'solvent_fraction': float(solvent_fraction),
            'mask_filter_radius': float(mask_radius),
            'anisotropy_correction': False,
            'cycles': total_cycles,
            'partial_model_used': False,
            'partial_model_note': ('in-house engine ignores the feedback model; '
                                   'experimental HL prior used directly'),
            'work_reflections': work_reflections, 'free_reflections_passed': 0,
            'free_reflections_in_input': free_rows,
            'early_stop_note': early_stop_note,
            'phase_source': 'original observation HL',
            'stats': stats.get('final', stats) if isinstance(stats, dict) else {},
            **(metadata or {}), **info}
