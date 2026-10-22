"""External-tool adapters; results must never be labelled as the native engine."""
import json
import subprocess
import time
from pathlib import Path
import gemmi
from .data import audit, fingerprint
from .runtime import resolve_runtime, runtime_environment


def buccaneer(mtz, output, ccp4=None, cycles=3, timeout=600):
    if cycles < 1:
        raise ValueError('Cycles must be positive')
    root, out = resolve_runtime(ccp4,required=('cbuccaneer',)), Path(output).resolve()
    binary = root / 'bin/cbuccaneer'
    library = root / 'lib/data/reference_structures/reference-1tqw'
    required = [binary, library.with_suffix('.pdb'), library.with_suffix('.mtz')]
    if not all(p.is_file() for p in required):
        raise ValueError('CCP4 Buccaneer binary or 1tqw training library not found')
    info = audit(mtz)
    labels = {c['label'] for c in info['columns']}
    if not {'F','SIGF','PHWT','FOM','FWT'}.issubset(labels):
        raise ValueError('This REFMAC-map adapter requires F,SIGF,PHWT,FOM,FWT in one MTZ')
    out.mkdir(parents=True, exist_ok=False)
    cmd = [str(binary), '-mtzin', str(Path(mtz).resolve()),
           '-mtzin-ref', str(library.with_suffix('.mtz')), '-pdbin-ref', str(library.with_suffix('.pdb')),
           '-colin-ref-fo', '/*/*/[FP.F_sigF.F,FP.F_sigF.sigF]',
           '-colin-ref-hl', '/*/*/[FC.ABCD.A,FC.ABCD.B,FC.ABCD.C,FC.ABCD.D]',
           '-colin-fo', '/*/*/[F,SIGF]', '-colin-phifom', '/*/*/[PHWT,FOM]',
           '-colin-fc', '/*/*/[FWT,PHWT]', '-pdbout', str(out / 'backbone.pdb'),
           '-xmlout', str(out / 'buccaneer.xml'), '-find', '-grow', '-join', '-link', '-prune',
           '-cycles', str(cycles), '-jobs', '1']
    started = time.perf_counter()
    report = {'engine': 'external-buccaneer', 'input': info['input'], 'command': cmd,
              'training_library': [fingerprint(p) for p in required[1:]],
              'binary': fingerprint(binary),
              'scope': 'Backbone-only CCP4 backend, not the autonomous PhaseBuild tracer. '
                       'Uses unrelated library 1tqw, no INS reference, no sequence or refinement.'}
    try:
        with (out / 'buccaneer.log').open('w') as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, timeout=timeout,env=runtime_environment(root))
        report['returncode'] = result.returncode
        if result.returncode:
            raise RuntimeError('Buccaneer failed; see ' + str(out / 'buccaneer.log'))
        pdb = out / 'backbone.pdb'
        if not pdb.exists():
            raise RuntimeError('Buccaneer did not create the output model')
        structure = gemmi.read_structure(str(pdb))
        structure.make_mmcif_document().write_file(str(out / 'backbone.cif'))
        report['residues'] = sum(1 for c in structure[0] for r in c if r.find_atom('CA','*') is not None)
        report['status'] = 'success'
    except Exception as error:
        report.update(status='failed', error=str(error))
        raise
    finally:
        report['total_seconds'] = time.perf_counter() - started
        (out / 'build_report.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
