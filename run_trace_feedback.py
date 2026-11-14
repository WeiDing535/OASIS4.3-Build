"""Execute the work-only density-modification + trace-feedback stage.

Reads the persisted pipeline configuration, uses only the prepared work data and
the persisted original observation likelihood, and writes the standard
``map_phases.mtz``/``map_phases.json`` contract consumed by the bridge stage.
"""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild import trace_feedback
from phasebuild.data import fingerprint


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    config_path=Path(a.config);cfg=json.loads(config_path.read_text())['config']
    root=config_path.parent
    out=Path(a.output)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError('Trace-feedback stage requires a new, empty directory: '+str(out))
    if cfg.get('auto_hand') and cfg['method']=='pplus_bundled':
        native=root/'02_sad/selected_native.mtz'
        sad=root/'02_sad/selected_sad.mtz'
    else:
        native=root/'01_inputs/native.mtz'
        sad=root/'01_inputs/sad.mtz'
    prior=root/'02_sad/observations.mtz'
    if not Path(native).is_file():
        # The phase stage always writes FP/PHIB/FOM to experimental.mtz, which
        # is all the density modifier needs; a phases-only run may have no
        # native amplitudes.
        from phasebuild.sad_bridge import FLAG
        out.mkdir(parents=True,exist_ok=True)
        native=Path(trace_feedback.work_from_experimental(root/'02_sad/experimental.mtz',sad,FLAG,
                                                         cfg['free_value'],out/'work_native.mtz'))
    for path in (native,prior):
        if not Path(path).is_file():raise ValueError('Missing stage input '+str(path))
    solvent=cfg.get('dm_solvent')
    if solvent is None:
        from phasebuild.trust_workflow import sample_info
        from phasebuild.sad_input import read_sad_pairs
        from phasebuild.sad_bridge import FLAG
        data,_,_=read_sad_pairs(str(root/'01_inputs/sad.mtz'),free_label=FLAG,free_value=cfg['free_value'])
        solvent=sample_info(cfg,data)['solvent_fraction']
    from phasebuild.runtime import resolve_runtime
    runtime=resolve_runtime(cfg.get('runtime') or None)
    report=trace_feedback.run(str(native),str(prior),cfg.get('sequence') or '',str(runtime),out,
        solvent_fraction=float(solvent),cycles=int(cfg.get('dm_cycles',20)),
        feedback_cycles=int(cfg.get('feedback_cycles',0)),
        trace_resolution=float(cfg.get('trace_resolution',2.0)),
        executable=cfg.get('dm_executable') or None,
        cmbwt=cfg.get('dm_cmbwt',0.4),
        holdout=float(cfg.get('dm_holdout',0.05)),
        hist_lib=cfg.get('dm_hist_lib') or None,
        phenix_bin=(cfg.get('phenix_bin') or None) or None,
        fcalc=cfg.get('fcalc','fc_all'),
        progress=lambda m:print(m,flush=True))
    report.update(solvent_source='user' if cfg.get('dm_solvent') is not None else 'sequence_estimate',
        original_prior=fingerprint(prior),
        map_formula='FWT=|FP|*FOM, PHWT=PHIB from the final handdm v2 HL quartet; work reflections only',
        free_reflections_in_map=0)
    (out/'map_phases.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
