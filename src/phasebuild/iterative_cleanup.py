"""Optional exact-sequence/work-map cleanup before each refinement round.

The wrapper never reads Rfree, generates coordinates or supplies a target
reference. Existing sequence/copy/coverage gates are retained. Refusal leaves
that round's original model intact and is explicitly recorded.
"""
import json
from pathlib import Path
import gemmi
from .data import fingerprint
from .oasis_options import read_sequences
from .partial_model_cleanup import work_density,cleanup_partial_single_copy


def wrap_refiner(refiner,dm_map,native,sequence,*,asu_copies,free_value=0):
    if isinstance(asu_copies,bool) or not isinstance(asu_copies,int) or asu_copies!=1:
        raise ValueError('Per-round cleanup requires explicit single-copy ASU')
    records=read_sequences(sequence)
    if len(records)!=1:raise ValueError('Per-round cleanup requires one sequence record')
    grid,map_audit=work_density(dm_map,native,free_value)
    references={name:fingerprint(path) for name,path in [('dm_map',dm_map),('native',native),('sequence',sequence)]}
    def refine(mtz,model,output,*args,**kwargs):
        for item in references.values():
            if fingerprint(item['path'])!=item:raise RuntimeError('Per-round cleanup input changed during recycling')
        if fingerprint(mtz)['sha256']!=references['native']['sha256']:raise ValueError('Per-round cleanup native differs from declared input')
        out=Path(output);stage=out.with_name(out.name+'_sequence_cleanup');stage.mkdir(parents=True,exist_ok=False)
        source=fingerprint(model);selected=Path(model)
        try:
            clean,audit=cleanup_partial_single_copy(gemmi.read_structure(str(model)),records[0][1],grid,asu_copies=asu_copies)
        except ValueError as error:
            audit=dict(status='refused',reason=str(error),source_model=source,selection_uses_Rfree=False)
        else:
            selected=stage/'cleaned.cif';clean.make_mmcif_document().write_file(str(selected))
            audit.update(status='completed',source_model=source,selected_model=fingerprint(selected),working_map_audit=map_audit)
        audit.update(preprocessor_source=fingerprint(__file__),fixed_inputs={k:dict(v) for k,v in references.items()})
        (stage/'cleanup.json').write_text(json.dumps(audit,indent=2)+'\n')
        report=refiner(mtz,selected,output,*args,**kwargs)
        if fingerprint(model)!=source:raise RuntimeError('Original constructed model changed during cleanup/refinement')
        for item in references.values():
            if fingerprint(item['path'])!=item:raise RuntimeError('Per-round cleanup input changed during refinement')
        report['pre_refinement_cleanup']=audit
        (out/'refinement_report.json').write_text(json.dumps(report,indent=2)+'\n')
        return report
    return refine
