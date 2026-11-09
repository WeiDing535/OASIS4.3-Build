"""Unified SAD -> density modification -> model building, with audited restart."""
import argparse
import json
from pathlib import Path
import sys
import os
for _thread_key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS']:
    os.environ[_thread_key]='1'
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['sad','sequence','heavy','output']:p.add_argument('--'+name,required=True)
    p.add_argument('--native',default='',help='Optional original mean observations; blank uses selected SAD observations')
    p.add_argument('--native-labels',default='')
    p.add_argument('--auto-hand',action=argparse.BooleanOptionalAction,default=None,help='Resolve both hands by masked work-only three-cycle sequence tracing; no Rfree selection')
    p.add_argument('--refinement-prior',choices=['original_sad','dm_phifom'],default=None)
    p.add_argument('--per-round-cleanup',action=argparse.BooleanOptionalAction,default=False,help='Opt-in single-copy exact-sequence/work-map cleanup before each refinement; original model retained on refusal')
    p.add_argument('--canonical-dm-seed',action=argparse.BooleanOptionalAction,default=False,help='Opt-in map-only phase0.001deg/FOM1e-6 precision contract; original P+ results unchanged')
    p.add_argument('--phase-blur-scale',type=float,default=None,help='Opt-in DM phase restraint confidence scale, paired with --phase-blur-b')
    p.add_argument('--phase-blur-b',type=float,default=None,help='Opt-in resolution-dependent DM phase restraint blurring B')
    p.add_argument('--building-recipe',choices=['legacy','ccp4','pipeline'],default=None,help='pipeline runs the same CCP4 buccaneer_pipeline builder as the reference route')
    p.add_argument('--input-columns',action=argparse.BooleanOptionalAction,default=False)
    p.add_argument('--fpp',default='',help='Legacy explicit factors for reproducibility; new runs can use --elements and mandatory --wavelength')
    p.add_argument('--elements',default='',help='Heavy-atom element symbols, e.g. Se or Se,S; computes theoretical factors offline from explicit wavelength')
    p.add_argument('--auto-fpp',action=argparse.BooleanOptionalAction,default=False)
    p.add_argument('--reference',default='');p.add_argument('--runtime');p.add_argument('--labels',default='')
    p.add_argument('--build-model',action=argparse.BooleanOptionalAction,default=True,
                   help='Enable Buccaneer/REFMAC (default); --no-build-model stops after DM without external engines')
    p.add_argument('--method',choices=['likelihood','pplus','pplus_trust','pplus_internal','pplus_bundled'],default=None,help='New-run default: bundled P+ Wang8 with work-weight search; explicit historical methods retained for reproducibility')
    p.add_argument('--solvent-fraction',type=float)
    p.add_argument('--aoe',type=float,default=.5)
    p.add_argument('--triplet-weight',type=float,default=.1)
    p.add_argument('--kmin',type=float,default=.03)
    p.add_argument('--max-neighbors',type=int,default=600)
    p.add_argument('--max-per-target',type=int,default=96)
    p.add_argument('--rank-fit',action='store_true')
    p.add_argument('--observation-prior',default='',help='Original pure Phaser SAD likelihood MTZ for pplus_trust')
    p.add_argument('--composition-sequence',default='',help='Optional matching sequence with explicit (MSE) residues')
    p.add_argument('--phenix-bin',default='')
    p.add_argument('--wavelength',type=float,default=None)
    p.add_argument('--trust-kl',type=float,default=.0005)
    p.add_argument('--trust-cleanup',action=argparse.BooleanOptionalAction,default=True)
    p.add_argument('--trust-work-weight',action=argparse.BooleanOptionalAction,default=None,help='Work-only pilots choose post-cleanup restraint weight; enabled for new default runs')
    p.add_argument('--trust-wang-radius',type=float,default=None,help='New default profile uses8 Angstrom; explicit historical routes retain14')
    p.add_argument('--trust-stability',action=argparse.BooleanOptionalAction,default=False,help='Opt-in disjoint-relation stability shrinkage of genuine P+ updates')
    p.add_argument('--trust-partial-cleanup',action=argparse.BooleanOptionalAction,default=None,help='Work-map/sequence-anchored cleanup; enabled for new default runs')
    p.add_argument('--dm-engine',choices=['resolve','parrot'],default=None,help='Density-modification engine; parrot adds anisotropic scaling correction and mask smoothing')
    p.add_argument('--dm-solvent',type=float,default=None,help='Override the sequence-estimated solvent content for density modification')
    p.add_argument('--dm-mask-radius',type=float,default=3.0,help='Solvent-mask filter radius in Angstrom')
    p.add_argument('--dm-anisotropy',action=argparse.BooleanOptionalAction,default=True,help='Anisotropic scaling correction during density modification')
    p.add_argument('--feedback-cycles',type=int,default=0,help='Work-only trace -> model Fcalc -> density-modification feedback cycles')
    p.add_argument('--trace-resolution',type=float,default=2.0,help='Resolution used for the partial-model tracing step')
    p.add_argument('--dm-executable',default='',help='Explicit cparrot path; blank searches the runtime and PATH')
    p.add_argument('--fcalc',choices=['fc','fc_all','refined'],default='fc_all',help='Partial-model structure factors used by the trace feedback: FC/PHIC or FC_ALL/PHIC_ALL (bulk solvent included)')
    p.add_argument('--build-semet',action='store_true')
    p.add_argument('--copies',default='auto',help='ASU 中序列文件的份数：正整数，或 auto（默认；先读 RCSB FASTA 头的 Chains 注释，否则取 Matthews 溶剂含量 20–65%% 的最小拷贝数）')
    for name,default in [('phase-cycles',5),('dm-cycles',5),('rounds',5),('cycles',10),('free-value',0),('seed',20261023)]:
        p.add_argument('--'+name,type=int,default=default)
    p.add_argument('--feedback-endpoints',default='',help='逗号分隔的 feedback 轮数（例如 2,3）：每个值各跑一遍固定协议，再按工作集 R 选择；不读 R-free')
    p.add_argument('--dm-solvents',default='',help='与 --feedback-endpoints 组合的溶剂含量，逗号分隔，auto 表示按序列估算（例如 auto,0.5）')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--check',action='store_true',help='Read-only preflight and plan; create no output')
    args=vars(p.parse_args(argv));resume=args.pop('resume');check=args.pop('check')
    endpoints=[int(x) for x in str(args.pop('feedback_endpoints') or '').replace(' ','').split(',') if x]
    solvents=[]
    for token in str(args.pop('dm_solvents') or '').replace(' ','').split(','):
        if not token:continue
        if token=='auto':solvents.append(None)
        else:
            try:solvents.append(float(token))
            except ValueError:p.error('--dm-solvents 必须是 auto 或 0 到 1 之间的数值')
    if args['copies']!='auto':
        try:args['copies']=int(args['copies'])
        except ValueError:p.error('--copies 必须是正整数或 auto')
    new_default=args['method'] is None
    if new_default:args['method']='pplus_bundled'
    profile=args['method']=='pplus_bundled' and (new_default or bool(args['elements']) or args['auto_fpp'])
    if new_default:args['input_columns']=True
    if args['auto_hand'] is None:args['auto_hand']=profile and args['build_model']
    if args['refinement_prior'] is None:args['refinement_prior']='dm_phifom' if profile else 'original_sad'
    if args['building_recipe'] is None:args['building_recipe']='pipeline' if profile else 'legacy'
    if args['dm_engine'] is None:args['dm_engine']='parrot' if profile else 'resolve'
    if profile and not args['fpp']:args['auto_fpp']=True
    for flag in ['trust_work_weight','trust_partial_cleanup']:
        if args[flag] is None:args[flag]=profile and args['build_model']
    if args['trust_wang_radius'] is None:args['trust_wang_radius']=8. if profile else 14.
    from phasebuild.sad_pipeline import validate,stage_plan,run_pipeline
    cfg=validate(args)
    if check:
        if Path(cfg['output']).exists() and not resume:p.error('Output exists; choose a new directory')
        print(json.dumps(dict(config=cfg,stages=stage_plan(cfg),output_created=False),indent=2));return 0
    if not endpoints:
        report=run_pipeline(cfg,resume=resume)
        print(json.dumps({'status':report['status'],'primary_model':report['primary_model'],
                          'dm_phases':report['dm_phases']},indent=2));return 0
    # Fixed-endpoint portfolio: every endpoint is a predeclared protocol, and the
    # retained one is chosen by the final work-set R only.
    root=Path(cfg['output']);runs=[]
    grid=[(fb,solvent) for fb in endpoints for solvent in (solvents or [cfg.get('dm_solvent')])]
    for feedback,solvent in grid:
        sub=dict(cfg);sub['feedback_cycles']=feedback
        if solvent is not None:sub['dm_solvent']=solvent
        tag='endpoint_fb%d%s'%(feedback,'' if solvent is None else '_s%s'%('%.2f'%solvent))
        sub['output']=str(root/tag)
        print(json.dumps({'event':'endpoint','feedback':feedback,'solvent':solvent,
                          'output':sub['output']}),flush=True)
        runs.append((feedback,solvent,tag,run_pipeline(sub,resume=False)))
    best_feedback,best_solvent,best_tag,best=min(runs,key=lambda item: item[3]['primary_metric']['r_work'])
    selection={'rule':'lowest final work-set R across the predeclared endpoint grid; R-free never read',
        'endpoints':[{'feedback':fb,'solvent':sol,'status':rep['status'],
                      'primary_metric':rep.get('primary_metric'),'output':str(root/tag)}
                     for fb,sol,tag,rep in runs],
        'selected_feedback':best_feedback,'selected_solvent':best_solvent,'selected_metric':best['primary_metric']}
    (root/'endpoint_selection.json').write_text(json.dumps(selection,indent=2)+'\n')
    print(json.dumps({'status':best['status'],'selected_feedback':best_feedback,'selected_solvent':best_solvent,
                      'primary_metric':best['primary_metric'],'primary_model':best['primary_model'],
                      'dm_phases':best['dm_phases'],'selection':str(root/'endpoint_selection.json')},indent=2))
    return 0


if __name__=='__main__':raise SystemExit(main())
