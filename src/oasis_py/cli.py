# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""A single entry point for heavy-atom search and protein phasing."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import platform
import sys
import time
import gemmi
import numpy as np
import scipy
from . import __product__, __version__
from .data import read_mtz,read_atoms,parse_composition,write_mtz
from .phasing import phase_native
from .substructure import search_native,write_atoms,invert_hand
from .reference import run_reference
from .backends import phaser_sad
from .sad import phase_sad_likelihood
from .density import density_modify_parrot
from .native_density import density_modify_native
from .fom import FomCalibration
from .ensemble import phase_sad_ensemble,mix_mtz_files
from .reference_data import bundled_reference,reference_info,fetch_reference


def pairs(text):
    if not text:return {}
    result={}
    for item in text.split(','):
        key,value=item.split('=',1);result[key.strip()]=value.strip()
    return result


def provenance():
    return {'product':__product__,'oasis_py':__version__,'python':platform.python_version(),
            'numpy':np.__version__,'scipy':scipy.__version__,'gemmi':gemmi.__version__}


def add_data_options(parser):
    parser.add_argument('hklin',type=Path)
    parser.add_argument('--mode',choices=['sad','sir','dmr'],default='dmr')
    parser.add_argument('--labels',help='f=FP,sigf=SIGFP,dano=DANO,...')
    parser.add_argument('--dmin',type=float);parser.add_argument('--dmax',type=float)
    parser.add_argument('--free-label');parser.add_argument('--free-value',type=int)


def add_phase_options(parser):
    parser.add_argument('--composition',help='ASU atom counts, e.g. C:500,N:120,O:150,S:5')
    parser.add_argument('--residues',type=int,help='Approximate ASU composition from residue count')
    parser.add_argument('--model',type=Path);parser.add_argument('--heavy',type=Path)
    parser.add_argument('--fpp',help='Explicit f-double-prime, e.g. S=0.56,Se=3.8')
    parser.add_argument('--cycles',type=int,default=2,help='Internal phase cycles, not model-building cycles')
    parser.add_argument('--seed',type=int,default=1)
    parser.add_argument('--pseudo-fraction',type=float,default=.05)
    parser.add_argument('--aoe',type=float,default=.5)
    parser.add_argument('--rank-fit',action='store_true',help='Enable legacy-like phase-difference rank fitting')
    parser.add_argument('--sad-method',choices=['pplus','likelihood','mixture','alternating'],default='pplus',
                        help='Experimental continuous SAD, equal mixture or fixed-prior alternating DM; opt-in')
    parser.add_argument('--triplet-weight',type=float,default=.1,
                        help='Empirical triplet tempering for SAD likelihood method')
    parser.add_argument('--kmin',type=float,default=.03)
    parser.add_argument('--max-neighbors',type=int,default=1000,help='0 disables cap')
    parser.add_argument('--max-per-target',type=int,default=256,help='0 disables cap')
    parser.add_argument('--fixed',type=Path,help='MTZ with PHIB and FOM seeds in the same cell/space-group setting')
    parser.add_argument('--fixed-percent',type=float,default=100)


def fixed_arrays(path,data,percent):
    if not 0<percent<=100:raise ValueError('fixed-percent must lie in (0,100]')
    mtz=gemmi.read_mtz_file(str(path));mtz.ensure_asu()
    if mtz.spacegroup.operations()!=data.spacegroup.operations() or not np.allclose(mtz.cell.parameters,data.cell.parameters,atol=1e-3,rtol=1e-4):
        raise ValueError('Fixed phases must use matching cell and space-group setting')
    hkl=mtz.make_miller_array()
    phase=np.deg2rad(np.array(mtz.column_with_label('PHIB').array,dtype=float))
    fom=np.array(mtz.column_with_label('FOM').array,dtype=float)
    valid=np.isfinite(phase)&np.isfinite(fom)&(fom>=0)&(fom<=1)
    candidates=np.flatnonzero(valid)
    chosen=candidates[np.argsort(-fom[candidates],kind='stable')[:int(np.ceil(len(candidates)*percent/100))]]
    lookup={tuple(h):i for i,h in enumerate(data.hkl)}
    fp=np.zeros(len(data.f));fm=np.zeros(len(data.f));mask=np.zeros(len(data.f),bool)
    for j in chosen:
        i=lookup.get(tuple(hkl[j]))
        if i is not None:fp[i]=phase[j];fm[i]=fom[j];mask[i]=True
    if not mask.any():raise ValueError('No fixed phases match the working reflection set')
    return fp,fm,mask


def load(args):
    return read_mtz(args.hklin,args.mode,pairs(args.labels),args.dmin,args.dmax,
                    args.free_label,args.free_value)


def phase(args,data,model=None,heavy=None,density_options=None):
    composition=parse_composition(args.composition,args.residues)
    model=model if model is not None else (read_atoms(args.model,data) if args.model else None)
    heavy=heavy if heavy is not None else (read_atoms(args.heavy,data) if args.heavy else None)
    fpp={gemmi.Element(el).name:float(v) for el,v in pairs(args.fpp).items()}
    if args.sad_method in {'likelihood','mixture','alternating'}:
        if args.mode!='sad' or model is not None or args.fixed or args.rank_fit:
            raise ValueError('SAD likelihood requires SAD heavy atoms; protein models, fixed phases and rank fit are unsupported')
        if heavy is None:raise ValueError('SAD likelihood requires heavy atoms')
        if args.sad_method=='alternating':
            if density_options is None:raise ValueError('Alternating SAD requires phase --native-dm')
            from .alternating import density_modify_alternating
            first=phase_native(data,'sad',composition,heavy=heavy,fpp=fpp,cycles=args.cycles,
                aoe=args.aoe,kmin=args.kmin,max_neighbors=args.max_neighbors,max_per_target=args.max_per_target)
            second=phase_sad_likelihood(data,composition,heavy,fpp,cycles=args.cycles,
                kmin=args.kmin,max_neighbors=args.max_neighbors,max_per_target=args.max_per_target,
                triplet_weight=args.triplet_weight)
            result=density_modify_alternating(data,first,second,**density_options)
            result.metadata['initial_component_metadata']=[first.metadata,second.metadata]
            return result
        if args.sad_method=='mixture':
            return phase_sad_ensemble(data,composition,heavy,fpp,cycles=args.cycles,
                kmin=args.kmin,max_neighbors=args.max_neighbors,max_per_target=args.max_per_target,
                triplet_weight=args.triplet_weight,aoe=args.aoe,density_options=density_options)
        return phase_sad_likelihood(data,composition,heavy,fpp,cycles=args.cycles,
            kmin=args.kmin,max_neighbors=args.max_neighbors,max_per_target=args.max_per_target,
            triplet_weight=args.triplet_weight)
    fixed=fixed_arrays(args.fixed,data,args.fixed_percent) if args.fixed else (None,None,None)
    return phase_native(data,args.mode,composition,model,heavy,fpp,args.cycles,args.seed,
                        args.pseudo_fraction,args.aoe,args.rank_fit,args.kmin,
                        args.max_neighbors,args.max_per_target,*fixed)


def no_overwrite(paths):
    for path in paths:
        if Path(path).exists():raise ValueError(f'Output exists: {path}; choose a new output path')
        Path(path).parent.mkdir(parents=True,exist_ok=True)


def main(argv=None):
    parser=argparse.ArgumentParser(description='OASIS4.3-GUI-0.1 protein phasing. Results require scientific validation.')
    parser.add_argument('--version',action='version',version=__version__)
    subs=parser.add_subparsers(dest='command',required=True)
    inspect=subs.add_parser('inspect',help='Inspect MTZ without modifying it');inspect.add_argument('hklin',type=Path)
    subs.add_parser('reference-info',help='Show default generic reference provenance')
    fetch=subs.add_parser('reference-fetch',help='Provision pinned upstream reference data (one-time network access)')
    fetch.add_argument('--source',type=Path,help='Use an existing exact reference file instead of downloading')
    fetch.add_argument('--force',action='store_true',help='Replace a corrupt/outdated cached reference')
    ph=subs.add_parser('phase',help='Native experimental SAD/SIR/DMR phasing')
    add_data_options(ph);add_phase_options(ph);ph.add_argument('--out',type=Path,required=True)
    ph.add_argument('--native-dm',action='store_true',help='Experimental independent Python density modification after native phasing')
    ph.add_argument('--solvent-fraction',type=float)
    ph.add_argument('--sequence',type=Path,help='Protein FASTA/PIR; derive ASU composition and Matthews solvent estimate')
    ph.add_argument('--copies',type=int,default=1,help='Copies of the complete sequence file in one ASU')
    ph.add_argument('--reference-mtz',type=Path,help='Generic reference MTZ; defaults to the provisioned 1TQW reference cache')
    ph.add_argument('--dm-cycles',type=int,default=5)
    ph.add_argument('--dm-engine',choices=['legacy','mlhl-best'],default='legacy',
        help='Optional experimental Clipper MLHL kernel + best-coefficient feedback; SAD likelihood only')
    ph.add_argument('--dm-pattern',choices=['A','B','AB','BA'],
        help='Fixed P+/SAD prior schedule; requires --sad-method alternating --native-dm')
    ph.add_argument('--fom-calibration',type=Path,help='Excluded-target empirical model for final map weights; omits HL/PPLUS')
    ph.add_argument('--target-id',help='Required calibration provenance id')
    ph.add_argument('--dm-anisotropy',action='store_true')
    ph.add_argument('--dm-ncs-operators',type=Path,help='JSON list of Cartesian rotation/translation operators from allowed partial models')
    search=subs.add_parser('substructure',help='Experimental native sparse-atom search')
    add_data_options(search);search.set_defaults(mode='sad')
    search.add_argument('--sites',type=int,required=True);search.add_argument('--element',default='S')
    search.add_argument('--trials',type=int,default=8);search.add_argument('--search-cycles',type=int,default=50)
    search.add_argument('--seed',type=int,default=1);search.add_argument('--out',type=Path,required=True)
    solve=subs.add_parser('solve-sad',help='Unified native substructure search followed by SAD phasing; evaluate both hands')
    add_data_options(solve);add_phase_options(solve);solve.set_defaults(mode='sad')
    solve.add_argument('--sites',type=int,required=True);solve.add_argument('--element',default='S')
    solve.add_argument('--trials',type=int,default=8);solve.add_argument('--search-cycles',type=int,default=50)
    solve.add_argument('--out-dir',type=Path,required=True)
    ref=subs.add_parser('reference',help='Adapter to a supplied Fortran executable; not a rewrite')
    ref.add_argument('--executable',type=Path,required=True);ref.add_argument('--keywords',type=Path,required=True)
    ref.add_argument('--work-dir',type=Path,required=True);ref.add_argument('--hklin',type=Path)
    ref.add_argument('--model',type=Path);ref.add_argument('--library-dir',type=Path);ref.add_argument('--timeout',type=float,default=600)
    ref.add_argument('--program',choices=['oasis','sapi'],default='oasis')
    ml=subs.add_parser('phaser-sad',help='External maximum-likelihood backend; optional EP_SSD substructure search')
    ml.add_argument('hklin',type=Path);ml.add_argument('--sequence',type=Path,required=True)
    ml.add_argument('--phenix-python',type=Path,required=True);ml.add_argument('--out-dir',type=Path,required=True)
    ml.add_argument('--wavelength',type=float,required=True);ml.add_argument('--heavy',type=Path)
    ml.add_argument('--sites',type=int);ml.add_argument('--element',default='S');ml.add_argument('--copies',type=int,default=1)
    ml.add_argument('--jobs',type=int,default=1);ml.add_argument('--dmin',type=float)
    ml.add_argument('--timeout',type=float,default=1800)
    ml.add_argument('--labels',help='Four comma-separated Friedel-pair column labels')
    dm=subs.add_parser('density-modify',help='External CCP4 Parrot solvent flattening and histogram matching')
    dm.add_argument('hklin',type=Path)
    dm.add_argument('--executable',type=Path,required=True)
    dm.add_argument('--out-dir',type=Path,required=True)
    dm.add_argument('--solvent-fraction',type=float,required=True)
    dm.add_argument('--reference-mtz',type=Path,required=True)
    dm.add_argument('--reference-pdb',type=Path,required=True)
    dm.add_argument('--cycles',type=int,default=20)
    dm.add_argument('--representation',choices=['hl','phifom'],default='hl')
    dm.add_argument('--canonical-centric-hl',action='store_true',
        help='Equivalent centric HL representation to avoid older Clipper underflow; HL input only')
    dm.add_argument('--timeout',type=float,default=300)
    mix=subs.add_parser('mix-phases',help='Equal first-moment mixture of two matching PHIB/FOM MTZ outputs')
    mix.add_argument('first',type=Path);mix.add_argument('second',type=Path);mix.add_argument('--out',type=Path,required=True)
    args=parser.parse_args(argv)
    start=time.perf_counter()
    try:
        if args.command=='mix-phases':
            print(json.dumps(mix_mtz_files(args.first,args.second,args.out),indent=2));return 0
        if args.command=='reference-fetch':
            path=fetch_reference(args.source,args.force)
            print(json.dumps({'path':str(path),'sha256':reference_info()['mtz_sha256']},indent=2));return 0
        if args.command=='reference-info':
            print(json.dumps(reference_info(),indent=2,allow_nan=False));return 0
        if args.command=='inspect':
            m=gemmi.read_mtz_file(str(args.hklin))
            print(json.dumps({'spacegroup':m.spacegroup.hm,'cell':list(m.cell.parameters),
                              'reflections':m.nreflections,'columns':[(c.label,c.type) for c in m.columns]},indent=2));return 0
        if args.command=='reference':
            print(json.dumps(run_reference(args.executable,args.keywords,args.work_dir,args.hklin,
                                           model=args.model,library_dir=args.library_dir,timeout=args.timeout,program=args.program),indent=2));return 0
        if args.command=='phaser-sad':
            labels=args.labels.split(',') if args.labels else None
            if labels is not None and len(labels)!=4:raise ValueError('Four Friedel-pair labels required')
            result=phaser_sad(args.phenix_python,args.hklin,args.sequence,args.out_dir,args.wavelength,
                              args.heavy,args.sites,args.element,args.copies,labels,args.jobs,args.dmin,args.timeout)
            print(json.dumps(result,indent=2));return 0
        if args.command=='density-modify':
            result=density_modify_parrot(args.executable,args.hklin,args.out_dir,args.solvent_fraction,
                args.reference_mtz,args.reference_pdb,args.cycles,args.representation,args.timeout,
                canonical_centric=args.canonical_centric_hl)
            print(json.dumps(result,indent=2));return 0
        if args.command=='solve-sad' and args.sad_method=='alternating':
            raise ValueError('Alternating DM is currently available only via phase --native-dm')
        data,input_meta=load(args)
        if args.command=='phase':
            sidecar=args.out.with_suffix('.json');no_overwrite([args.out,sidecar])
            sequence_meta=None
            if args.sequence:
                if args.composition or args.residues:
                    raise ValueError('--sequence cannot be combined with --composition or --residues')
                from .sequence import sequence_asu
                sequence_meta=sequence_asu(args.sequence,args.copies,data.cell,data.spacegroup)
                args.composition=','.join(f'{k}:{v}' for k,v in sequence_meta['asu_composition'].items())
                if args.solvent_fraction is None:args.solvent_fraction=sequence_meta['solvent_fraction']
            elif args.copies!=1:
                raise ValueError('--copies requires --sequence')
            if args.native_dm and args.solvent_fraction is None:
                raise ValueError('--native-dm requires --solvent-fraction')
            if args.native_dm and args.reference_mtz is None:
                args.reference_mtz=bundled_reference()
            if (args.dm_anisotropy or args.dm_ncs_operators) and not args.native_dm:
                raise ValueError('Density modification options require --native-dm')
            if args.dm_engine=='mlhl-best':
                if not args.native_dm or args.mode!='sad' or args.sad_method!='likelihood':
                    raise ValueError('--dm-engine mlhl-best requires --native-dm --mode sad --sad-method likelihood')
                if args.dm_anisotropy or args.dm_ncs_operators or args.fom_calibration:
                    raise ValueError('MLHL-best has no validated anisotropy, NCS or FOM calibration combination')
            if args.fom_calibration and not args.target_id:raise ValueError('--fom-calibration requires --target-id')
            density_options=None
            if args.dm_pattern and args.sad_method!='alternating':
                raise ValueError('--dm-pattern requires --sad-method alternating')
            if args.sad_method=='alternating':
                if not args.native_dm:raise ValueError('--sad-method alternating requires --native-dm')
                if args.dm_anisotropy or args.dm_ncs_operators or args.fom_calibration:
                    raise ValueError('Alternating DM has no validated anisotropy, NCS or FOM calibration combination')
                density_options=dict(solvent_fraction=args.solvent_fraction,reference_mtz=args.reference_mtz,
                    cycles=args.dm_cycles,seed=args.seed,pattern=args.dm_pattern or 'AB')
            if args.native_dm and args.sad_method=='mixture':
                density_options=dict(solvent_fraction=args.solvent_fraction,reference_mtz=args.reference_mtz,
                    cycles=args.dm_cycles,seed=args.seed,anisotropy=args.dm_anisotropy,
                    ncs_operators=json.loads(args.dm_ncs_operators.read_text()) if args.dm_ncs_operators else None)
            result=phase(args,data,density_options=density_options)
            if args.native_dm and args.sad_method not in {'mixture','alternating'}:
                initial_meta=result.metadata
                if args.dm_engine=='mlhl-best':
                    from .density_mlhl import density_modify_mlhl
                    result=density_modify_mlhl(data,result,args.solvent_fraction,args.reference_mtz,
                        cycles=args.dm_cycles,seed=args.seed,feedback='best')
                else:
                    result=density_modify_native(data,result,args.solvent_fraction,args.reference_mtz,args.dm_cycles,args.seed,anisotropy=args.dm_anisotropy,ncs_operators=json.loads(args.dm_ncs_operators.read_text()) if args.dm_ncs_operators else None)
                result.metadata['initial_phasing']=initial_meta
            if args.fom_calibration:
                result=FomCalibration.read(args.fom_calibration).apply(data,result,args.target_id)
            write_mtz(args.out,data,result.phase,result.fom,result.hl,result.probability)
            meta={**result.metadata,'input':input_meta,'dependencies':provenance(),
                  'seconds':time.perf_counter()-start,'source':data.source}
            if sequence_meta is not None:
                meta['sequence_asu']=sequence_meta
                if args.native_dm:meta['sequence_asu']['solvent_fraction_used']=args.solvent_fraction
            if args.native_dm:
                import hashlib
                meta['reference_map']={'path':str(args.reference_mtz.resolve()),
                    'sha256':hashlib.sha256(args.reference_mtz.read_bytes()).hexdigest(),
                    'bundled':args.reference_mtz.resolve()==bundled_reference().resolve()}
            sidecar.write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
        elif args.command=='substructure':
            no_overwrite([args.out,args.out.with_suffix('.json')])
            result=search_native(data,args.sites,args.element,args.trials,args.search_cycles,args.seed)
            write_atoms(args.out,result.atoms,data)
            meta={**result.metadata,'proxy_correlation':result.correlation,'trials':result.trials,
                  'input':input_meta,'dependencies':provenance()}
            args.out.with_suffix('.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
        else:
            if args.mode!='sad':raise ValueError('solve-sad only supports SAD data')
            if not args.fpp:raise ValueError('Explicit --fpp is required before starting search')
            if args.out_dir.exists() and any(args.out_dir.iterdir()):raise ValueError('Use an empty output directory')
            args.out_dir.mkdir(parents=True,exist_ok=True)
            found=search_native(data,args.sites,args.element,args.trials,args.search_cycles,args.seed)
            other,sg=invert_hand(found.atoms,data.spacegroup)
            if sg.operations()!=data.spacegroup.operations():
                raise ValueError('Enantiomorphic space-group change requires consistent data reindexing; run hands separately')
            outputs=[]
            for name,atoms in [('hand1',found.atoms),('hand2',other)]:
                write_atoms(args.out_dir/f'{name}.pdb',atoms,data)
                result=phase(args,data,heavy=atoms)
                write_mtz(args.out_dir/f'{name}.mtz',data,result.phase,result.fom,result.hl,result.probability)
                (args.out_dir/f'{name}.json').write_text(json.dumps(result.metadata,indent=2,allow_nan=False)+'\n')
                outputs.append(name)
            meta={'status':'experimental','hands':outputs,'hand_selected':False,
                  'selection_requires':'density modification and model-building validation',
                  'substructure':found.metadata,'proxy_correlation':found.correlation,
                  'input':input_meta,'dependencies':provenance(),'seconds':time.perf_counter()-start}
            (args.out_dir/'run.json').write_text(json.dumps(meta,indent=2,allow_nan=False)+'\n')
        print(json.dumps(meta,indent=2,allow_nan=False));return 0
    except (ValueError,RuntimeError,OSError) as exc:
        print(f'Error: {exc}',file=sys.stderr);return 2

if __name__=='__main__':raise SystemExit(main())
