"""Audited external REFMAC refinement with an existing, unchanged free set."""
import json
import os
import re
import subprocess
import time
import xml.etree.ElementTree as ET
from pathlib import Path
import gemmi
import numpy as np
from . import __version__
from .atom_validation import atom_geometry_report
from .chemistry import monomer
from .data import fingerprint
from .runtime import resolve_runtime, runtime_environment


def reflection_partition(path,fp='F',sigfp='SIGF',free='FreeR_flag',free_value=0):
    mtz=gemmi.read_mtz_file(str(path))
    expected={fp:'F',sigfp:'Q',free:'I'}
    for label,kind in expected.items():
        if mtz.column_labels().count(label)>1:
            raise ValueError('Ambiguous duplicate MTZ column: '+label)
        column=mtz.column_with_label(label)
        if column is None or column.type!=kind:
            raise ValueError('Required MTZ column/type missing: '+label+' / '+kind)
    f=np.asarray(mtz.column_with_label(fp));sigma=np.asarray(mtz.column_with_label(sigfp))
    flags=np.asarray(mtz.column_with_label(free))
    observed=np.isfinite(f)&np.isfinite(sigma)&(f>=0)&(sigma>0)
    if np.any(~np.isfinite(flags[observed])) or np.any(flags[observed]!=np.rint(flags[observed])):
        raise ValueError('Every usable observation must have an integer free-set flag')
    test=observed&(flags==free_value);work=observed&(flags!=free_value)
    if not test.any() or not work.any():
        raise ValueError('Existing free and work sets must both contain usable observations')
    values,counts=np.unique(flags[np.isfinite(flags)],return_counts=True)
    distances=mtz.make_d_array()[observed]
    return mtz,{'input':fingerprint(path),'fp':fp,'sigfp':sigfp,'free_column':free,
                'free_value':int(free_value),'flag_counts':{str(int(v)):int(n) for v,n in zip(values,counts)},
                'mtz_rows':mtz.nreflections,'usable_observations':int(observed.sum()),
                'usable_work':int(work.sum()),'usable_free':int(test.sum()),
                'resolution_high':mtz.resolution_high(),
                'observed_resolution_high':float(distances.min()),
                'observed_resolution_low':float(distances.max()),'flags_regenerated':False}


def parse_refmac(log_text,xml_path):
    root=ET.parse(str(xml_path)).getroot()
    if int(root.findtext('err_level','999'))!=0:
        raise RuntimeError('REFMAC XML reports an error')
    # The log cycle table has four R-factor decimals; XML has only three.
    sections=log_text.split('Ncyc    Rfact    Rfree')
    if len(sections)<2:
        raise RuntimeError('REFMAC cycle table missing; refusing to guess R factors')
    cycles=[]
    for line in sections[-1].split('$TEXT:')[0].splitlines():
        fields=line.split()
        if len(fields)<11 or not fields[0].isdigit():
            continue
        try:
            values=[float(x) for x in fields[1:11]]
        except ValueError:
            continue
        if not all(np.isfinite(values)) or not all(0.<=v<=1. for v in values[:2]):
            raise RuntimeError('Invalid R-factor row')
        cycles.append({'cycle':int(fields[0]),'r_work':values[0],'r_free':values[1],
                       'bond_rmsd_angstrom':values[5],'angle_rmsd_degrees':values[7]})
    xml_cycles=root.findall('./Overall_stats/stats_vs_cycle/new_cycle')
    if not cycles or len(cycles)!=len(xml_cycles):
        raise RuntimeError('REFMAC log/XML cycle counts disagree')
    for row,x in zip(cycles,xml_cycles):
        if row['cycle']!=int(x.findtext('cycle')) or any(
            abs(row[k]-float(x.findtext(tag)))>.0011 for k,tag in [('r_work','r_factor'),('r_free','r_free')]):
            raise RuntimeError('REFMAC log/XML values disagree')
    all_count=int(root.findtext('./Overall_stats/n_reflections_all'))
    free_count=int(root.findtext('./Overall_stats/n_reflections_free'))
    return {'initial':cycles[0],'final':cycles[-1],'cycles':cycles,
            'reflections_used_total':all_count,'reflections_used_free':free_count,
            'reflections_used_work':all_count-free_count,
            'resolution_low':float(root.findtext('./Overall_stats/resolution_low')),
            'resolution_high':float(root.findtext('./Overall_stats/resolution_high')),
            'r_factor_precision_decimals':4}


def verify_free_flags(before_mtz,after_mtz,free,fp,sigfp):
    def rows(mtz):
        labels=mtz.column_labels()
        result={tuple(r[:3].astype(int)):r for r in np.array(mtz,copy=False)}
        if len(result)!=mtz.nreflections:
            raise ValueError('Duplicate HKLs are not supported by this audit')
        return labels,result
    labels,before=rows(before_mtz);after_labels,after=rows(after_mtz)
    if after.keys()-before.keys():
        raise RuntimeError('REFMAC added unexpected HKLs')
    missing=before.keys()-after.keys()
    for key in missing:
        row=before[key];f=row[labels.index(fp)];sigma=row[labels.index(sigfp)]
        if np.isfinite(f) and np.isfinite(sigma) and f>=0 and sigma>0:
            raise RuntimeError('REFMAC removed an observed reflection')
    for key in after:
        a=before[key][labels.index(free)];b=after[key][after_labels.index(free)]
        if not (a==b or np.isnan(a) and np.isnan(b)):
            raise RuntimeError('REFMAC changed an original free-set flag')
    return {'all_retained_flags_unchanged':True,'removed_unobserved_hkls':[[int(v) for v in k] for k in sorted(missing)],
            'scope':'Unobserved rows may be omitted; observed reflections may not be omitted by this adapter.'}


def refmac_refine(mtz_path,model_path,output,cycles=10,ccp4=None,
                  fp='F',sigfp='SIGF',free='FreeR_flag',free_value=0,timeout=300):
    if not isinstance(cycles,int) or cycles<0:
        raise ValueError('Refinement cycles must be a nonnegative integer')
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+',s) for s in [fp,sigfp,free]):
        raise ValueError('Unsupported MTZ label')
    ccp4=resolve_runtime(ccp4,required=('refmac5',))
    binary=ccp4/'bin/refmac5'
    if not binary.is_file():
        raise ValueError('REFMAC executable unavailable: '+str(binary))
    mtz,partition=reflection_partition(mtz_path,fp,sigfp,free,free_value)
    structure=gemmi.read_structure(str(model_path))
    if len(structure)!=1 or not structure[0]:
        raise ValueError('Exactly one populated model required')
    sg=structure.find_spacegroup()
    if sg is None or sg.number!=mtz.spacegroup.number or not np.allclose(
            structure.cell.parameters,mtz.cell.parameters,rtol=0.,atol=.01):
        raise ValueError('Model/data crystal cell or space group differs')
    out=Path(output).resolve();out.mkdir(parents=True,exist_ok=False)
    structure.write_pdb(str(out/'input.pdb'))
    # Expand bounds slightly to avoid dropping the highest-resolution observations
    # through rounding between MTZ and REFMAC's internal resolution limits.
    resolution_low=partition['observed_resolution_low']*1.0001
    resolution_high=partition['observed_resolution_high']*.9999
    keywords=('LABIN FP=%s SIGFP=%s FREE=%s\nFREE %d\n'%(fp,sigfp,free,free_value)+
              'REFI TYPE REST RESO %.6f %.6f\n'%(resolution_low,resolution_high)+
              'REFI RESI MLKF\nREFI BREF ISOT\n'
              'MAKE HYDR N HOUT N BUIL N\nSCAL TYPE BULK\nSCAL MLSC WORK\n'
              'WEIGHT AUTO\nNCYC %d\nEND\n'%cycles)
    (out/'keywords.txt').write_text(keywords)
    command=[str(binary),'HKLIN',str(Path(mtz_path).resolve()),'XYZIN',str(out/'input.pdb'),
             'HKLOUT',str(out/'refined.mtz'),'XYZOUT',str(out/'refined.pdb'),
             'XMLOUT',str(out/'refmac.xml')]
    report={'version':__version__,'engine':'external-REFMAC5','model_input':fingerprint(model_path),
            'reflections':partition,'binary':fingerprint(binary),'source':fingerprint(__file__),
            'command':command,'keywords':keywords,'requested_cycles':cycles,
            'requested_resolution_limits':[resolution_low,resolution_high],
            'selection':'Fixed final cycle; Rfree never used for selection or early stopping.',
            'scope':'Native amplitudes, original free flags, no target reference or phase restraints. '
                    'ML model-error parameters estimated from work reflections (SCAL MLSC WORK); '
                    'least-squares scaling uses the default working set. '
                    'Prior REFMAC-map building may already contain free-set information; '
                    'these values are not a certified independent blind-test Rfree.'}
    started=time.perf_counter()
    environment=runtime_environment(ccp4,scratch=out/'scratch')
    report['scratch_directory']=environment['CCP4_SCR']
    try:
        with (out/'refmac.log').open('w') as log:
            result=subprocess.run(command,input=keywords,text=True,stdout=log,stderr=subprocess.STDOUT,
                                  cwd=str(out),env=environment,timeout=timeout)
        report['returncode']=result.returncode
        if result.returncode!=0 or not (out/'refined.pdb').is_file() or not (out/'refined.mtz').is_file():
            raise RuntimeError('REFMAC failed; see '+str(out/'refmac.log'))
        stats=parse_refmac((out/'refmac.log').read_text(),out/'refmac.xml')
        stats['usable_observations_not_in_refmac_statistics']=partition['usable_observations']-stats['reflections_used_total']
        if stats['final']['cycle']!=cycles:
            raise RuntimeError('REFMAC did not complete the requested cycle count')
        output_mtz=gemmi.read_mtz_file(str(out/'refined.mtz'))
        flag_audit=verify_free_flags(mtz,output_mtz,free,fp,sigfp)
        refined=gemmi.read_structure(str(out/'refined.pdb'))
        refined.make_mmcif_document().write_file(str(out/'model.cif'))
        components={r.name:monomer(r.name,Path(ccp4)/'lib/data/monomers')
                    for c in refined[0] for r in c if r.name!='UNK'}
        atom_geometry=atom_geometry_report(refined,{'entries':{name:{'component':c} for name,c in components.items()}})
        report['monomer_restraints']=[fingerprint(c['path']) for c in components.values()]
        report.update(status='success',statistics=stats,free_flags_preserved=True,free_flag_audit=flag_audit,
                      atom_geometry=atom_geometry,
                      output_model=fingerprint(out/'model.cif'),output_mtz=fingerprint(out/'refined.mtz'))
    except Exception as error:
        report.update(status='failed',error=str(error))
        raise
    finally:
        report['seconds']=time.perf_counter()-started
        (out/'refinement_report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report
