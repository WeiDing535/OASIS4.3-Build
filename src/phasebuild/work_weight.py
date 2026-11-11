"""Opt-in restraint-weight calibration without original test observations.

Pilots contain no measured test amplitudes. Selection uses fixed geometry
limits only. Final refinement must restart from the pre-pilot model and use
untouched original native data. This does not make a previously built model
an end-to-end blind test.
"""
from pathlib import Path
import json,subprocess,shutil,xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
import gemmi
import numpy as np
from .runtime import runtime_environment
from .data import fingerprint


def mask_test_observations(native_file,output,free_value=0):
    from .sad_bridge import columns
    m=gemmi.read_mtz_file(str(native_file));labels=columns(m,{'F':'F','SIGF':'Q','FreeR_flag':'I'});a=np.array(m,copy=True);original=a.copy();free=a[:,labels.index('FreeR_flag')]==free_value
    if not free.any():raise ValueError('An explicit original test set is required')
    for j,label in enumerate(labels):
        if label not in ('H','K','L','FreeR_flag'):a[free,j]=np.nan
    m.set_data(a);m.write_to_file(str(output))
    if not np.array_equal(a[~free],original[~free],equal_nan=True) or not np.array_equal(a[:,labels.index('FreeR_flag')],original[:,labels.index('FreeR_flag')],equal_nan=True):raise RuntimeError('Pilot working data or original flags changed')
    return {'original_test_rows_masked':int(free.sum()),'working_rows_unchanged':True,'original_flags_unchanged':True,'all_nonindex_nonflag_test_columns_masked':True,'source':fingerprint(native_file),'output':fingerprint(output)}


def select_geometry_weight(reports,bond_limit=.010,angle_limit=1.3,default=.005):
    eligible=[]
    for r in reports:
        values=[r.get('weight'),r.get('bond_rmsd_angstrom'),r.get('angle_rmsd_degrees')]
        if r.get('status')=='completed' and r.get('cycle')==30 and all(isinstance(v,(int,float)) and np.isfinite(v) for v in values) and values[0]>0 and 0<=values[1]<=bond_limit and 0<=values[2]<=angle_limit:eligible.append(values[0])
    return max(eligible) if eligible else default


def parse_work_geometry(log_file,xml_file):
    root=ET.parse(xml_file).getroot()
    if int(root.findtext('err_level','999'))!=0:raise RuntimeError('Work-only pilot failed')
    if int(root.findtext('./Overall_stats/n_reflections_free','999'))>0:raise RuntimeError('Work pilot unexpectedly contains observed test reflections')
    cycles=root.findall('./Overall_stats/stats_vs_cycle/new_cycle')
    if [int(c.findtext('cycle')) for c in cycles]!=list(range(31)):raise RuntimeError('Work pilot did not complete fixed30 cycles')
    last=cycles[-1];result={'cycle':30,'r_work':float(last.findtext('r_factor')),'bond_rmsd_angstrom':float(last.findtext('rmsBOND')),'angle_rmsd_degrees':float(last.findtext('rmsANGLE'))}
    if not np.isfinite(list(result.values())).all():raise RuntimeError('Nonfinite work pilot geometry')
    rows=[]
    sections=Path(log_file).read_text().split('Ncyc    Rfact    Rfree')
    if len(sections)<2:raise RuntimeError('Missing work pilot cycle table')
    for line in sections[-1].split('$TEXT:')[0].splitlines():
        fields=line.split()
        if len(fields)>=11 and fields[0].isdigit():
            try:rows.append((int(fields[0]),float(fields[1]),float(fields[6]),float(fields[8])))
            except ValueError:continue
    if not rows or rows[-1][0]!=30 or abs(rows[-1][1]-result['r_work'])>.0011:raise RuntimeError('Work pilot XML/log statistics disagree')
    # Deliberately never parse or return the absent-test Rfree sentinel.
    return result


def calibrate_work_weight(native_file,model_file,runtime,output,free_value=0):
    from .refinement import reflection_partition
    out=Path(output);out.mkdir(parents=True,exist_ok=False);masked=out/'work_only.mtz';mask=mask_test_observations(native_file,masked,free_value)
    _,part=reflection_partition(native_file,free_value=free_value)
    def pilot(weight):
        folder=out/('w'+str(weight).replace('.','p'));folder.mkdir();shutil.copy2(masked,folder/'input.mtz');gemmi.read_structure(str(model_file)).write_pdb(str(folder/'input.pdb'))
        text=('LABIN FP=F SIGFP=SIGF FREE=FreeR_flag\nFREE %d\nREFI TYPE REST RESO %.6f %.6f\nREFI RESI MLKF\nREFI BREF ISOT\nMAKE HYDR N HOUT N BUIL N\nSCAL TYPE BULK\nSCAL MLSC WORK\nWEIGHT MATRIX %.6g\nNCYC 30\nEND\n')%(free_value,part['observed_resolution_low']*1.0001,part['observed_resolution_high']*.9999,weight)
        (folder/'keywords.txt').write_text(text);cmd=[str(Path(runtime)/'bin/refmac5'),'HKLIN','input.mtz','XYZIN','input.pdb','HKLOUT','refined.mtz','XYZOUT','refined.pdb','XMLOUT','refmac.xml'];(folder/'command.json').write_text(json.dumps(cmd,indent=2)+'\n');env=runtime_environment(runtime,scratch=folder/'scratch');env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
        with (folder/'refmac.log').open('x') as f:r=subprocess.run(cmd,cwd=folder,env=env,input=text,text=True,stdout=f,stderr=subprocess.STDOUT,timeout=1800)
        if r.returncode:raise RuntimeError('Work-only weight pilot failed: '+str(weight))
        result={'status':'completed','weight':weight,**parse_work_geometry(folder/'refmac.log',folder/'refmac.xml')};(folder/'GEOMETRY_REPORT.json').write_text(json.dumps(result,indent=2)+'\n');return result
    with ThreadPoolExecutor(max_workers=2) as pool:reports=list(pool.map(pilot,[.005,.01,.02,.04]))
    selected=select_geometry_weight(reports);report={'selected_weight':selected,'policy':'Largest of0.005,0.01,0.02,0.04 meeting bondRMSD<=0.010A and angleRMSD<=1.3deg at fixed30 cycles; fallback0.005','pilots':reports,'mask_audit':mask,'prepilot_model':fingerprint(model_file),'pilot_models_reused':False,'selection_uses_original_Rfree':False,'scope':'No original test observations supplied to pilots. Previously generated model/maps are not certified end-to-end work-only.'};(out/'SELECTION.json').write_text(json.dumps(report,indent=2)+'\n');return selected,report
