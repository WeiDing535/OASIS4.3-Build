"""Assemble local evaluation bundle with explicitly attributed SAD/DM engines.
Copies existing authorized components; does not install or download software.
"""
from pathlib import Path
import argparse,importlib.util,json,shutil,subprocess,os
from build_installation_free import main as base_build,sha,manifest

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['python-source','runtime-source','frontend-source','density-source','output']:p.add_argument('--'+k,type=Path,required=True)
 p.add_argument('--resume-assembly',action='store_true',help='Resume this builder’s newly created incomplete local assembly')
 a=p.parse_args();out=a.output.resolve()
 if not a.resume_assembly:base_build(['--python-source',str(a.python_source),'--runtime-source',str(a.runtime_source),'--output',str(out)])
 runtime=out/'runtime';rows=json.loads((a.frontend_source/'MANIFEST.json').read_text())['files'];lookup={r['path']:r for r in rows};needed=set()
 def visit(name):
  if name in needed:return
  needed.add(name)
  for dep in lookup[name].get('dependencies',[]):
   if dep['target']!='macOS system':visit(dep['target'])
 visit('bin/phaser')
 for name in sorted(needed):
  src=a.frontend_source/name;dst=runtime/name
  if dst.exists():
   if sha(src)!=sha(dst):raise ValueError('Conflicting native library: '+name)
  else:dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
 if not (runtime/'density').exists():shutil.copytree(a.density_source,runtime/'density',symlinks=True)
 (runtime/'density/COMPONENT_SCOPE.json').write_text(json.dumps({'origin':'Existing locally authorized distribution; Phenix-derived RESOLVE/CCTBX and Python runtime components','full_Phenix_suite_bundled':False,'requires_separate_installation':False,'public_redistribution_ready':False},indent=2)+'\n')
 m=json.loads((runtime/'runtime_manifest.json').read_text());m['files']=[r for r in manifest(runtime) if r['relative_path']!='runtime_manifest.json'];m.update(bundled_external_engines=['Phaser','RESOLVE','CCTBX','Buccaneer','REFMAC'],public_redistribution_ready=False)
 (runtime/'runtime_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
 doc=out/'app/portable_doctor.py';s=doc.read_text();s=s.replace("cfg['method'] in ('pplus_internal','pplus_bundled')","cfg['method']=='pplus_bundled'")
 s=s.replace("['Buccaneer','REFMAC']","['Phaser','RESOLVE','CCTBX','Buccaneer','REFMAC']")
 s=s.replace('print(json.dumps(',"from phasebuild.bundled_pplus import engine_layout\nengine_layout(base/'runtime')\nprint(json.dumps(")
 doc.write_text(s)
 density_probe='import sys,libtbx.load_env,cctbx,mmtbx.command_line.density_modification;from solve_resolve import resolve;print(sys.executable);print(libtbx.env.build_path);print(cctbx.__file__);print(resolve.__file__);from cctbx import crystal,miller;from cctbx.array_family import flex;ma=miller.array(miller.set(crystal.symmetry(unit_cell=(20,20,30,90,90,90),space_group_symbol="P 41"),flex.miller_index([(1,0,0)]),False),flex.double([1.]));ma.as_mtz_dataset(column_root_label="F")'
 run=subprocess.run([str(runtime/'density/bin/oasis-density-python'),'-c',density_probe],env={**os.environ,'CLIBD':str(runtime/'lib/data')},capture_output=True,text=True,timeout=120)
 (out/'DENSITY_DOCTOR.txt').write_text(run.stdout+run.stderr)
 if run.returncode:raise RuntimeError('Density runtime doctor failed')
 run=subprocess.run([str(out/'oasis'),'doctor'],capture_output=True,text=True,timeout=120)
 (out/'DOCTOR.json').write_text(run.stdout)
 if run.returncode:raise RuntimeError('Bundle doctor failed: '+run.stderr)
 report=json.loads((out/'bundle.json').read_text());report.update(default_method='pplus_bundled',default_execution='pplus_w8_weight',GUI_portfolio=False,automatic_fpp=True,explicit_wavelength_required=True,Phenix_bundled=False,Phenix_derived_components_bundled=True,full_Phenix_suite_bundled=False,bundled_external_engines=['Phaser','RESOLVE','CCTBX','Buccaneer','REFMAC'],scientific_scope='Fixed bundled P+ Wang8/classic40 with sequence cleanup and work-only geometry weight search; theoretical fpp from elements and explicit mandatory wavelength',engine_manifest_SHA256=sha(runtime/'runtime_manifest.json'),doctor=json.loads(run.stdout))
 (out/'bundle.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'使用说明.txt').write_text('OASIS4.3 '+report['version']+' 本机授权评估包（macOS Intel）\n双击 OASIS4.3.command，保持目录完整。无需另装 Python、Phenix 或 CCP4。\n包含第三方 Phaser、RESOLVE、CCTBX、Buccaneer、REFMAC 及其运行库，并非全部算法由 OASIS 原创。RESOLVE/CCTBX 组件来源于已有授权发行版，并非完整 Phenix 套件。\n新建 GUI 任务仅运行固定 P+ Wang8＋工作集精修权重搜索；填写元素和必填实验波长，程序自动计算理论 f″。\n请填写实验波长，输入原始工作/自由集、序列和重原子子结构。不会自动猜测波长或输入目标蛋白模型。\n仅供本机授权研究评估，公开再分发条件尚未落实。单个开发案例成功不保证其他样品均达到相同 Rfree。\n')
 app=out/'app';(app/'app_manifest.json').write_text(json.dumps({'files':[r for r in manifest(app) if r['relative_path']!='app_manifest.json']},indent=2)+'\n')
 print(json.dumps({'status':'assembled','bundle':str(out)},indent=2))
if __name__=='__main__':main()
