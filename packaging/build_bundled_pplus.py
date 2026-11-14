"""Assemble local evaluation bundle with explicitly attributed SAD/DM engines.
Copies existing authorized components; does not install or download software.
"""
from pathlib import Path
import argparse,importlib.util,json,shutil,subprocess,os
from build_installation_free import main as base_build,sha,manifest

def main():
 p=argparse.ArgumentParser(description=__doc__)
 for k in ['python-source','runtime-source','frontend-source','output']:p.add_argument('--'+k,type=Path,required=True)
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
 # No density/ tree: P+ density-modifies with the in-tree handdm v2 engine, so
 # the bundle holds no Phenix-derived or external density-modification component.
 programs={name:(runtime/'bin'/name).is_file() for name in ('phaser','cbuccaneer','refmac5')}
 if not all(programs.values()):raise RuntimeError('Runtime is missing P+ programs: '+', '.join(sorted(k for k,v in programs.items() if not v)))
 labels=[label for name,label in (('phaser','Phaser'),('cbuccaneer','Buccaneer'),('refmac5','REFMAC')) if (runtime/'bin'/name).is_file()]
 m=json.loads((runtime/'runtime_manifest.json').read_text());m['files']=[r for r in manifest(runtime) if r['relative_path']!='runtime_manifest.json'];m.update(bundled_external_engines=labels,public_redistribution_ready=False)
 (runtime/'runtime_manifest.json').write_text(json.dumps(m,indent=2)+'\n')
 doc=out/'app/portable_doctor.py';s=doc.read_text();s=s.replace("cfg['method'] in ('pplus_internal','pplus_bundled')","cfg['method']=='pplus_bundled'")
 s=s.replace("['Buccaneer','REFMAC']","['Phaser','Buccaneer','REFMAC']")
 s=s.replace('print(json.dumps(',"from phasebuild.bundled_pplus import engine_layout\nengine_layout(base/'runtime')\nprint(json.dumps(")
 doc.write_text(s)
 run=subprocess.run([str(out/'oasis'),'doctor'],capture_output=True,text=True,timeout=120)
 (out/'DOCTOR.json').write_text(run.stdout)
 if run.returncode:raise RuntimeError('Bundle doctor failed: '+run.stderr)
 report=json.loads((out/'bundle.json').read_text());report.update(default_method='pplus_bundled',default_execution='pplus_w8_weight',GUI_portfolio=False,automatic_fpp=True,explicit_wavelength_required=True,Phenix_bundled=False,Phenix_derived_components_bundled=False,full_Phenix_suite_bundled=False,bundled_external_engines=['Phaser','Buccaneer','REFMAC'],ccp4_programs=programs,scientific_scope='Bundled P+ with in-tree handdm v2 density modification, sequence cleanup and work-only geometry weight search; theoretical fpp from elements and explicit mandatory wavelength',engine_manifest_SHA256=sha(runtime/'runtime_manifest.json'),doctor=json.loads(run.stdout))
 (out/'bundle.json').write_text(json.dumps(report,indent=2)+'\n')
 (out/'使用说明.txt').write_text('OASIS4.3 '+report['version']+' 本机授权评估包（macOS Intel）\n双击 OASIS4.3.command，保持目录完整。无需另装 Python、Phenix 或 CCP4。\n包含第三方 Phaser、Buccaneer、REFMAC 及其运行库，并非全部算法由 OASIS 原创。密度修饰由随程序的 handdm v2 完成，包内不含 Phenix、Phenix 派生组件或外部密度修饰程序。\n新建 GUI 任务仅运行固定 P+ Wang8＋工作集精修权重搜索；填写元素和必填实验波长，程序自动计算理论 f″。\n请填写实验波长，输入原始工作/自由集、序列和重原子子结构。不会自动猜测波长或输入目标蛋白模型。\n仅供本机授权研究评估，公开再分发条件尚未落实。单个开发案例成功不保证其他样品均达到相同 Rfree。\n')
 app=out/'app';(app/'app_manifest.json').write_text(json.dumps({'files':[r for r in manifest(app) if r['relative_path']!='app_manifest.json']},indent=2)+'\n')
 print(json.dumps({'status':'assembled','bundle':str(out)},indent=2))
if __name__=='__main__':main()
