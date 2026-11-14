"""Execute one auditable trusted P+ GUI stage from its persisted configuration."""
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild import trust_workflow

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['phase','dm','bridge','finish']);p.add_argument('--config',required=True);p.add_argument('--output',required=True)
    a=p.parse_args();cfg=json.loads(Path(a.config).read_text())['config'];cfg['prepared_sad']=str(Path(a.config).parent/'01_inputs/sad.mtz')
    if cfg.get('input_columns'):
        cfg['native']=str(Path(a.config).parent/'01_inputs/native.mtz');cfg['labels']=''
    if cfg.get('auto_hand') and cfg['method']=='pplus_internal':
        selected=Path(a.config).parent/'02_hand'
        cfg['native']=str(selected/'selected_native.mtz')
        cfg['prepared_sad']=str(selected/'selected_sad.mtz')
        cfg['heavy']=str(selected/'refined_heavy.pdb')
        cfg['labels']=''
    if cfg.get('auto_hand') and cfg['method']=='pplus_bundled' and a.stage!='phase':
        cfg['native']=str(Path(a.config).parent/'02_sad/selected_native.mtz')
        cfg['prepared_sad']=str(Path(a.config).parent/'02_sad/selected_sad.mtz')
    out=Path(a.output)
    if not out.is_dir() or any(out.iterdir()):raise FileExistsError('Trusted P+ stage requires a new, empty stage directory: '+str(out))
    module=trust_workflow
    if cfg['method']=='pplus_internal' and a.stage in ('phase','dm'):
        from phasebuild import internal_pplus
        module=internal_pplus
    if cfg['method']=='pplus_bundled' and a.stage in ('phase','dm'):
        from phasebuild import bundled_pplus
        module=bundled_pplus
    if a.stage=='dm' and module is trust_workflow:
        # P+ density-modifies with the in-tree handdm v2, which sad_pipeline runs
        # through run_trace_feedback.py instead.
        raise ValueError('P+ 两条流程的密度修饰由树内 handdm v2 完成（--dm-engine handdm），'
                         '没有可单独调用的内置 dm 阶段')
    getattr(module,a.stage)(cfg,out)

if __name__=='__main__':main()
