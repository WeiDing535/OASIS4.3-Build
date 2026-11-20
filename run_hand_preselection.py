"""Select the anomalous hand with bundled work-only sequence tracing.

The selected substructure and indexed SAD/native files are then passed to the
user's chosen phasing algorithm. This scout never supplies its protein model.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parent/'src'))
from phasebuild.bundled_pplus import phase

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    state=Path(args.config)
    cfg=json.loads(state.read_text())['config']
    if cfg['method'] not in ('pplus_internal','likelihood','pplus') or not cfg['auto_hand']:
        raise ValueError('Hand preselection requires a supported non-bundled auto-hand workflow')
    prepared=state.parent/'01_inputs'
    cfg.update(method='pplus_bundled',prepared_sad=str(prepared/'sad.mtz'),
               native=str(prepared/'native.mtz'),labels='',auto_hand=True,
               building_recipe='ccp4',trust_wang_radius=8.0)
    out=Path(args.output)
    if not out.is_dir() or any(out.iterdir()):raise FileExistsError(out)
    phase(cfg,out)
    for name in ('hand_selection.json','selected_sad.mtz','selected_native.mtz','refined_heavy.pdb'):
        if not (out/name).is_file():raise RuntimeError('Hand preselection missing '+name)

if __name__=='__main__':main()
