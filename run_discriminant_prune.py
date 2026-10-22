"""Opt-in synthetic-decoy-calibrated backbone pruning; no reference input."""
import argparse
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
from phasebuild.data import Density, fingerprint
from phasebuild.model import read_backbone
from phasebuild.trace import write_backbone
from phasebuild.backbone_discriminant import BackboneDiscriminant, supported_segments


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('mtz')
    p.add_argument('--model', required=True)
    p.add_argument('--discriminant', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--calibration-decoy-rate', choices=['0.01', '0.05', '0.1', '0.2'], default='0.05')
    p.add_argument('--minimum-length', type=int, default=4)
    p.add_argument('--bridge', type=int, default=1)
    a = p.parse_args()
    out = Path(a.output)
    out.mkdir(parents=True, exist_ok=False)
    started = time.time()
    density = Density(a.mtz)
    model = BackboneDiscriminant.load(a.discriminant)
    chains = read_backbone(a.model, density)
    threshold = model.metadata['calibration']['thresholds'][a.calibration_decoy_rate]['logit']
    scores, retained, intervals = [], [], []
    for chain in chains:
        logits = model.score(density, chain)
        pieces = supported_segments(logits, threshold, a.minimum_length, a.bridge)
        scores.append(logits.tolist())
        intervals.append(pieces)
        retained.extend(chain[start:stop] for start, stop in pieces)
    if not retained:
        raise ValueError('No supported backbone segments remain; inputs were preserved')
    write_backbone(out / 'backbone.cif', retained, density)
    report = {'engine': 'phasebuild-experimental-discriminant-prune', 'parameters': vars(a),
        'input': density.source, 'model_input': fingerprint(a.model),
        'discriminant_input': fingerprint(a.discriminant),
        'threshold': threshold, 'scores': scores, 'retained_intervals': intervals,
        'input_residues': sum(map(len, chains)), 'output_residues': sum(map(len, retained)),
        'output_fragments': len(retained), 'seconds': time.time() - started,
        'warning': 'Backbone-only output requires fresh sequence assignment; synthetic decoy rate is not an experimental error rate.'}
    (out / 'build_report.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({key: report[key] for key in ['input_residues', 'output_residues', 'output_fragments', 'seconds']}), flush=True)


if __name__ == '__main__':
    main()
