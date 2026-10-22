"""Select from autonomous candidates plus an explicitly supplied external model.

The output is labelled hybrid; it is not an autonomous-tracer benchmark.
"""
import argparse
import json
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
import gemmi
import numpy as np
from phasebuild.data import Density, fingerprint
from phasebuild.profile import DensityProfile
from phasebuild.selection import select_fragments
from phasebuild.trace import write_backbone
from phasebuild.geometry import oxygen


def model_fragments(path, density):
    structure = gemmi.read_structure(str(path))
    if not np.allclose(structure.cell.parameters, density.cell.parameters, rtol=1e-4, atol=1e-3):
        raise ValueError('External model has a different unit cell')
    sg = gemmi.find_spacegroup_by_name(structure.spacegroup_hm)
    if sg is None or sg.number != density.sg.number:
        raise ValueError('External model has a different space group')
    pieces = []
    for chain in structure[0]:
        current = []
        for residue in chain:
            atoms = [residue.find_atom(name,'*') for name in ['N','CA','C']]
            xyz = None
            if all(a is not None for a in atoms):
                o = residue.find_atom('O','*')
                xyz = np.array([list(a.pos) for a in atoms] +
                               [list(o.pos) if o is not None else [np.nan]*3])
            broken = xyz is None or (current and not 1.1 < np.linalg.norm(current[-1][2]-xyz[0]) < 1.9)
            if broken:
                if len(current) >= 3:
                    pieces.append(np.array(current))
                current = []
            if xyz is not None:
                current.append(xyz)
        if len(current) >= 3:
            pieces.append(np.array(current))
    # Buccaneer's backbone-only mode commonly writes N/CA/C without O.
    for piece in pieces:
        for i, residue in enumerate(piece):
            if np.isfinite(residue[3]).all():
                continue
            if i + 1 < len(piece):
                points = [gemmi.Position(*p) for p in [*residue[:3], piece[i+1,0]]]
                psi = np.rad2deg(gemmi.calculate_dihedral(*points))
            else:
                # Terminal oxygen is provisional: choose its density-supported angle.
                angles = np.arange(-180.,180.,10.)
                probes = np.array([oxygen(*residue[:3], angle) for angle in angles])
                psi = angles[np.argmax(density.sample(probes))]
            residue[3] = oxygen(*residue[:3], psi)
    return pieces


def main():
    p = argparse.ArgumentParser()
    p.add_argument('mtz')
    p.add_argument('--profile', required=True)
    p.add_argument('--candidates', required=True)
    p.add_argument('--external-model', required=True)
    p.add_argument('--output', required=True)
    p.add_argument('--cutoff', type=float, default=2.5)
    p.add_argument('--max-residues', type=int, default=65)
    args = p.parse_args()
    started = time.perf_counter()
    density = Density(args.mtz)
    density.profile = DensityProfile.load(args.profile)
    with np.load(args.candidates, allow_pickle=False) as data:
        if json.loads(str(data['metadata']))['input']['sha256'] != density.source['sha256']:
            raise ValueError('Candidate archive belongs to another MTZ')
        candidates = list(data['seeds'])
        start = 0
        for length in data['grown_lengths']:
            candidates.append(data['grown'][start:start+length])
            start += length
    additional = model_fragments(args.external_model, density)
    if not additional:
        raise ValueError('External model contains no usable backbone fragments; refusing silent native-only fallback')
    chains, report = select_fragments(density, candidates + additional, args.max_residues, args.cutoff)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    write_backbone(out / 'backbone.cif', chains, density)
    report.update(engine='hybrid-phasebuild-and-external-model', input=density.source,
                  external_model=fingerprint(args.external_model), external_fragments=len(additional),
                  missing_oxygen_policy='Internal O reconstructed from next N; terminal O provisionally density-fitted. No coordinate refinement.',
                  candidates=fingerprint(args.candidates), profile=fingerprint(args.profile),
                  runner=fingerprint(__file__), parameters=vars(args),
                  scope='Includes an external model. Not an autonomous-tracer result. No INS evaluation reference is read.',
                  total_seconds=time.perf_counter()-started)
    (out / 'build_report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
