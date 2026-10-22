"""Experiment-declared SeMet conversion; no coordinates are inferred or added."""
import json
from pathlib import Path
import gemmi
from .data import fingerprint


def apply_semet(structure, enabled=False):
    if not isinstance(enabled, bool):
        raise ValueError('SeMet conversion requires an explicit boolean experiment flag')
    result = structure.clone()
    changed, atoms, incomplete = 0, 0, 0
    for model in result:
        for chain in model:
            for residue in chain:
                if not enabled or residue.name != 'MET':
                    continue
                residue.name = 'MSE'
                changed += 1
                count = 0
                for atom in residue:
                    if atom.name == 'SD':
                        atom.name = 'SE'
                        atom.element = gemmi.Element('Se')
                        count += 1
                atoms += count
                incomplete += int(count == 0)
    return result, {'experiment_build_semet': enabled,
                    'converted_met_residues': changed, 'converted_sd_atoms': atoms,
                    'converted_residues_without_sd': incomplete,
                    'scope': 'Only an explicit experiment flag changes MET/SD to MSE/SE; positions, occupancies and B factors are preserved. Missing atoms are never created.'}


def prepare_semet_model(model_path, output, enabled=False):
    out = Path(output)
    if out.exists():
        raise ValueError('Chemistry output already exists')
    original = gemmi.read_structure(str(model_path))
    result, report = apply_semet(original, enabled)
    out.parent.mkdir(parents=True, exist_ok=True)
    result.make_mmcif_document().write_file(str(out))
    report.update(input_model=fingerprint(model_path), output_model=fingerprint(out),
                  source=fingerprint(__file__))
    out.with_suffix('.chemistry.json').write_text(json.dumps(report, indent=2)+'\n')
    return report
