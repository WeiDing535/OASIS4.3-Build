"""Post-build reference evaluation and fixed-stage comparison; never a builder input."""
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent / 'src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from phasebuild.data import fingerprint
from phasebuild.evaluate import evaluate
from phasebuild.sidechain_evaluate import evaluate_sidechains

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'runs/comparison_v0p8'
REFERENCE = Path('/Volumes/For_data/IPCAS-test/rebuilt/INS/ins_reference.pdb')


def read(path):
    return json.loads(Path(path).read_text())


def row(label, cycles, report_path):
    report_path = Path(report_path)
    report = read(report_path)
    model = report_path.parent / 'model.cif'
    backbone = evaluate(model, REFERENCE)
    sidechains = evaluate_sidechains(model, REFERENCE)
    # Apply the same current chemistry audit to every comparator, including Pro.
    import gemmi
    from phasebuild.atom_validation import atom_geometry_report
    from phasebuild.rotamers import load_library
    geometry = atom_geometry_report(gemmi.read_structure(str(model)),
                                   load_library(ROOT / 'models/rotamers_1tqw_v08.json'))
    assert report['status'] == 'success' and report['free_flags_preserved']
    assert report['reflections']['free_value'] == 0
    assert report['statistics']['reflections_used_free'] == 299
    return dict(label=label, cycles=cycles, model=fingerprint(model),
                refinement_report=fingerprint(report_path),
                r_work=report['statistics']['final']['r_work'],
                r_free=report['statistics']['final']['r_free'],
                backbone=backbone, sidechains=sidechains, geometry=geometry)


def main():
    OUT.mkdir(exist_ok=True)
    baseline = ROOT / 'runs/ins_v025_phasebuild_v07_refined/refinement_report.json'
    rows = [row('PhaseBuild 0.7 baseline', 20, baseline)]
    curves = {}
    build_audits = []
    for folder, label in [('ins_v038_v081_primary', '0.8.1 original-map restart'),
                          ('ins_v039_v081_updated_only', '0.8.1 updated-map only')]:
        stages = [rows[0]]
        for cycles, part in [(30, 'cycles_30'), (40, 'continuation/cycles_40'),
                             (50, 'continuation/cycles_50')]:
            root = ROOT / 'runs' / folder / part
            item = row(label, cycles, root / 'refine/refinement_report.json')
            rows.append(item)
            stages.append(item)
            report = read(root / 'build/build_report.json')
            assert not report['source_changed_during_run']
            assert all(fingerprint(source['path']) == source for source in report['sources'])
            build_audits.append(dict(report=fingerprint(root / 'build/build_report.json'),
                                     input_map=report['input'], extension=report['extension'],
                                     seconds=report['seconds']))
        assert not read(ROOT / 'runs' / folder / 'continuation/experiment_report.json')['source_changed_during_run']
        curves[label] = stages
    old = read(ROOT / 'runs/ins_v031_updated_map_50/experiment_report.json')
    curves['0.7 updated-map only'] = old['stages']
    rows.append(row('PhaseBuild 0.7 updated-map only', 50,
                    old['final_refinement_report']['path']))
    bc = read(ROOT / 'runs/ins_v028_buccaneer_pipeline/pipeline_report.json')
    curves['Buccaneer correlation pipeline'] = [dict(s, cycles=10*s['round']) for s in bc['stages']]
    rows.append(row('Buccaneer correlation pipeline', 50, bc['stages'][-1]['refinement_report']['path']))
    protocols = [read(r['refinement_report']['path']) for r in rows]
    assert len({r['reflections']['input']['sha256'] for r in protocols}) == 1
    assert len({tuple(r['keywords']) for r in protocols}) == 1
    result = dict(rows=rows, build_audits=build_audits,
                  map_diagnosis=read(ROOT / 'runs/diagnosis_v08/map_diagnosis.json'),
                  reference=fingerprint(REFERENCE),
                  scope='Development set, model-derived input phases. Fixed final stage; no Rfree selection. '
                        'Main arm restarts from original supplied map at cumulative cycle 20; '
                        'Buccaneer subsequently uses its own updated maps. Not identical map schedules. '
                        'Reference access is post-build evaluation only.')
    (OUT / 'comparison.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    for label, stages in curves.items():
        x = [s.get('cycles', s.get('cumulative_cycles')) for s in stages]
        ax.plot(x, [s['r_free']*100 for s in stages], 'o-', label=label)
    ax.set(xlabel='Cumulative REFMAC cycles', ylabel='Rfree (%)', xticks=[10,20,30,40,50],
           title='INS: fixed-stage comparison (model-derived phases)')
    ax.grid(alpha=.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    for suffix in ['png','svg']:
        fig.savefig(str(OUT / ('rfree_cycles.'+suffix)), dpi=170)
    plt.close(fig)
    table = ['| Workflow | Cumulative refinement rounds | Rwork | Rfree | Built residues | N/CA/C all <1 Å | Correct side-chain atoms/built |',
             '|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        b,s = r['backbone'],r['sidechains']
        table.append('| {} | {} | {:.2f}% | {:.2f}% | {} | {} | {}/{} |'.format(
            r['label'],r['cycles'],r['r_work']*100,r['r_free']*100,b['built_residues'],
            b['n_ca_c_correct_residues'],s['sidechain_atoms_within_cutoff'],s['sidechain_atoms_built']))
    report = '''# PhaseBuild 0.8.1: terminal extension versus phase-map comparison

On the INS development data, the new workflow lowers the fixed 50-round Rfree to **27.05%**, versus **26.91%** for the Buccaneer comparison, a difference of 0.14 percentage points. The new model has 49/49 residues with strict backbone correctness, covering 96.1% of the 51 reference residues; this does not justify claiming parity with or superiority to Buccaneer on independent data.

## Results

TABLE

Strict correctness means N, CA and C are each within <1 Å of the reference under the same crystal symmetry image; side chains use the same threshold and allow chemically equivalent atom names to be interchanged. All residues are matched one-to-one by CA, and the post-match RMSD must not be treated as whole-model accuracy. In the final main run, all 49 residue types match the reference, covering A1–A21 and B2–B29, with only B1 and B30 missing. The two chains correspond respectively to the same reference chain and the same symmetry image, with matched residues in increasing order. The current comparison uses the full correlation-mode Buccaneer construction workflow; other historical modes are in comparison_v0p7, and the weaker default mode was not used alone for comparison.

![Fixed-stage Rfree curve](rfree_cycles.png)

## Where the gap was and what was fixed

1. The strictly correct portion of the old model was concentrated in A1–A21 and B4–B19, with both ends of chain B missing. The old local cycles could not generate new terminal backbone; increasing the refinement rounds still stopped at 38 residues. We have now added sequence-anchored bidirectional extension, multipath torsion-angle search, repair of one terminal residue, continuous optimization and density/geometry/symmetry-clash screening. Internal anchors are fixed during extension, and REFMAC then optimizes the whole model. The official Buccaneer workflow also includes multiple steps such as find/grow/join/link/sequence/correct; iterative refinement is not the same as iterative completion; see [official documentation](https://www.ccp4.ac.uk/html/cbuccaneer.html). This program has not ported or decompiled its code.
2. An incomplete model's own updated map weakens the information for regions not yet built. In this case, the reference B20–B30 positions were used only for post-hoc diagnostic sampling; the mean backbone atom densities of the originally supplied map and the old model's updated map are 3.10 and 1.97, with profiles of 3.49 and 2.89. The main workflow reused the originally supplied map once at the round-20 starting point, then used only its own updated map twice, finally reaching 49 residues and 27.05%. The all-self-updated-map branch gave 42 residues and 35.65%. This is strong evidence of workflow dependence; the extension algorithm and the input phase information jointly contributed to the improvement.
3. Synthetic validation of positive Gly φ exposed insufficient sampling of the terminal O direction, which has now been added. The Pro monomer dictionary uses a different bond-type column name, and the old parser missed ring linkages; it now supports both formats and reports an error when the bond graph is missing. A regenerated homologous 1TQW conformer library was used, not the INS reference conformers. Fixing the dictionary does not mean all Pro side chains were successfully completed.

## Protocol, quality and reproducibility

All refinements in the table use the same native F/SIGF, original FreeR_flag=0, 299 free-set reflections and the same REFMAC keywords; the final stage is fixed and no round is chosen by minimum Rfree. Of 6148 available observations, REFMAC statistics used 6146 (5847 work + 299 free); the reason the other two were excluded has not been fully established, but the statistics agree across groups. The free-set flags of the retained HKLs were verified unchanged.

The final refinement of the main workflow has a bond-length RMSD of 0.0112 Å and a bond-angle RMSD of 1.794°; the current simplified geometry audit found no severe nonbonded heavy-atom contacts below 1.8 Å. This is not equivalent to full MolProbity/Ramachandran validation. Of the 49 residues, 19 built side-chain atoms still fail to reach the reference 1 Å threshold; termini and side chains still need improvement. All comparison objects were re-audited with the same repaired chemical library; see comparison.json for details.

51 tests pass; see ../verification_v0p8_tests.xml. Each construction saves the actual input map hash, source hashes, sequence anchors and new events; the summary checks that the current source matches the run record. The pipeline_report.json of run directories ins_v038_v081_primary and ins_v039_v081_updated_only are both success. The 0.8.0 exploratory runs ins_v032–ins_v037 are retained and labelled and must not be merged with the final 0.8.1 into a single trajectory. From the project root, run `ccp4-python summarize_extension_experiment.py` to rebuild this report and the figures.

The input phases come from an existing REFMAC model and may already contain free-set/model information; these data are also used for development. The builder does not read the target reference or Buccaneer coordinates; the reference is accessed only by the post-hoc evaluator. The map-use schedules of the main workflow and Buccaneer differ; 50 cumulative rounds mean only the same number of refinements, not the same total computational budget. A 0.14-percentage-point Rfree difference is insufficient for a statistical judgement of superiority.

## Next steps

The priority is to test the "retain original phase constraints" iterative strategy on a multi-protein test set with independent experimental phases, then proceed to gap connection, whole-map residue discovery and side-chain optimization. Phase constraints should come from OASIS experimental phase/quality information with explicit free-set separation, rather than repeatedly relying on the refined map of some reference model. What this round has verified is the completion capability and the importance of preserving map information; generalization from SAD to a complete model has not yet been validated.
'''.replace('TABLE','\n'.join(table))
    (OUT / 'REPORT.md').write_text(report)
    print('\n'.join(table))


if __name__ == '__main__':
    main()
