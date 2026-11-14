"""Inventory the ten campaign SAD cases for building a REFMAC-like refiner and a
Buccaneer-like builder.  Writes work/tencase_inventory.json.

Source of truth for the runs is the calibrated ten-case campaign tree
(work/exp/tencase_gui_20261110/<PDB>/00_sad_to_model/endpoint_*/); the older
cases/ tree is only a fallback.

Run:  ./.conda-unified/bin/python work/bundle_tencase.py
"""
import json, os, re
from pathlib import Path
import gemmi
import numpy as np

CAMPAIGN = Path('/Volumes/Fast_4T/Deepseek-Oasis/work/exp/tencase_gui_20261110')
CASES_ROOT = Path('/Volumes/Fast_4T/Deepseek-Oasis/cases')
SRC_ROOT = Path('/Volumes/For_data/IPCAS-test/100_SAD_data')
OUT = Path(__file__).with_name('tencase_inventory.json')


def source_dir(pdbid):
    for d in SRC_ROOT.iterdir():
        if d.is_dir() and pdbid.lower() in d.name.lower():
            return d
    return None


def mtz_meta(path):
    m = gemmi.read_mtz_file(str(path))
    return dict(path=str(path), spacegroup=m.spacegroup.hm,
                cell=[round(x, 3) for x in m.cell.parameters], nref=m.nreflections,
                dmin=round(m.resolution_high(), 3), dmax=round(m.resolution_low(), 3),
                columns=[c.label for c in m.columns])


def structure_stats(path):
    st = gemmi.read_structure(str(path))
    st.remove_alternative_conformations(); st.remove_hydrogens(); st.remove_empty_chains()
    chains = [ch for ch in st[0] if len(ch)]
    residues = [r for ch in chains for r in ch]
    mass = sum(a.element.weight for ch in st[0] for r in ch for a in r)
    z = len(gemmi.SpaceGroup(st.spacegroup_hm).operations())
    prot_vol = 1.2288 * mass * z          # ASU content x Z copies per cell
    return dict(path=str(path), chains=len(chains), residues=len(residues),
                ca=sum(1 for r in residues if r.find_atom('CA', '*') is not None),
                solvent=round(1.0 - prot_vol / st.cell.volume, 3), z=z,
                mass_da=round(mass),
                cell=[round(x, 3) for x in st.cell.parameters], spacegroup=st.spacegroup_hm)


def sfcalc_phases(truth, mtz):
    st = gemmi.read_structure(str(truth))
    st.remove_alternative_conformations(); st.remove_hydrogens()
    st.remove_ligands_and_waters(); st.remove_empty_chains()
    m = gemmi.read_mtz_file(str(mtz))
    calc = gemmi.StructureFactorCalculatorX(st.cell)
    a = np.array(m, copy=False)
    out = np.zeros(len(a))
    for i in range(len(a)):
        out[i] = np.angle(calc.calculate_sf_from_model(st[0], [int(a[i, 0]), int(a[i, 1]), int(a[i, 2])]))
    return out


def phase_cc(mtz, truth):
    m = gemmi.read_mtz_file(str(mtz))
    labels = {c.label.upper(): i for i, c in enumerate(m.columns)}
    phi_idx = next((labels[l] for l in ('PHIB', 'PHIC', 'PHI') if l in labels), None)
    fom_idx = next((labels[l] for l in ('FOM', 'FOMM') if l in labels), None)
    if phi_idx is None:
        return None
    a = np.array(m, copy=False)
    phi = np.deg2rad(a[:, phi_idx])          # MTZ stores phases in degrees
    w = a[:, fom_idx] if fom_idx is not None else np.ones(len(m))
    fcph = sfcalc_phases(truth, mtz)
    d = np.angle(np.exp(1j * (phi - fcph)))
    c = float(np.sum(w * np.cos(d)) / np.sum(w))
    return dict(n=int(len(a)), cc=round(c, 4), cc_abs=round(abs(c), 4),
                within45=round(float(np.sum(w * (np.abs(d) < np.pi / 4)) / np.sum(w)), 4),
                hand='same' if c >= 0 else 'inverted')


def model_vs_truth(model, truth):
    a = gemmi.read_structure(str(model)); b = gemmi.read_structure(str(truth))
    for s in (a, b):
        s.remove_alternative_conformations(); s.remove_hydrogens(); s.remove_empty_chains()
    ca_a = [r.find_atom('CA', '*').pos for ch in a[0] for r in ch if r.find_atom('CA', '*') is not None]
    ca_b = [r.find_atom('CA', '*').pos for ch in b[0] for r in ch if r.find_atom('CA', '*') is not None]
    if not ca_a or not ca_b:
        return dict(built_ca=len(ca_a), truth_ca=len(ca_b), completeness=0.0)
    pa = np.array([[p.x, p.y, p.z] for p in ca_a]); pb = np.array([[p.x, p.y, p.z] for p in ca_b])
    nn = np.sqrt(((pa[:, None, :] - pb[None, :, :]) ** 2).sum(-1)).min(1)
    return dict(built_ca=len(ca_a), truth_ca=len(ca_b),
                completeness=round(len(ca_a) / len(ca_b), 3),
                ca_within_1a=round(float((nn < 1.0).mean()), 3),
                ca_within_2a=round(float((nn < 2.0).mean()), 3),
                ca_median=round(float(np.median(nn)), 2))


def campaign_endpoint(case_dir):
    """The endpoint the campaign selected: lowest work-set R, R-free never read."""
    root = case_dir / '00_sad_to_model'
    sel = root / 'endpoint_selection.json'
    if sel.exists():
        s = json.loads(sel.read_text())
        best = None
        for e in s.get('endpoints', []):
            if e.get('status') != 'completed':
                continue
            if best is None or (e['primary_metric']['r_work'] < best['primary_metric']['r_work']):
                best = e
        if best:
            return Path(best['output']), best['primary_metric']
    # fall back: any completed endpoint with the lowest work R
    best = None
    for d in sorted(root.glob('endpoint_*')):
        p = d / 'pipeline_report.json'
        if not p.exists():
            continue
        r = json.loads(p.read_text())
        if r.get('status') != 'completed':
            continue
        m = r.get('primary_metric') or {}
        if m.get('r_work') is None:
            continue
        if best is None or m['r_work'] < best[1]['r_work']:
            best = (d, m)
    return best if best else (None, None)


def main():
    cases = []
    for case_dir in sorted(CAMPAIGN.iterdir()):
        if not case_dir.is_dir():
            continue
        pdbid = case_dir.name
        run, metric = campaign_endpoint(case_dir)
        src = source_dir(pdbid)
        truth = None
        fasta = None
        if src:
            for cand in (src / f'{pdbid.lower()}.pdb', src / f'{pdbid}.pdb'):
                if cand.exists():
                    truth = cand
            f = src / f'rcsb_pdb_{pdbid.upper()}.fasta'
            fasta = f if f.exists() else None
        entry = dict(case=case_dir.name, pdb=pdbid, campaign=str(case_dir),
                     run=str(run) if run else None, metric=metric,
                     truth=str(truth) if truth else None, fasta=str(fasta) if fasta else None,
                     source=str(src) if src else None)
        if run:
            rep = json.loads((run / 'pipeline_report.json').read_text())
            entry['config'] = {k: rep['config'].get(k) for k in
                               ('method', 'copies', 'elements', 'wavelength', 'dm_solvent',
                                'dm_cycles', 'phase_cycles', 'rounds', 'cycles', 'feedback_cycles',
                                'building_recipe', 'refinement_prior', 'dm_engine', 'trust_wang_radius')}
            for key, fname in (('native', '01_inputs/native.mtz'), ('sad', '01_inputs/selected_sad.mtz'),
                               ('heavy', '02_sad/refined_heavy.pdb')):
                p = run / fname
                entry[key] = str(p) if p.exists() else None
            dm = run / '03_dm' / 'map_phases.mtz'
            entry['dm'] = str(dm) if dm.exists() else None
            model = run / 'primary_model.cif'
            entry['model'] = str(model) if model.exists() else None
            entry['data'] = mtz_meta(run / '01_inputs/native.mtz')
        if truth:
            entry['truth_stats'] = structure_stats(truth)
        cases.append(entry)
        q = entry.get('build_quality') or {}
        print(f"{pdbid:5s} {entry.get('data',{}).get('dmin','?'):>5} A  "
              f"R={metric['r_work'] if metric else '?':>6}/{metric['r_free'] if metric else '?':>6}  "
              f"solv={(entry.get('truth_stats') or {}).get('solvent','?')}", flush=True)
    OUT.write_text(json.dumps(cases, indent=1))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
