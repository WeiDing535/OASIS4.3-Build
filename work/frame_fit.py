"""Put the deposited model into the frame of each campaign DM map.

Two independent ambiguities have to be resolved:

  * origin: the DM map and the deposited model differ by a translation, found by
    FFT density correlation of FWT/PHWT against the model density;
  * symmetry mate: the density is invariant under the space group, so the ASU
    placement is only fixed up to a symmetry operation.  We test every operator
    and keep the one that puts the truth closest to the campaign build.

Output: work/tencase_frames.json + work/frames/<PDB>_truth_in_map_frame.pdb
"""
import json
from pathlib import Path
import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
WORK = ROOT / 'work'
FRAMES = WORK / 'frames'


def strip(x):
    st = x.clone() if isinstance(x, gemmi.Structure) else gemmi.read_structure(str(x))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    try:
        st.remove_ligands_and_waters()
    except RuntimeError:
        pass
    st.remove_empty_chains()
    return st


def density(st, dmin):
    calc = gemmi.DensityCalculatorX()
    calc.d_min = dmin
    calc.set_grid_cell_and_spacegroup(st)
    calc.put_model_density_on_grid(st[0])
    return np.asarray(calc.grid, dtype=np.float64)


def cc(A, B):
    A0 = A - A.mean()
    B0 = B - B.mean()
    return float((A0 * B0).sum() / np.sqrt((A0 * A0).sum() * (B0 * B0).sum()))


def move(st, shift, op=None, wrap=True):
    out = st.clone()
    cell = out.cell
    for ch in out[0]:
        for r in ch:
            for at in r:
                f = cell.fractionalize(at.pos)
                xyz = [f.x, f.y, f.z]
                if op is not None:
                    xyz = op.apply_to_xyz(xyz)
                xyz = [xyz[k] + shift[k] for k in range(3)]
                if wrap:
                    xyz = [x % 1.0 for x in xyz]
                at.pos = cell.orthogonalize(gemmi.Fractional(*xyz))
    return out


def ca_array(st, chain=None):
    return np.array([[at.pos.x, at.pos.y, at.pos.z] for ch in st[0]
                     for r in ch for at in r
                     if at.name == 'CA' and (chain is None or ch.name == chain)])


def nn(a, b):
    d = np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1))
    return d.min(1)


def quality(built, truth):
    a, b = ca_array(built), ca_array(truth)
    if len(a) == 0 or len(b) == 0:
        return None
    v = nn(a, b)
    return dict(built_ca=len(a), truth_ca=len(b),
                ca_within_1a=round(float((v < 1.0).mean()), 3),
                ca_within_2a=round(float((v < 2.0).mean()), 3),
                median=round(float(np.median(v)), 2))


def main():
    FRAMES.mkdir(exist_ok=True)
    inv = json.loads((WORK / 'tencase_inventory.json').read_text())
    out = []
    for c in inv:
        pdb = c['pdb']
        dm = gemmi.read_mtz_file(c['dm'])
        dmin = dm.resolution_high()
        truth = strip(c['truth'])
        built = strip(c['model']) if c.get('model') else None
        B = density(truth, dmin)
        A = np.asarray(dm.transform_f_phi_to_map('FWT', 'PHWT', exact_size=list(B.shape)),
                       dtype=np.float64)
        A0, B0 = A - A.mean(), B - B.mean()
        corr = np.fft.ifftn(np.fft.fftn(A0) * np.conj(np.fft.fftn(B0))).real
        idx = np.unravel_index(int(np.argmax(corr)), corr.shape)
        shape = B.shape
        shift = [i / s for i, s in zip(idx, shape)]
        best_shift, best_cc = None, -2.0
        for sign in (+1, -1):
            s = [sign * x for x in shift]
            val = cc(A, density(move(truth, s), dmin))
            if val > best_cc:
                best_shift, best_cc = s, val
        ops = list(gemmi.SpaceGroup(truth.spacegroup_hm).operations())
        cands = []
        for i, op in enumerate(ops):
            st = move(truth, best_shift, op)
            q = quality(built, st) if built is not None else None
            score = q['median'] if q else 0.0
            cands.append((score, i, op.triplet(), st, q))
        cands.sort(key=lambda x: x[0])
        score, i, triplet, st, q = cands[0]

        # Per-chain placement: a multi-copy build can assign its copies to
        # different symmetry mates than the deposited ASU, which is still the
        # same crystal and matches the map equally well.  Lattice translations
        # are kept so that atoms land where the reference build wrote them.
        chain_rows = []
        if built is not None:
            matches = []
            for bch in built[0]:
                bca = ca_array(built, bch.name)
                if len(bca) == 0:
                    continue
                best = None
                for tch in truth[0]:
                    if len(ca_array(truth, tch.name)) == 0:
                        continue
                    for j, op in enumerate(ops):
                        cand = move(truth, best_shift, op, wrap=False)
                        cca = ca_array(cand, tch.name)
                        if len(cca) == 0:
                            continue
                        fb = cand.cell.fractionalize(gemmi.Position(*bca.mean(0)))
                        ft = cand.cell.fractionalize(gemmi.Position(*cca.mean(0)))
                        lat = [round(fb.x - ft.x), round(fb.y - ft.y), round(fb.z - ft.z)]
                        cand = move(cand, lat, None, wrap=False)
                        m = float(np.median(nn(bca, ca_array(cand, tch.name))))
                        if best is None or m < best[0]:
                            best = (m, j, tch.name, cand, lat)
                matches.append((bch.name, len(bca), best))

            # one transform per truth chain: trust the longest matching build
            chosen = {}
            for bname, nca, best in matches:
                if best is None or best[0] >= 2.0:
                    chain_rows.append(dict(built_chain=bname, truth_chain=None,
                                           symmetry_op=None, lattice=None,
                                           median=None if best is None else round(best[0], 2)))
                    continue
                m, j, tname, cand, lat = best
                chain_rows.append(dict(built_chain=bname, truth_chain=tname,
                                       symmetry_op=j, lattice=lat, median=round(m, 2)))
                if tname not in chosen or nca > chosen[tname][0]:
                    chosen[tname] = (nca, m, j, cand, lat)
            st = move(truth, best_shift, ops[i], wrap=False)
            for tname, (nca, m, j, cand, lat) in chosen.items():
                src = [ch for ch in cand[0] if ch.name == tname][0]
                dst = [ch for ch in st[0] if ch.name == tname][0]
                for ra, rb in zip(src, dst):
                    for aa, ab in zip(ra, rb):
                        ab.pos = aa.pos
            q = quality(built, st)
        out_pdb = FRAMES / f'{pdb}_truth_in_map_frame.pdb'
        st.write_pdb(str(out_pdb))
        rec = dict(pdb=pdb, hand='same', cc=round(best_cc, 3), idx=list(map(int, idx)),
                   shape=list(shape), shift_fractional=[round(x, 4) for x in best_shift],
                   symmetry_op=i, symmetry_triplet=triplet, chain_map=chain_rows,
                   truth_in_map_frame=str(out_pdb),
                   build_quality_map_frame=q)
        out.append(rec)
        print(f"{pdb}: CC={best_cc:6.3f} shift={[round(x,3) for x in best_shift]} "
              f"op{i}={triplet} median={q['median'] if q else '?'} "
              f"<1A={q['ca_within_1a'] if q else '?'}", flush=True)
    (WORK / 'tencase_frames.json').write_text(json.dumps(out, indent=1))
    print('wrote', WORK / 'tencase_frames.json')


if __name__ == '__main__':
    main()
