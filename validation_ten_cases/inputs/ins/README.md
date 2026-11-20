# Insulin S-SAD Example Inputs (`ins/`)

A complete sulfur SAD (S-SAD) mini-example: space group I213, ~2.1 Å, wavelength 2.0 Å, two chains A/B with 51 residues total, and 6 cysteine SG sulfur atoms as heavy-atom sites. It reproduces the full "SAD observations → protein model" flow without external data.

## File → input parameter map

| File | Content | CLI parameter | Note |
| --- | --- | --- | --- |
| `ins_wave2.0_reso2.1A_sad.mtz` | I213, ~2.1 Å, SAD columns F+/F− (anomalous pair) | `--sad` | Main phasing input; column names auto-detected |
| `3INS.pir` | 51-residue sequence (chain A 21 + chain B 30), pir/FASTA header | `--sequence` | `--copies auto` infers ASU copy number from it |
| `sg.pdb` | 6 SG sulfur atoms (A6/A7/A11/A20, B7, B19) | `--heavy` | Heavy-atom/substructure input; I213 cell matches the MTZ |
| `ipcas_test_915.1_refmac1.pdb` | Reference model | External comparison only, **not passed to the command** | `pplus_bundled` neither reads nor requires `--reference`; do not pass it, to avoid contaminating phasing |

This example has no native MTZ; `--native` is left empty (the selected SAD observations are used by default), so `--free-value` uses the default 0. The free-set value is **not** an end-to-end blind test.

## Why element S, wavelength 2.0 Å

Sulfur's intrinsic K edge is around 2.07 Å (~5.99 keV). At 2.0 Å its *f''* is ~0.56 e — weak but measurable anomalous signal, a typical **S-SAD (intrinsic sulfur SAD)** case. Therefore:

- `--elements S --wavelength 2.0 --auto-fpp`: theoretical f'/f'' are computed offline from the explicit wavelength and element; no `--fpp` needed.
- The heavy-atom input is just the 6 SG atoms in `sg.pdb`; few sites and weak signal, so `--seed 20261023` is fixed with longer `phase-cycles 100 / dm-cycles 40` and a 5-round × 10-cycle build/refine budget.

## Full runnable command (paths relative to package root)

Run from the project root; `work/` is a local output directory ignored by `.gitignore`. OASIS bundles no CCP4 component — point `--ccp4-root` at a **user-supplied, licensed** install:

```bash
python3.12 -u run_sad_to_model.py \
  --sad validation_ten_cases/inputs/ins/ins_wave2.0_reso2.1A_sad.mtz \
  --sequence validation_ten_cases/inputs/ins/3INS.pir \
  --heavy validation_ten_cases/inputs/ins/sg.pdb \
  --output work/reproduce/ins_external \
  --ccp4-root /Applications/ccp4-9 \
  --dm-executable /Applications/ccp4-9/bin/cparrot \
  --method pplus_bundled \
  --seed 20261023 \
  --copies auto \
  --phase-cycles 100 --dm-cycles 40 --rounds 5 --cycles 10 \
  --dm-engine parrot --dm-solvent 0.5 --dm-mask-radius 3.0 --dm-anisotropy \
  --feedback-cycles 2 --feedback-endpoints 2,3 --trace-resolution 2.0 \
  --elements S --wavelength 2.0 --auto-fpp \
  --build-model
```

Requires Python 3.12 + gemmi 0.7.5 and a user-supplied CCP4 (`phaser` / `cparrot` / `cbuccaneer` / `refmac5` plus `lib/data`). The CCP4 path may also be given via `CCP4_ROOT`; the interpreter via `PHASEBUILD_PYTHON`.

## Expected results (measured locally)

`--feedback-endpoints 2,3` declares two fixed-protocol endpoints in parallel, then selects **by ascending R-work** (R-free is not read):

| Endpoint | R-work | R-free |
| --- | ---: | ---: |
| fb2 (2 trace-feedback rounds) | 0.3166 | 0.3779 |
| fb3 (3 trace-feedback rounds) | 0.2141 | 0.2343 |

The R-work rule selects **fb3**, so `primary_model.cif` corresponds to fb3. Final model: **51 residues, 401 atoms**.

> These numbers are local measurements tied to `--seed 20261023` and all parameters above; changing machine, runtime version, seed or any hyperparameter changes them. Result-reading guide: `docs/READING_OUTPUT.md`.
