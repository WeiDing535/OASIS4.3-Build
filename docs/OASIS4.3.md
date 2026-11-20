# OASIS4.3 — User & Technical Guide

OASIS4.3 is a pure-CCP4 SAD phasing pipeline: sequence → phases → density map → model. It chains SAD phasing, density modification (DM) and optional automated Buccaneer / REFMAC building into one fixed protocol. Version `4.3.8.6` (internal PEP440 `4.3.8+refmac-timeout-isolation`).

All external engines are user-supplied and selected via `--ccp4-root`: `phaser`, `cparrot`, `cbuccaneer`, `refmac5` from one CCP4 installation. Density modification uses CCP4 Parrot; trace/Fcalc feedback uses REFMAC zero-cycle (`NCYC 0`) structure factors. The source distribution ships no third-party runtime. Python 3.12 only; gemmi 0.7.5.

## Method paths

Each phasing method is paired with a fixed DM engine; method-specific parameters have no effect under other methods.

| Method | Phasing | Density modification | Building / refinement |
| --- | --- | --- | --- |
| `likelihood` | Continuous SAD | MLHL-best | Seedless Buccaneer + REFMAC observed-HL refinement |
| `pplus_bundled` | Direct-method P+ with bundled Phaser SAD likelihood | Parrot | Same |
| `pplus_trust` | Direct-method P+ with imported raw Phaser SAD likelihood and KL trust region | Parrot | Same |
| `pplus_internal` | Internal anomalous-difference likelihood + KL trust-region P+ | Internal NumPy MLHL DM | Same |

- `likelihood`: continuous SAD; requires `--reference`.
- `pplus_bundled`: P+ using the bundled Phaser SAD likelihood; no prior likelihood file supplied.
- `pplus_trust`: imports a raw Phaser SAD pure-observation likelihood MTZ (FP/SIGFP, PHIB/FOM, HLA–HLD) from the same cell, space group, hand and heavy-atom model. No automatic Phaser run. The P+ update is a real, non-zero update, not a "no-update observation" control.
- `pplus_internal`: internal anomalous-difference likelihood + internal MLHL DM; no Phaser prior input and no Parrot float entry point, so it is cross-platform consistent (7E1D R-free 0.3325 on Linux vs 0.3409 on macOS). It additionally requires `--reference` (sha256-locked). Its internal likelihood and DM are approximations, not equivalent to a full Phaser SAD likelihood.

## SAD phasing flow

1. **Input verification** — SAD MTZ, raw native observations with original free set, heavy-atom PDB / mmCIF, sequence, cell / space-group / hand consistency, column types and finite values.
2. **Phasing** — the selected method produces phases and FOM. Outputs `02_sad/experimental.mtz` (PHIB/FOM/HL; P+ also PPLUS) and `02_sad/observations.mtz`, a byte-copy of the pure observation likelihood used as the fixed refinement prior.
3. **Density modification** — Parrot (or internal MLHL DM) produces `03_dm/map_phases.mtz` with DM phases / FOM and a per-cycle mean-FOM trace. The original free set is excluded.
4. **Map handoff** — `04_bridge/building.mtz`: native amplitudes plus `FWT = F × final DM FOM` (applied once) and PHWT. Free-reflection map coefficients remain absent at the initial handoff.
5. **Build / refine** — per-cycle Buccaneer building followed by REFMAC refinement (`05_build/`), fixed 5 cycles × 10 refinement cycles by default; each cycle reuses the original pure-observation HL.
6. **Finish** — single-copy cleanup, work-set weight calibration, post-refinement (`06_finish/`); the final model is the predeclared last cycle, never chosen by R-free.

f″ is computed theoretically from element + wavelength by Gemmi's Cromer–Liberman algorithm (`--auto-fpp`, default on, no manual entry). Manual `--fpp` is retained only to reproduce legacy records. Invalid, non-finite or out-of-range wavelengths are rejected; MTZ default wavelengths are not used.

## Endpoint selection protocol

The feedback grid is defined by `--feedback-endpoints` (e.g. `2,3`) and `--dm-solvents` (e.g. `auto,0.5`); an empty value means a single endpoint. Every endpoint runs the same protocol. Selection rule, recorded verbatim in `endpoint_selection.json`:

> *lowest final work-set R across the predeclared endpoint grid; R-free never read*

`endpoint_selection.json` lists all endpoints' `r_work` / `r_free`, the selected endpoint and the rule. R-work fits the working set; R-free is an in-protocol reference only — not comparable across methods or software, and not an end-to-end blind score.

## Fcalc normalisation (trace/Fcalc feedback)

After the 4.3.0 engine migration, trace/Fcalc feedback uses `refmac5 NCYC 0` to compute FC/PHIC without moving coordinates. REFMAC FC carries an overall scale ≈1.56 (~1.66× legacy fmodel amplitudes), which over-weights the truncated model in Parrot and can flip feedback from positive to negative (7E1D R-free 42.5% → 59.17%).

`parrot_dm._attach_fc()` normalises resolution-dependently using REFMAC's `Overall B` from the log:

```
FC_normalised = FC × exp(B · s² / 4)
```

B < 0, so high-resolution terms are mainly down-weighted. This restores the legacy fmodel scale (per-reflection Pearson 0.99995, median amplitude ratio 0.996). `fc` is the default; `fc_all` / `refined` select FC_ALL / FC_ALL_LS. Parse failures fall back to un-normalised behaviour and record metadata; `parrot_with_model.mtz` stores the chosen column, B value and normalisation flag under `fcalc_engine`.

## Inputs

- **SAD MTZ**: mean amplitudes + sigmas and anomalous differences + sigmas, or full F+/F− Friedel pairs. No phases needed.
- **Raw-observation MTZ**: `F` (F type), `SIGF` (Q), `FreeR_flag` (I); original free set preserved; no FWT / PHWT / FOM / SAD_HL columns. May be the same file as SAD if columns are complete.
- **Heavy-atom PDB / mmCIF**: existing sites with correct hand / setting; no search or re-indexing.
- **Sequence**: standard 20-amino-acid FASTA / PIR, one record per chain; ASU copies supplied by the user.
- **f″**: computed automatically (`--auto-fpp`).
- **Output directory**: must not already exist.

SAD and raw observations must share cell, space-group Hall setting and merged reciprocal ASU indexing; the original free set is matched by exact HKL and never regenerated.

## Key parameters

| Group | Option | Default / range |
| --- | --- | --- |
| Method | `likelihood` / `pplus_bundled` / `pplus_trust` / `pplus_internal` | `likelihood` |
| Sample | ASU copies of the whole sequence file | 1; 1–100 |
| Sample | Solvent fraction | auto; else strictly 0–1 |
| Phasing | Iteration cycles | 5; 1–100 |
| DM | Cycles | 5; 1–100 |
| Build / refine | Enable Buccaneer / REFMAC | on; can be disabled |
| Build / refine | Rounds | 5; 1–20 |
| Refine | Cycles per round | 10; 1–200 |
| P+ | AOE | 0.5; strictly 0–1 |
| P+ | KL limit | 0.0005 by default |
| Common | kmin | 0.03; 0–1 |
| Common | Neighbour cap | 600; 0–10000 |
| Common | Per-target relation cap | 96; 0–4096 |
| Runtime | Random seed | 20261023; 0–2147483647 |
| Free set | Original free-set value | 0; 0–999 |
| Columns | SAD mapping | auto; `f=FP,sigf=SIGFP,dano=DANO,sigdano=SIGDANO` or `fplus=F(+),sigfplus=SIGF(+),fminus=F(-),sigfminus=SIGF(-)` |
| Resources | CCP4 root | `--ccp4-root` (GUI auto-detect) |
| Resources | DM reference MTZ | bundled reference; generic, not a target structure |
| Engines | Density modification | `--dm-engine parrot` (default, external CCP4 Parrot) or `handdm` (in-house) |
| Engines | Refinement | `--refine-engine refmac` (default, external CCP4 REFMAC) or `minirefmac` (in-house) |
| Engines | In-house engine paths | `--handdm-script`, `--handdm-hist-lib`, `--minirefmac-root`; blank searches the standard deployment locations |

Relation caps of 0 may increase time and memory; defaults 600 / 96 are finite search settings, not exhaustive computation. Fixed settings: Buccaneer 3 cycles in round 1 and 2 thereafter; single thread; no protein seed in round 1; final model from the fixed last cycle.

### Engine substitution

The external CCP4 engines remain the default. The in-house engines are opt-in and never take over a run silently: they are selected by name, and every record states which engine produced the result.

| Stage | External (default) | In-house | Switch |
| --- | --- | --- | --- |
| Density modification | CCP4 Parrot (`cparrot`) | `handdm` | `--dm-engine handdm` |
| Refinement | CCP4 REFMAC (`refmac5`) | `minirefmac` | `--refine-engine minirefmac` |

Both in-house engines run in the OASIS interpreter, so they need no separate installation; `minirefmac` carries its own monomer-restraint dictionary and needs no CCP4 monlib.

Protocol differences are refused explicitly rather than silently weakened:

- `handdm` cannot combine a partial model with its map, so `--feedback-cycles` must be `0`. A non-zero request is rejected at validation time; it does not run a loop whose result would differ from what the option says.
- `minirefmac` has no equivalent of the audited fixed-phase-prior refiner. With `--refinement-prior dm_phifom` it uses the DM phase/FOM columns as an MLHL phase prior instead. That is an approximation, not the same protocol: the in-house engine refines the phases with the prior, while the REFMAC route restrains them. The refinement report records `phase_prior.mode: mlhl_prior`, `approximation: true` and `equivalent_to_fixed_phase_restraint: false`.
- The REFMAC `WEIGHT MATRIX` pilot reads geometry RMSD; `minirefmac` takes a direct xray/geometry balance coefficient, so the pilot is skipped, the engine's declared default restraint weight is used, and `SELECTION.json` records `weight_calibration: not_supported_by_engine`.
- `--per-round-cleanup` installs its wrapper on the external refiner, so it is rejected with `--refine-engine minirefmac`.
- In the `pipeline` building recipe the later rounds build from FWT/PHWT only, because the in-house engine does not emit the combined HL columns that the REFMAC route produces; `hybrid_report.json` lists the differences. The adapter derives those columns (`PHCOMB` and `HLACOMB`..`HLDCOMB`) from the refined phase and figure of merit, because CCP4 Buccaneer's correlation mode cannot build its map without them; the derivation is recorded as `combined_phase_quartet` in the refinement report.

Phasing (CCP4 Phaser) and model building (CCP4 Buccaneer) still require the user-supplied CCP4 installation.

## Emitted files

The `--output` directory (`00_sad_to_model/` in the ten cases):

| Path | Contents |
| --- | --- |
| `01_inputs/` | SAD verification, original free set; expanded modelling sequence when building is on |
| `02_sad/experimental.mtz` | Pre-DM PHIB/FOM/HL (P+ also PPLUS) |
| `02_sad/observations.mtz` | Pure observation likelihood HL and moments; fixed refinement prior |
| `03_dm/map_phases.mtz` | DM phases, FOM and metadata |
| `04_bridge/building.mtz` | Native amplitudes + `FWT = F × final DM FOM` + PHWT |
| `05_build/` | Per-round Buccaneer, REFMAC, refined maps and logs |
| `06_finish/` | Cleanup, work-set weight calibration, post-refinement; final model copy |
| `primary_model.cif` | Final model of the fixed last cycle (not R-free-selected) |
| `pipeline_report.json` | Auditable record; `status: "completed"` means all planned stages finished and hashes verified, otherwise check `error` |
| `endpoint_selection.json` | Multi-endpoint grid and selection rule |

Three MTZ classes are kept strictly separate: the original phasing result (`experimental.mtz`), the pure observation likelihood (`observations.mtz`) and the DM posterior (`map_phases.mtz`). REFMAC uses `observations.mtz` every round; it is never replaced by phasing, DM or model posteriors. With `--no-build-model`, only `01_inputs`, `02_sad`, `03_dm` and the pipeline record are produced; `primary_model` is `null` and `dm_phases` is the final phase file — a successful completion, not a building failure.

## GUI

Native Tk desktop app (`OASIS4.3.command` / `OASIS4.3.app`, or `python -B run_desktop.py --project .`): no web server, no port, no network — data stays local. Runs on macOS and Linux (headless: `xvfb-run -a python -B run_desktop.py --project .`; Apple Silicon needs Rosetta).

- Two full-width tabs: "Configuration / Input" and "Run / Results". Actions auto-switch tabs; the mouse wheel is isolated per page.
- "Check configuration" is read-only. Advanced parameters are hidden by default behind "Show advanced parameters".
- "Refinement results" plots R-work (solid) / R-free (dashed) per cycle; the curves and table share one per-cycle report. "Save final model" and "Export results" copy the full output directory. "3D structure" shows the Cα trace in a Canvas (rotate / zoom, optional all atoms).
- "Job history" persists runs in `<project>/.gui-oasis/jobs`; "Reuse parameters" refills the form with a new output directory; imported results are view-only. Up to four jobs run concurrently, each with its own output directory, log and subprocess group.
- CLI equivalent: `./OASIS4.3.command run_sad_to_model.py --help`.

## Resume and version retention

Stopping cleans up subprocesses but keeps stage files. CLI resumes with the same arguments plus `--resume`, verifying inputs, resources, environment, code and completed-file hashes before skipping finished stages; interrupted logs are saved as `<stage>.interrupted-<id>.log`. The GUI has no resume button. A new output directory is required when changing the build switch, algorithm, parameters, code or inputs; checkpoints are never mixed across versions.

## Validation scope and limitations

Validated on ten cases under a fixed protocol: 3O2E, 5D9D, 5HHK, 6HF7, 6ZPP, 7E1D, 7ENM, 7TBS, 7TFQ, 7TJ1. "10/10" means all ten complete the full pipeline with R-free within ±5 pp of an independent reference; some absolute R-free values remain high (e.g. 3O2E 42.64%, 7E1D 42.33%). Defaults are not the validation protocol — the ten-case settings were tuned per case (e.g. 7TFQ uses `pplus_trust` with observational priors).

- 4.3.8: a single diverging REFMAC pilot in finish no longer collapses the case; pilots are isolated and the timeout cut from 1800 s to 600 s; 6ZPP now completes (R-free 32.78%).
- 4.3.7: Fcalc amplitude normalisation restored 7E1D R-free 59.17% → 41.67%.
- 7E1D on the bundled path differs by ~15 pp between Linux and macOS (a last-digit Phaser SAD float difference amplified by Parrot near-singularity), so 41.67% is **not** a cross-platform hard metric; the `pplus_internal` path is consistent (0.3325 vs 0.3409).
- FOM is an internal confidence value, not true phase accuracy. SAD/DM exclude the original free set; downstream external map handoff is not certified work-set-only, so R-free is not an end-to-end blind metric. Cα previews do not replace chemical / structural validation. This is a research implementation requiring experimental validation.
