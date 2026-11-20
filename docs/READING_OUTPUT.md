# Reading Pipeline Output

Applies to `run_sad_to_model.py` output directories (e.g. `work/reproduce/<case_id>`, or `work/reproduce/<PDB>` written by the reproduce scripts). Follow this page top to bottom: first decide "did it complete", then "is it trustworthy".

## 1. Output directories by stage

| Directory | Stage | Key artifacts | Check first on failure |
| --- | --- | --- | --- |
| `01_inputs/` | Input parsing and free-set verification | `sad.mtz`, `sad.json`, `building.fasta`, `native.mtz` | `sad.json`: correct SAD columns, space group, free-set flag |
| `02_sad/` | P+/SAD phasing; phasing result and observed HL saved separately | `experimental.mtz(.json)`, `observations.mtz(.json)` | `experimental.json` FOM and site hand; `observations.mtz` is the observation-only likelihood handed to refinement |
| `02_hand/` | Optional: automatic hand determination (only when not bundled and `--auto-hand` is set) | `selected_sad.mtz`, `selected_native.mtz`, `refined_heavy.pdb` | Whether hand/sites were filtered correctly |
| `03_dm/` | Density modification, producing model-building phases | `map_phases.mtz`, `map_phases.json` | `map_phases.json`, Parrot/cparrot logs; the `dm_phases` fingerprint points here |
| `04_bridge/` | Generate work-set building map by HKL (free set masked) | `building.mtz`, `building.audit.json` | `building.audit.json`: HKL consistency and free-set masking audit |
| `05_build/` | Buccaneer building + refinement with fixed observed HL | `primary_model.cif`, `hybrid_report.json` | `hybrid_report.json` per-cycle R, Cα count, `refmac.log/xml` |
| `06_finish/` | Finalization/endpoint cleanup (pplus_trust / pplus_internal / pplus_bundled) | `primary_model.cif`, `endpoint.json` | Whether cleanup was accepted, final refinement geometry |
| output root | Summary layer | `pipeline_report.json`, `primary_model.cif`, `endpoint_selection.json` | Overall status and endpoint-selection conclusion |

> Final location of `primary_model.cif`: bundled/trust/internal in `06_finish/`, other methods in `05_build/`;
> the root usually also holds a copy for direct GUI/summary reading. Trust the `primary_model` fingerprint in `pipeline_report.json`.

## 2. Key fields in `pipeline_report.json`

| Field | Meaning | How to read |
| --- | --- | --- |
| `status` | Overall status: `completed` / `partial` / `failed` | Only `completed` counts as success |
| `primary_metric.cycle` | Final refinement cycle used | Fixed protocol usually takes the last cycle (e.g. 5), not reselected by R-free |
| `primary_metric.r_work` / `r_free` | Work/free R factors | Typical magnitude ~0.2–0.6; a large gap suggests overfitting or data problems |
| `primary_metric.bond_rmsd_angstrom` | Bond-length RMSD (Å) | Finite and ~0.005–0.02 Å; `null` means REFMAC geometry was not read |
| `primary_metric.angle_rmsd_degrees` | Bond-angle RMSD (°) | Finite and ~1–3°; `null` as above, check `refmac.log/xml` |
| `primary_model` | sha256 fingerprint of the final model | Must match the hash of on-disk `primary_model.cif` to be complete |
| `dm_phases` | sha256 fingerprint of `03_dm/map_phases.mtz` | Confirms the map was not replaced; the file must exist |
| `experimental_prior` / `pre_dm_result` / `inputs` / `environment` | Fingerprints of observed prior, phasing result, inputs and environment | For reproduction checks; cross-machine numeric identity is not required |

> `bond_rmsd_angstrom` / `angle_rmsd_degrees` may be `null` in intermediate endpoint-scan reports — this is normal;
> only investigate the geometry logs if still `null` in the final completed state.

## 3. `endpoint_selection.json` (endpoint selection)

| Field | Meaning |
| --- | --- |
| `rule` | `lowest final work-set R across the predeclared endpoint grid; R-free never read` |
| `endpoints[]` | `feedback`, `solvent`, `status`, `primary_metric`, `output` for each predeclared endpoint |
| `selected_feedback` / `selected_solvent` | Selected endpoint parameters |
| `selected_metric` | `r_work` / `r_free` / geometry of the selected endpoint |

**Rule: select only by work-set R, lowest to highest; never read R-free at any point.** Do not manually substitute an endpoint because its R-free is lower —
that breaks the predeclared protocol and introduces selection bias.

## 4. Did it succeed (all three must hold)

| # | Criterion | On failure |
| --- | --- | --- |
| 1 | `status == "completed"` in `pipeline_report.json` | Check the `error` field and the last-stage log |
| 2 | `primary_model.cif` exists, non-empty, and its hash equals the `primary_model` fingerprint | Missing/hash mismatch → artifact truncated or replaced, rerun |
| 3 | `r_work` / `r_free` ~0.2–0.6 with no large gap, geometry RMSD finite and reasonable | R > 0.6 or diverging geometry → treat as model-building failure, do not cite as a success |

Auxiliary checks: the `03_dm/map_phases.mtz` referenced by `dm_phases` exists; when `endpoint_selection.json` is present, its
`selected_metric` matches `primary_metric`.

## 5. Common failure stages and where to look

| Symptom | Common cause | Where to look |
| --- | --- | --- |
| Stops at `01_inputs` | Wrong MTZ column names/space group/free flag; sequence vs ASU copy count mismatch | Check `sad.json`; confirm SAD columns F+/F−; verify `--copies`, `--elements`, `--wavelength` |
| Stops at `02_sad` | Wrong heavy-atom sites/hand; unreasonable f'/f''; P+ not converging | Check `experimental.json` FOM and the `sg.pdb` cell; verify `--fpp` or `--auto-fpp` |
| Stops at `03_dm` | `cparrot` missing or crashed; wrong solvent content | Confirm `runtime/bin/cparrot` is executable and `--dm-executable` is correct; check `map_phases.json`, the Parrot log and `--dm-solvent` |
| Stops at `04_bridge` | Inconsistent HKL sets, free set not masked | Check `building.audit.json`; do not continue building with a map that failed the audit |
| Stops at `05_build` | Empty build/diverging refinement; external engine missing | Check `hybrid_report.json` per-cycle R and Cα count, `refmac.log/xml`; confirm `--rounds`, `--cycles` and the runtime |
| Some endpoints fail in `06_finish` | Some endpoint protocols did not converge | Check each `endpoints[].status` in `endpoint_selection.json`; at least one `completed` is required to be selected |
| `status != completed` but a model file exists | Interrupted/partial write | Do not trust it; locate the cause via the `error` field and last-stage log, and use a new output directory before rerunning |

## 6. Important: R-free is not an end-to-end blind test

- R-free here comes from the free set flagged in the input MTZ; phasing, DM, building and refinement share the same observations and map handoff
  (`observations.mtz` / `map_phases.mtz`). It reflects consistency and overfitting trends, **but it is not an external blind test independent of the whole pipeline**.
- Do not use R-free to pick endpoints, cycles or models; selection follows §3, by work-set R.
- Low R factors, high FOM, many residual bonds — no single metric alone proves a structure is correct; check sequence coverage, geometry
  and the work/free gap together, and perform an independent external validation when needed.
