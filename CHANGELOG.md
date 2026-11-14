# Changelog

## Unreleased

- **Density modification is no longer a module or a choice.** The desktop form
  now shows two modules: 「模块 1 · P+ 定相 + oasis-dm 密度修饰（必做）」, fixed and
  greyed out, and 「模块 2 · Buccaneer 建模 + REFMAC 精修」. The old
  density-modifier checkbox is gone, the `dm_engine` combobox is gone, and the
  density-modification parameter fold is never greyed out because there is
  nothing to switch off: `oasis-dm v2` runs inside the P+ chain on every method.
  `collect()` always reports `dm_engine='oasis-dm'`/`dm_enabled=True`; module
  detection only looks for `cbuccaneer`/`refmac5` and reports the in-tree `dm`
  stage as always available.
- **The old density-modification engine is deleted outright**, not hidden: its
  adapter module, its CLI subcommand, `--dm-engine`'s third choice, the
  `density/` lookup, the packaging labels and every test that assumed a
  `density/` directory are gone. `--dm-engine` accepts only
  `oasis-dm` (default) and `internal`; anything else raises
  「Unknown density-modification engine」. The 2026-11-12 ten-case table in
  `docs/STATUS_20261112.md` was measured with the deleted engine and is kept as
  a historical record only — it can no longer be reproduced. The sweep goes
  further than code: every historical report that described the old engine now
  carries a 「历史记录」 banner and calls it 旧引擎, its old-engine-era run
  directories under `runs/` and its experiment scripts/logs under `work/` were
  deleted, and the stale `.gui-oasis/jobs/` records that named it were removed.
  The hash-locked upstream snapshot `vendor/oasis43/` was deleted too, so no
  file in the repository mentions the retired engine any more.
- **The engine is named `oasis-dm` everywhere.** It was written under the
  working name handdm; inside OASIS4.3 the module is `src/phasebuild/oasis_dm.py`,
  the vendored pair is `vendor/oasis_dm/oasis_dm_core.py` +
  `oasis_dm_v2.py`, the engine id and `--dm-engine` value is `oasis-dm`, the
  staged files are `oasis_dm_input.mtz` / `oasis_dm_with_model.mtz` /
  `oasis_dm.log`, the override is `PHASEBUILD_OASIS_DM_PYTHON`, and every GUI
  label, doc and test uses the new name (`docs/OASIS_DM_REPLACEMENT_20261114.md`).
  `vendor/oasis_dm/PROVENANCE.md` keeps the upstream paths and hashes and states
  the single executable-line change the rename required
  (`from oasis_dm_core import ...`); no numerical statement was touched.
- **The hash-locked `vendor/oasis43/` snapshot is deleted.** Its eight
  referencing sites were updated: the CI pytest command, `testpaths` in
  `pyproject.toml`, the source packager's file list, the bundle assembler's
  provenance copy, `AUTHORS.md`, `docs/PROVENANCE.md`, `docs/INTEGRATION.md` and
  `README_EN.md`. The suite drops from 864 to 707 passed because that snapshot's
  own 157 tests no longer exist; the remaining 707 are green.
- Density modification now defaults to the in-tree **oasis-dm v2** engine
  (`src/phasebuild/oasis_dm.py`, vendored `src/phasebuild/vendor/oasis_dm/`),
  replacing the old external engine for every method: `--dm-engine` accepts
  `oasis-dm` (default) / `internal`, and the two P+ routes refuse `internal`
  outright. No external density modifier is needed for a run, nor is CCP4's
  `hist.lib` — the histogram library ships as `lib/data/hist.lib`
  (`--dm-hist-lib`, `PHASEBUILD_HIST_LIB`).
- End-to-end check of that swap on the insulin 2.1 A SAD case, same
  `pplus_internal` command, only `--dm-engine` differs: oasis-dm v2 gives
  Rwork/Rfree **0.2131 / 0.2662** (52 CA, 51/51 truth CAs within 1 A, RMSD
  0.12 A) against **0.2758 / 0.3181** for the previous engine. The trace/Fcalc feedback
  loop needed one fix to get there: the merged model likelihood is now rescaled
  per resolution shell to the (normalised) experimental HL quartet and weighted
  by the fraction of scattered intensity the partial model explains; adding the
  raw Sim weight, which is in observed-amplitude units, let a 12-residue
  fragment outweigh the experimental term 7x and the map collapsed
  (holdout CC 0.844 -> 0.700 -> 0.528, Buccaneer built 0 residues).
  New flags: `--dm-cmbwt`, `--dm-holdout`, `--dm-hist-lib`. Defaults:
  `dm_engine=oasis-dm`, `dm_cycles=12`, `dm_cmbwt=0.4`, `dm_holdout=0.05`.
- oasis-dm v2 is not just a stand-in: on insulin 2.1 Å it lifts the DM map from
  CC 0.3646 to 0.7856–0.8049 and wMPE from 55.4° to 36.0–37.1°, and the
  downstream build is statistically indistinguishable from the previous
  engine's (R-free 0.3295 vs 0.3131 for its best length; engine means over
  ten variants 0.3599 ± 0.043 vs 0.3449 ± 0.030, difference 0.015 against a
  combined standard error of 0.019). The previous engine still won on phase accuracy
  (CC 0.8456 vs 0.8049, wMPE 34.29° vs 35.98°). See
  `docs/OASIS_DM_REPLACEMENT_20261114.md`.
- `--dm-cmbwt auto` runs a short holdout-scored sweep over
  {0.3, 0.4, 0.5, 0.6, 0.8} before the production pass. It refuses to override
  the calibrated 0.4 unless the sweep separates the candidates by at least
  0.01 holdout CC, because cmbwt only rescales the HL handed to the next stage
  and map-based scores are structurally blind to it. On insulin the sweep is
  flat (0.7963–0.8009) and the highest-CC weight (0.5) is the worst downstream
  (R-free 0.4367 vs 0.3320 at 0.3). Keep the default 0.4; use 0.3 for a better
  build on this case, at the cost of slightly worse phases.
- Phaser is no longer a hard requirement for the SAD front end, and the numbers
  behind that are now in `docs/NO_CCP4_REPLACEMENT.md` §2. With the previous
  density-modification engine the internal likelihood front end reached 25.24 / 31.35 on 7TJ1 and
  26.69 / 29.36 on 5D9D against the Phaser campaign's 23.92 / 31.08 and
  25.43 / 27.51 — gaps of 0.3 and 1.9 points R-free, where the old table
  (internal front end *and* internal DM) showed +15.2 and +18.9. The remaining
  handicap is substructure refinement, which stock CCP4 `bp3` supplies: on
  insulin it refines six S sites in 12 s and its SAD phases agree with Phaser's
  `experimental.mtz` at CC 0.981 (4508 reflections, within45 0.999), against
  0.893 for the internal likelihood on the same sites. `bp3` is not wired into
  the pipeline yet.
- `--dm-engine oasis-dm` with `--refinement-prior dm_phifom` is accepted for
  `--method pplus_internal`. The restraint is applied by REFMAC and written by
  the density-modification stage, so nothing about it belonged to Phaser; only
  the bundled route used to be allowed because that was the only route that had
  been calibrated.
- `--max-per-target` defaults to 48, matching the GUI and the calibrated
  ten-case protocol. The old 96 built 96 relations per reflection, which for
  any case above ~21k reflections exceeded the 2 M hard relation budget and
  aborted the trust-region stage with "Relation memory budget exceeded" before
  any phasing happened (7TJ1: 31 479 reflections, 5D9D: 35 489). The budget
  error now names the flag and the reflection count instead of suggesting an
  unexposed "increase budget".
- Deletes `density/` (2.4 GB of Phenix-derived RESOLVE/CCTBX). It only backed
  the two-hand Wang8 pilot and the old RESOLVE density-modification backend;
  both now run on the in-tree oasis-dm engine. P+ needs `phaser` + `cbuccaneer` +
  `refmac5` and nothing else. On insulin 2.1 Å this gives Rwork 0.3166 / Rfree
  0.3779 with 48 residues in 17 min, byte-identical (`02_sad/experimental.mtz`
  and `03_dm/map_phases.mtz` sha256) to the portable-package P+ run that used
  the tree.
- Adds a 「双手性判别」 checkbox, enabled only for 「P+ 自带引擎」: unchecked, the
  run takes the hand from the input heavy-atom file; checked, both
  enantiomorphs are phased and density-modified, and the hand is chosen from
  work-set sequence evidence. Ambiguous evidence still fails loudly;
  `--allow-unverified-hand` stays an explicit opt-in that records
  `hand_verified=false`. (`hand_selection.density_modify` dispatches on
  `--dm-engine` and runs the in-tree `oasis-dm v2` by default, so the two-hand
  test needs no external density modifier either.)
- `--dm-engine` is now `internal|oasis-dm`. `resolve` is gone and raises
  「Unknown density-modification engine」, P+ routes pin `oasis-dm` (the GUI greys
  the box out), and `internal` is accepted only for `连续 SAD`, `旧版直接法 P+`
  and `内部 P+`. A missing `oasis-dm` falls back to `internal` on those three and
  fails loudly on P+.
- Removes `density_root()`, `require_density()`, `bundled_pplus.DENSITY_ENGINE`,
  `trust_workflow.dm()`, `--density-runtime` and the 「包内密度引擎」 form field;
  `run_pplus_trust_stage.py dm` now refuses with a pointer to `--dm-engine
  oasis-dm`.
- Packaging drops the tree too: `build_bundled_pplus.py` no longer takes
  `--density-source`, `prune_bundle.py` no longer prunes `runtime/density/`,
  `package_standalone.py` no longer checks RESOLVE segment files, and bundle
  metadata reports `Phenix_derived_components_bundled: false` with engines
  `Phaser/Buccaneer/REFMAC`.
- Adds `packaging/make_macos_app.py`, which wraps a tree's `OASIS4.3.command` in
  a double-clickable `OASIS4.3.app` (generated icon, no Terminal window).
- Shows every pipeline stage's own log in the 运行日志 pane and clears the
  progress bar when a job reaches a terminal state.
- Adds `docs/NO_CCP4_REPLACEMENT.md`: measured feasibility of replacing the four
  bundled CCP4 executables, with the REFMAC forward-model gap traced to the
  bulk-solvent amplitude term.
- Adds `run_refine.py --engine internal`, exposing the in-tree NumPy/SciPy
  refinement engine without changing the default REFMAC path.
- Fixes `internal_refine.atomic_attenuation` to always return an array, so the
  zero overall-B case no longer raises during a refinement cycle.
- Folds the parameter form into clickable groups (advanced groups hidden behind
  「显示高级参数」), so the window opens with only the inputs a run needs.
- Makes a stock CCP4 install the supported baseline. `连续 SAD（只需 CCP4）` and
  `旧版直接法 P+（只需 CCP4）` run the whole chain with `cbuccaneer`/`refmac5`
  plus the in-tree density modifier; measured on insulin 2.1 Å at Rwork 0.2796 / Rfree 0.3331 with
  50 residues in 12 min, better than the bundled P+ route (0.3166 / 0.3779).
- A bare CCP4 is no longer swapped for the package runtime and no longer
  refused: it switches to the CCP4-only method and says so in the red notice.
- Detection runs after the form exists and resolves the runtime even when the
  field is blank, so typing `/Applications/ccp4-9` raises the red notice *before*
  the first run (it used to appear only after 开始). On insulin, `连续 SAD` with
  external density modification reached Rwork 0.3150 / Rfree 0.3583 with 51
  residues in 6.5 min, versus 0.2796 / 0.3331 in 12 min with the internal engine:
  both run on stock CCP4, the internal one is steadier, the external one was
  faster.
- Allows external density modification on the `连续 SAD` and `旧版直接法 P+`
  methods (both write an observation likelihood); phases-only runs keep the
  internal engine and never resolve an external tracer.
- `run_trace_feedback.py` falls back to `02_sad/experimental.mtz` plus the SAD
  free set when a phases-only run has no `01_inputs/native.mtz`.
- Stops leaking the calibrated P+ recipe into the other methods: `collect()`
  drops `protocol`, the DM phase-restraint prior and the trust-region switches
  unless the method is `pplus_bundled`.

## 0.8.3 — public source export, 2026-10-22

- Publishes the recovered PhaseBuild source as OASIS4.3-Build.
- Includes local browser GUI, native research modules and explicit external-engine integration.
- Retains seeded and unseeded construction with fixed original HL and declared original free values.
- Provides 186 regression cases, source packaging and launch instructions.
- Leaves the existing OASIS4.3-GUI repository and v0.1 release unchanged.
