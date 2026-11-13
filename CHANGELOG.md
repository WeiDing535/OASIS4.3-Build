# Changelog

## Unreleased

- Phaser is no longer a hard requirement for the SAD front end, and the numbers
  behind that are now in `docs/NO_CCP4_REPLACEMENT.md` §2. With `--dm-engine
  parrot` the internal likelihood front end reaches 25.24 / 31.35 on 7TJ1 and
  26.69 / 29.36 on 5D9D against the Phaser campaign's 23.92 / 31.08 and
  25.43 / 27.51 — gaps of 0.3 and 1.9 points R-free, where the old table
  (internal front end *and* internal DM) showed +15.2 and +18.9. The remaining
  handicap is substructure refinement, which stock CCP4 `bp3` supplies: on
  insulin it refines six S sites in 12 s and its SAD phases agree with Phaser's
  `experimental.mtz` at CC 0.981 (4508 reflections, within45 0.999), against
  0.893 for the internal likelihood on the same sites. `bp3` is not wired into
  the pipeline yet.
- `--dm-engine parrot` with `--refinement-prior dm_phifom` is now accepted for
  `--method pplus_internal`. The restraint is applied by REFMAC and written by
  Parrot, so nothing about it belonged to Phaser; only the bundled route was
  allowed to use it because that was the only route that had been calibrated.
- `--max-per-target` defaults to 48, matching the GUI and the calibrated
  ten-case protocol. The old 96 built 96 relations per reflection, which for
  any case above ~21k reflections exceeded the 2 M hard relation budget and
  aborted the trust-region stage with "Relation memory budget exceeded" before
  any phasing happened (7TJ1: 31 479 reflections, 5D9D: 35 489). The budget
  error now names the flag and the reflection count instead of suggesting an
  unexposed "increase budget".
- Deletes `density/` (2.4 GB of Phenix-derived RESOLVE/CCTBX). It only backed
  the two-hand Wang8 pilot and the old RESOLVE density-modification backend;
  both now run on CCP4 `cparrot`. P+ needs `phaser` + `cparrot` + `cbuccaneer` +
  `refmac5` and nothing else. On insulin 2.1 Å this gives Rwork 0.3166 / Rfree
  0.3779 with 48 residues in 17 min, byte-identical (`02_sad/experimental.mtz`
  and `03_dm/map_phases.mtz` sha256) to the portable-package P+ run that used
  the tree.
- Adds a 「双手性判别」 checkbox, enabled only for 「P+ 自带引擎」: unchecked, the
  run takes the hand from the input heavy-atom file; checked, both
  enantiomorphs are phased and density-modified with `cparrot`, and the hand is
  chosen from work-set sequence evidence. Ambiguous evidence still fails loudly;
  `--allow-unverified-hand` stays an explicit opt-in that records
  `hand_verified=false`.
- `--dm-engine` is now `internal|parrot`. `resolve` is gone and raises
  「已移除」, P+ routes pin `parrot` (the GUI greys the box out), and `internal`
  is accepted only for `连续 SAD`, `旧版直接法 P+` and `内部 P+`. A missing
  `cparrot` falls back to `internal` on those three and fails loudly on P+.
- Removes `density_root()`, `require_density()`, `bundled_pplus.DENSITY_ENGINE`,
  `trust_workflow.dm()`, `--density-runtime` and the 「包内密度引擎」 form field;
  `run_pplus_trust_stage.py dm` now refuses with a pointer to `--dm-engine
  parrot`.
- Packaging drops the tree too: `build_bundled_pplus.py` no longer takes
  `--density-source`, `prune_bundle.py` no longer prunes `runtime/density/`,
  `package_standalone.py` no longer checks RESOLVE segment files, and bundle
  metadata reports `Phenix_derived_components_bundled: false` with engines
  `Phaser/Parrot/Buccaneer/REFMAC`.
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
  `旧版直接法 P+（只需 CCP4）` run the whole chain with `cparrot`/`cbuccaneer`/
  `refmac5` only; measured on insulin 2.1 Å at Rwork 0.2796 / Rfree 0.3331 with
  50 residues in 12 min, better than the bundled P+ route (0.3166 / 0.3779).
- A bare CCP4 is no longer swapped for the package runtime and no longer
  refused: it switches to the CCP4-only method and says so in the red notice.
- Detection runs after the form exists and resolves the runtime even when the
  field is blank, so typing `/Applications/ccp4-9` raises the red notice *before*
  the first run (it used to appear only after 开始). On insulin, `连续 SAD` with
  Parrot density modification reaches Rwork 0.3150 / Rfree 0.3583 with 51
  residues in 6.5 min, versus 0.2796 / 0.3331 in 12 min with the internal engine:
  both run on stock CCP4, the internal one is steadier, Parrot is faster.
- Allows Parrot density modification on the `连续 SAD` and `旧版直接法 P+`
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
