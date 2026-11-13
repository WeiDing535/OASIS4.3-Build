# Changelog

## Unreleased

- P+ no longer requires `density/`. Only the two-hand Wang8 pilot and the old
  RESOLVE/CCTBX density-modification backend lived there; with `--dm-engine
  parrot --no-auto-hand` the P+ chain runs on `phaser` + `cparrot` +
  `cbuccaneer` + `refmac5` alone. On insulin 2.1 Å this gives Rwork 0.3166 /
  Rfree 0.3779 with 48 residues in 17 min, byte-identical (`experimental.mtz`
  and `map_phases.mtz` sha256) to the portable-package P+ run that does use
  `density/`. The GUI keeps the P+ method on a bare CCP4 when `phaser` and
  `cparrot` are present, and states in the red notice that the hand now comes
  from the input heavy-atom file and was not verified.
- `require_density()` is only called when the run really needs the tree, and a
  `density/` borrowed from a neighbouring portable tree is now printed instead
  of being used silently.
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
- The form gains 「包内密度引擎（density/）」 and `run_sad_to_model.py` gains
  `--density-runtime`, but a bare CCP4 is no longer swapped for the package
  runtime and no longer refused: it switches to the CCP4-only method and says so
  in the red notice. `--density-runtime` now only names an explicit `density/`.
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
