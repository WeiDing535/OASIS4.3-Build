# Changelog

## Unreleased

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
- Separates the P+ engine from CCP4: the form gains 「包内密度引擎（density/）」,
  `run_sad_to_model.py` gains `--density-runtime`, and a bare CCP4 now either
  borrows a reachable `density/` tree or is switched to the package runtime —
  with a message naming `density/bin/oasis-density-python` instead of failing in
  Phaser. `docs/INSTALL_BUNDLE.md` now states why the runtime cannot be
  redistributed and what a user with only CCP4 can do.
- Stops leaking the calibrated P+ recipe into the other methods: `collect()`
  drops `protocol`, the DM phase-restraint prior and the trust-region switches
  unless the method is `pplus_bundled`.

## 0.8.3 — public source export, 2026-10-22

- Publishes the recovered PhaseBuild source as OASIS4.3-Build.
- Includes local browser GUI, native research modules and explicit external-engine integration.
- Retains seeded and unseeded construction with fixed original HL and declared original free values.
- Provides 186 regression cases, source packaging and launch instructions.
- Leaves the existing OASIS4.3-GUI repository and v0.1 release unchanged.
