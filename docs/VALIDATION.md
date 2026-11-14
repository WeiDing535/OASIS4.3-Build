# Current validation pointers

Current research workspace: 417 passed, one reference-data test deselected.
Eight additional [dual-channel and conditional-branch checks with fixed-budget experiments](DUAL_CHANNEL_PPLUS_20261024.md) bring the current total to 417. Three previous [explicit DM-prior handoff tests and 18 matched-prior trials](DM_PRIOR_AUDIT_20261024.md) correct the previous external-engine attribution. Twelve targeted P+ checks and the [three-case experimental validation with DM-engine follow-ups](PPLUS_MULTICASE_VALIDATION_20261024.md) cover circular moments, float32 MTZ persistence, centric singletons, calibration reuse and HKL-order alignment. The prior [controlled insulin trials](PPLUS_OPTIMIZATION_20261023.md) remain historical evidence. The GUI and default product protocol remain 4.3.0.dev2; the prior product-only regression result was 394 passed.
The [insulin diagnosis](INSULIN_DIAGNOSIS_20261023.md) records independent hand controls, matched-model prior ablations, the observation-only handoff fix, and fixed-budget real-data results.

# Publication validation — 2026-10-22 (historical baseline)

For the 2026-10-23 recovery baseline, current regression results, limited line
coverage measurement and actual external-engine smoke results, see
[the development recovery report](../development/RECOVERY_REPORT_20261023.md).
The 186-case result below describes the original publication only.

The public source export passed 186 tests (0 failures, 0 errors, 0 skipped) in 20.67 seconds on the existing macOS Intel environment: Python 3.9.18, Gemmi 0.7.5, NumPy 1.19.5, SciPy 1.7.3, pytest 6.2.5.

Command: `python -B -m pytest -q -p no:cacheprovider --basetemp validation/tmp --junitxml validation/tests.xml tests`.

All 47 Python modules and GUI assets match the tested local trial copy byte-for-byte. Source files parse successfully, the launcher passes `bash -n`, source packaging passes its CRC check, and the public documentation links resolve locally. A credential-pattern scan found no matches. Local logs and test scratch files are excluded from Git.

The desktop window and the plan/start controller calls were checked separately. No real-data building/refinement task was started. External scientific-engine tests include mock processes. Line/branch coverage was not measured. GitHub Actions results are separate from this local record.
