# Current validation pointers

This release is the ten-case validation build `4.3.0.dev17-ten-case` (commit `3fa5b86`, dev17 portable runtime).
The full test suite measured **776 passed / 0 failed / 1 deselected** (2026-11-16, Python 3.12.15 / gemmi 0.7.5).
Ten-case solving results and reproduce commands are in the root-level deployment notes and `validation_ten_cases/verify_3fa5b86_20261116/`;
stage-artifact reading guidance is in [`READING_OUTPUT.md`](READING_OUTPUT.md).

Historical research records (some reports are not distributed with this package; only still-available documents are listed):
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
