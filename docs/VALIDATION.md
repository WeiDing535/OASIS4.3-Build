# Publication validation — 2026-10-22

The public source export passed 186 tests (0 failures, 0 errors, 0 skipped) in 20.67 seconds on the existing macOS Intel environment: Python 3.9.18, Gemmi 0.7.5, NumPy 1.19.5, SciPy 1.7.3, pytest 6.2.5.

Command: `python -B -m pytest -q -p no:cacheprovider --basetemp validation/tmp --junitxml validation/tests.xml tests`.

All 47 Python modules and GUI assets match the tested local trial copy byte-for-byte. Source files parse successfully, the launcher passes `bash -n`, source packaging passes its CRC check, and the public documentation links resolve locally. A credential-pattern scan found no matches. Local logs and test scratch files are excluded from Git.

The trial GUI and token-protected localhost defaults API were checked separately. No real-data building/refinement task was started. External scientific-engine tests include mock processes. Line/branch coverage was not measured. GitHub Actions results are separate from this local record.
