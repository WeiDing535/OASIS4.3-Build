# Downstream test-only portability adjustment

The repeated-phase assertion in `tests/test_science.py::test_native_modes_return_finite_phases_and_reproducible_results` compares circular phase differences and FOM with absolute tolerance `4 * sqrt(float64 epsilon)` (about 5.96e-8 radians for phase, relative tolerance zero), instead of bitwise phase-array equality. The existing phase finiteness, FOM bounds and error-weight checks remain.

In the unified conda macOS environment, the DMR case initially differed at 1 of 364 phase entries by 4.44e-16 radians. A later P+ SAD repeat differed by up to 1.99e-8 radians. These old doublet tests include arccos arguments at +/-1; its endpoint conditioning amplifies O(epsilon) argument rounding to O(sqrt(epsilon)) phase rounding. Both were observed with numerical-library thread counts fixed to one. The tolerance follows that conditioning scale and adds a FOM-repeat check; it is not a structure-solving success threshold.

No `oasis_py` numerical source is modified. Source hashes remain those of the upstream commit recorded in UPSTREAM.json. This is an explicit downstream test portability adjustment, not a claim of bitwise reproducibility across math-library builds.
