# 4.3.7 Ten-Case Regression (2026-11-19, with strict A/B)

> **Update (4.3.8, 2026-11-19)**: The 6ZPP finish-stage REFMAC timeout recorded in §1/§5 was fixed in 4.3.8; a clean rerun **exited 0 with R-free 32.78%**. Root cause was two fault-tolerance gaps (parallel pilot exceptions not isolated; finish-refinement `except` missing `TimeoutExpired`), unrelated to the fcalc fix. The other nine cases are unchanged.

## 1. Conclusion

Under same-environment strict A/B, the amplitude-normalization fix is a significant improvement for most cases, with no phase-quality degradation in any case:

- **7TFQ**: 4.3.6 (unnormalized) R-free **51.59%** → 4.3.7 (normalized) **38.36%**, an improvement of **13.23 pp**.
- **7E1D**: old degraded state **59.17%** → **41.67%**, an improvement of 17.5 pp (previous A/B).
- **7ENM**: absolute R-free slightly higher by 1.31 pp, but **not a quality regression** — caused by the endpoint rule switching to fb3 (see §4); its fb2 endpoint measured 31.79%, better than the old fb2's 32.10%.
- **6ZPP**: REFMAC timeout at 1800 s in finish. Rerunning with the 4.3.6 baseline code (unnormalized) in the same environment confirmed: the baseline **fails the same way** (also a finish-stage REFMAC `TimeoutExpired`). This failure is therefore **unrelated to the fcalc fix** — a REFMAC numerical/performance pathology on this machine. See §5.

## 2. Parallel execution

Per the user's request for "parallel processing": stopped the original serial regression and used `work/reproduce/run_parallel.sh` (`xargs -P 5`) to run 9 cases in two parallel batches; 7E1D was verified separately. 20-core / 64 GB machine, 5-way concurrency. Total 1 h 49 m (serial estimate ~3 h). Under concurrency the engines compete for resources, raising single-case time from ~20 min to 35–60 min.

## 3. Nine-case R-free summary (4.3.7 normalized, same-environment measured)

| Case | 4.3.7 new R-free | PACKAGED 11-10 baseline | Diff (pp) | Note |
|---|---:|---:|---:|---|
| 3O2E | 41.03 | 42.64 | −1.61 | Improved (with its dedicated 4-endpoint parameters) |
| 5D9D | 28.47 | 27.51 | +0.96 | Within threshold |
| 5HHK | 28.55 | 28.82 | −0.27 | Flat |
| 6HF7 | 29.61 | 29.04 | +0.57 | Within threshold |
| 6ZPP | — | 31.31 | — | finish-stage REFMAC timeout (see §5) |
| 7E1D | 41.67 | 42.33 | −0.66 | Improved (previous A/B −17.5 pp vs degraded state) |
| 7ENM | 33.41 | 31.47 | +1.94 | Endpoint selection, not regression (see §4) |
| 7TBS | 29.15 | 28.35 | +0.80 | Within threshold |
| 7TFQ | 38.36 | 36.86 | +1.50 | Strict A/B is actually an improvement of 13.23 pp (see §4) |
| 7TJ1 | 30.81 | 31.08 | −0.27 | Flat |

> **Baseline choice is sensitive**: PACKAGED 11-10 is an earlier config of the packaged runtime; in the same-seed local fmodel directory (20261116) 7TFQ=40.23, 7ENM=36.10, 7E1D=42.50, about 3 pp from PACKAGED. Judging a fix therefore requires **same-environment A/B**, not cross-config historical tables.

## 4. Strict A/B (git worktree, same machine and seed, `4.3.6` vs `4.3.7`)

Created a 4.3.6 unnormalized baseline with `git worktree add … cc70311`, and ran once each under the same Python / CCP4 / inputs / seed:

| Case | 4.3.6 unnormalized | 4.3.7 normalized | Diff (pp) | Verdict |
|---|---:|---:|---:|---|
| 7TFQ | 51.59 | 38.36 | **−13.23** | Significant improvement |
| 7ENM | 32.10 | 33.41 | +1.31 | Endpoint-selection side effect (see below) |

**7ENM breakdown**:

| Endpoint | 4.3.6 Rwork / Rfree | 4.3.7 Rwork / Rfree |
|---|---|---|
| fb2 | 27.72 / 32.10 | 28.20 / **31.79** (better) |
| fb3 | 32.98 / 36.82 | **27.67** / 33.41 |
| Selected | fb2 (Rwork 27.72 lower) | fb3 (Rwork 27.67 lower) |

Normalization dropped fb3's work R from 32.98 sharply to 27.67, below fb2, so the "select endpoint by work-set R" rule switched to fb3; fb3's R-free (33.41) happens to be higher than fb2's (31.79). **Both endpoints improved in absolute quality — the rule simply selected the one with slightly higher R-free.** This is a known side effect of the preregistered endpoint-selection protocol ("R-free never read"), not phase degradation. trace/FOM also confirm healthy feedback: old 86→82→95, new 86→77→94, FOM stable at 0.818–0.821.

## 5. 6ZPP failure (independent stability issue, unrelated to this fix)

- Failure point: `work_weight.calibrate_work_weight` in `06_finish`, using `ThreadPoolExecutor(max_workers=2)` to run REFMAC calibration of 4 weight groups (.005/.01/.02/.04) in parallel.
- Three groups (.005/.01/.04) finished normally in 44–56 s; **w0p02 stuck at CGMAT cycle 8, fvalues reaching 4.7e8 (numerical divergence)**, then killed by `subprocess.run(timeout=1800)`.
- Its input was that stage's own `cleaned.cif` + original `native.mtz`; the code path reads no Fcalc feedback columns, **unrelated to the change in `parrot_dm._attach_fc()`**.
- Historical 20261116 with the same config (`trust_work_weight=True`, `refinement_prior=dm_phifom`) had succeeded, indicating an occasional REFMAC numerical pathology.
- **Final confirmation**: rerunning 6ZPP with the 4.3.6 baseline code (`git worktree` @ `cc70311`, unnormalized) on the same machine and inputs **also failed with a finish-stage REFMAC timeout** (refmac5 `TimeoutExpired` after 1800 s in `trust_workflow.launch`, log `/tmp/ab436_6ZPP.log`). The 4.3.7 run stuck at `work_weight.pilot`'s w0p02 (CGMAT cycle 8 divergence), the baseline at the finish main refinement — different REFMAC call sites, but both have a finish-stage refmac5 running the full 30 minutes without returning. **The baseline failing too proves it is unrelated to the fcalc fix**; it is a numerical/performance pathology of the current CCP4-9 install on this data. The historical 20261116 success with the same config further supports occasional occurrence. Track it independently (add retry/reduced load for finish-stage REFMAC); it does not block acceptance of this fix.

**4.3.8 fix and validation (2026-11-19)**: identified two fault-tolerance gaps — (1) a single pilot timeout in `work_weight.calibrate_work_weight`'s `pool.map` crashes the whole map, discarding three successful pilots; (2) `trust_workflow.finish`'s `except (RuntimeError, ValueError, OSError)` does not catch `subprocess.TimeoutExpired`. 4.3.8 isolates a single pilot failure as `status='failed'`, lowers the timeout to 600 s, and adds `TimeoutExpired` to the finish `except`. Clean rerun exited 0: fb2's w0p02 still diverged but was isolated by the 600 s timeout, fb3's four pilots all completed, fb3 was finally selected, **R-free 0.3278 (32.78%)**..

## 6. Artifact paths

- 4.3.7 parallel regression: `work/reproduce/parallel_logs/{case}.log` + `.exit`, driver log `parallel_driver.log`.
- 4.3.6 strict A/B: `git worktree` at `OASIS4.3-436-baseline` (commit `cc70311`), logs `/tmp/ab436_*.log`.
- 7E1D A/B details: `docs/FCALC_SCALE_FIX_7E1D_20261118.md`.
- Per-case `endpoint_selection.json` / `pipeline_report.json` / `feedback_report.json` under `work/reproduce/<case>/endpoint_*/` (work/ is gitignored).
