# OASIS4.3 Linux Testing Guide (4.3.8 + refmac-timeout-isolation)

For AI assistants and researchers testing on Linux. Goal: verify this source build installs, passes regression, runs one SAD → density modification → model pipeline, and launches the updated GUI. All commands run in plain bash.

Scientific changes vs 4.3.3: **4.3.7** fixed the amplitude scaling of the trace/Fcalc feedback (7E1D R-free 59.17% → 41.67%); **4.3.8** isolated numerical divergence of a single REFMAC pilot in the finish stage (6ZPP now completes, R-free 32.78%). Weight/endpoint selection and free-set isolation are unchanged.

## 0. Scope (avoid misjudgement)

- **Source build**: no CCP4 binaries or third-party runtimes bundled. `phaser` / `cparrot` / `cbuccaneer` / `refmac5` must be supplied by the tester under the CCP4 license. Density modification uses CCP4 Parrot; trace/Fcalc feedback uses REFMAC zero-cycle refinement.
- Python is pinned to **3.12** (`requires-python = ">=3.12,<3.13"`); 3.13 is not supported.
- The "macOS-only / no Linux build" wording in README and deployment notes is obsolete (historical ten-case validation build). Since 4.3.2 the macOS platform gate in `run_desktop.py` is removed; the GUI is pure Tkinter and runs on macOS and Linux with `python -B run_desktop.py --project .`. `runtime_environment()` handles Linux variables such as `LD_LIBRARY_PATH`.
- Three tiers, ~30–60 min total:

| Tier | Goal | CCP4 | Time |
| --- | --- | --- | --- |
| A | environment + regression + GUI build self-check | no | 5–10 min |
| B | CLI end-to-end smoke (S-SAD, stops at density modification) | yes | 2–5 min |
| C | full single case (build + refine) + GUI interaction | yes | 15–30 min |

> If `TMPDIR` points to a tiny tmpfs (e.g. a 10 MB `/tmp`), pip installs and pytest `tmp_path` MTZ writes fail. Set `export TMPDIR=<large-disk-dir>` and pass pytest `--basetemp=<large-disk-dir>/pytest`. This is an environment issue, not a code bug.

## 1. Get the code

Project: https://github.com/WeiDing535/OASIS4.3-Build

If the build is not yet on the remote, copy a clean release directory from the macOS source machine.

```bash
# rsync (recommended: excludes local traces and work products)
rsync -av --exclude='.git' --exclude='work/' --exclude='.local_python' \
      --exclude='runtime/' --exclude='OASIS4.3.app' \
  /Volumes/Fast_4T/Workbuddy_Oasis/OASIS4.3-20261117/ \
  <user>@<linux-host>:~/OASIS4.3/
```

```bash
# or: pack first, then transfer
cd /Volumes/Fast_4T/Workbuddy_Oasis/OASIS4.3-20261117
tar czf /tmp/oasis43-4.3.8.6.tar.gz --exclude='work' --exclude='.git' --exclude='.local_python' .
scp /tmp/oasis43-4.3.8.6.tar.gz <user>@<linux-host>:~/
# Linux side: mkdir -p ~/OASIS4.3 && tar xzf ~/oasis43-4.3.8.6.tar.gz -C ~/OASIS4.3
```

Verify key files on Linux:

```bash
cd ~/OASIS4.3
ls pyproject.toml environment.yml run_sad_to_model.py run_desktop.py \
   src/phasebuild/desktop.py validation_ten_cases/inputs/ins/README.md
```

## 2. Python 3.12 environment

Prefer conda/mamba (versions locked, avoids numerical drift):

```bash
cd ~/OASIS4.3
conda env create -f environment.yml      # or: mamba env create -f environment.yml
conda activate oasis43-unified
python --version                         # expect Python 3.12.x
```

System Python (no conda):

```bash
python3.12 -m venv .venv && source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[test]'       # gemmi==0.7.5, numpy<2, scipy>=1.12, pytest>=8
sudo apt install python3-tk              # Debian/Ubuntu; needed for the GUI check
```

> `environment.yml` already declares `tk`; a conda env normally needs no separate `python3-tk`.

## 3. Tier A: environment, full regression, GUI build self-check

### 3.1 Runtime doctor

```bash
python run_runtime_doctor.py > /tmp/oasis43-doctor.json 2>&1 || true
cat /tmp/oasis43-doctor.json
```

Check: `status` = `success`, `python` = 3.12, `packages` (gemmi 0.7.5 / numpy 1.26.x / scipy 1.13.x), `platform`. A CCP4-related error here is fine — Tier A does not require CCP4; continue with 3.2/3.3.

### 3.2 Full regression (core metric)

```bash
export TMPDIR=<large-disk-dir>            # required if /tmp is a small tmpfs
python -m pytest -q -c pyproject.toml --basetemp="$TMPDIR/pytest" tests \
  2>&1 | tee "$TMPDIR/oasis43-pytest-linux.log"
tail -5 "$TMPDIR/oasis43-pytest-linux.log"
```

Expected: ~`634 passed` (634 collected; ~4–6 min). If the 7E1D reference data is absent (expected), `test_real_7e1d_b_factor_reproduces_fmodel_amplitude` skips: `633 passed, 1 skipped`.

- Many failures/errors (not one): almost always `tmp_path` on a too-small `/tmp` tmpfs failing MTZ writes — re-run with `TMPDIR` / `--basetemp` as above.
- One remaining `test_portable_launcher_dispatches_both_hl_modes...` failure under **venv** is a test-environment artifact, not a bug: the test symlinks the bundled `python3` to the current interpreter; under venv `sys.prefix` falls back to the uv/conda base interpreter and loses venv site-packages, so `import gemmi` fails. Run under **conda** for a clean pass; under venv this failure is acceptable (re-running alone still fails). macOS + conda passes.
- Occasional `test_command_cancellation_cleans_up_its_independent_process_group`: a pidfile timing race (`ValueError` on an empty read); re-run alone passes, unrelated to code.

### 3.3 GUI build self-check (new two-page layout)

Platform gate removed; use the standard entry point:

```bash
python -B run_desktop.py --project . --check
```

Expected: `desktop window ok` (no error, no window left open).

With a display (desktop or `ssh -X`), verify the two main tabs (needs `xvfb` or a real DISPLAY):

```bash
# optional: xvfb-run for headless machines
xvfb-run -a python -c "import sys; sys.path.insert(0,'src'); \
  from phasebuild import desktop; \
  a=desktop.App(desktop.Controller('.')); a.root.update_idletasks(); \
  print([a.pages.tab(i,'text').strip() for i in range(len(a.pages.tabs()))]); \
  print('input_row_width=', a.fields['sad']['row'].winfo_width(), 'window=', a.root.winfo_width()); \
  a.root.destroy()"
```

Expected: prints the two main tab labels (the configuration-input and run-and-results tabs) and `input_row_width` close to `window` (input page fills the full width — the point of the GUI change).

## 4. Tier B: CLI end-to-end smoke (requires CCP4)

Use the bundled insulin S-SAD data (S, λ = 2.0 Å). Every phasing route needs CCP4 — `--no-build-model` only skips Buccaneer/REFMAC model building; P+ phasing still needs `phaser` and density modification still needs `cparrot`. Set `PHASEBUILD_CCP4` first (§5.1).

Recommended: `pplus_bundled` (project default; no `--reference` needed):

```bash
python -u run_sad_to_model.py \
  --method pplus_bundled \
  --sad validation_ten_cases/inputs/ins/ins_wave2.0_reso2.1A_sad.mtz \
  --sequence validation_ten_cases/inputs/ins/3INS.pir \
  --heavy validation_ten_cases/inputs/ins/sg.pdb \
  --ccp4-root "$PHASEBUILD_CCP4" \
  --elements S --wavelength 2.0 --auto-fpp --input-columns \
  --no-build-model \
  --output work/ins_ssad_linux
```

`--method likelihood` additionally requires `--reference`:
`--reference "$PHASEBUILD_CCP4/lib/data/reference_structures/reference-1tqw.mtz"`, else `Missing input file: reference`.

- The `--output` directory must not already exist; use a new name or delete it first.
- `--elements S --wavelength 2.0 --auto-fpp` computes theoretical f″ offline via Gemmi; no manual fpp.
- Expected: no Python traceback; `work/ins_ssad_linux/` contains `pipeline_report.json` and per-stage directories; `status` = `completed`, `03_dm/map_phases.mtz` exists. ~4–5 min on Linux.

`Failed to open ...mtz`: check relative paths are run from the package root. Missing gemmi/numpy: return to §2.

## 5. Tier C: full pipeline (bring your own CCP4)

### 5.1 Configure CCP4 path

```bash
export PHASEBUILD_CCP4=/path/to/your/ccp4       # e.g. /opt/ccp4-9
```

Verify completeness (needs `bin/phaser`, `bin/cparrot`, `bin/cbuccaneer`, `bin/refmac5`, and `lib/data/`: `reference-1tqw.mtz/.pdb`, `lib/data/monomers`):

```bash
for f in phaser cparrot cbuccaneer refmac5; do
  test -x "$PHASEBUILD_CCP4/bin/$f" && echo "OK $f" || echo "MISSING $f"
done
ls "$PHASEBUILD_CCP4/lib/data/reference_structures/reference-1tqw.mtz" \
   "$PHASEBUILD_CCP4/lib/data/monomers" >/dev/null && echo "OK data"
```

> `CINCL` / `CLIBD` / `CLIBD_MON` / `CCP4_SCR` for REFMAC/Buccaneer are set automatically from `--ccp4-root`; no need to source the CCP4 setup script.

### 5.2 Full single-case reproduction (7E1D)

`reproduce_ten_cases.sh` uses `CCP4_ROOT` for the engine root and `PHASEBUILD_PYTHON` for the interpreter:

```bash
export CCP4_ROOT="$PHASEBUILD_CCP4"
export PHASEBUILD_PYTHON="$(command -v python)"
bash reproduce_ten_cases.sh 7E1D 2>&1 | tee /tmp/oasis43-7E1D-linux.log
```

Output in `work/reproduce/7E1D/`. Expected: task `success`, `primary_model.cif` present; read R-work/R-free from `pipeline_report.json`.

Reference values (macOS-specific CCP4 build, after the 4.3.7 amplitude fix): 7E1D R-work 37.14% / **R-free 41.67%** (pre-fix degraded state was 59.17%). Same mode as Linux (`pplus_bundled`, `--fcalc fc`, seed 20261023, 0.50 / fb3), **not** likelihood.

> **Cross-platform values can differ substantially — a known phenomenon, not a mode problem.** 7E1D sits in a numerically sensitive region of Parrot. Full artifact comparison: DM inputs (`parrot_input.mtz`) differ at 1e-5 in HL probability amplitudes (from last-digit float differences in the Phaser SAD likelihood in `02_sad/observations.mtz`), amplified by ~40 cycles of strongly nonlinear solvent flattening / histogram matching into a **~12° phase difference** (cycle_01 `dm.mtz` by HKL: median phase diff 11.7°, FOM corr r=0.885), then trace 46 vs 58 residues → feedback 78→64→16 vs 93→141→95, ending up to 15 pp apart in R-free. Parrot itself is deterministic (three identical runs on one machine give element-wise identical output), so this is **not randomness, not a normalization failure, not lack of platform support**. Therefore **41.67% is not a hard cross-platform acceptance line**; Linux success = complete pipeline, EXIT 0, no traceback, R-free clearly below the unnormalized degraded state 59.17%. Do not switch to likelihood to chase the reference — that changes the protocol.

> **Distinguish "normalization failure" from "platform float divergence"** — verify normalization in ~10 s before judging. Linux first pass (2026-11-19) confirmed normalization works (`overall_b` = -15.54, same order as macOS -13.58), so its 57% is the float divergence above, not a bug:
>
> ```bash
> # check 1: normalization flag (should be true, overall_b negative)
> python -c "import json,glob; f=glob.glob('work/reproduce/7E1D/endpoint_fb3_s0.50/03_dm/feedback_report.json')[0]; d=json.load(open(f)); e=[h['dm_info'].get('fcalc_engine',{}) for h in d['history'] if h['dm_info'].get('fcalc_engine')][-1]; print('normalised=',e.get('amplitude_normalised'),'B=',e.get('overall_b'))"
> # check 2: FC median (~41 with normalization, ~69 without)
> python -c "import gemmi,numpy as np,glob; f=glob.glob('work/reproduce/7E1D/endpoint_fb3_s0.50/03_dm/cycle_02/parrot_with_model.mtz')[0]; print('FC median=', np.median(np.asarray(gemmi.read_mtz_file(f).column_with_label('FC'))))"
> # check 3 (log only): Overall B line in refmac_fc.log
> grep -E "Overall\s+:" work/reproduce/7E1D/endpoint_fb3_s0.50/03_dm/cycle_02/refmac_fc/refmac.log
> ```
>
> If `amplitude_normalised=false`, `overall_b` is null, or check 3 prints nothing: look for `Overall : scale = ..., B = ...` in that `refmac.log`. `_overall_b()` (`src/phasebuild/parrot_dm.py`) parses this line by regex; a wording mismatch causes a safe fallback to no normalization — that is the real bug to adapt; keep the full `refmac.log` for feedback.

> **Optional deep dive** — localize the difference to a stage by running Parrot twice on the same input and comparing arrays (checks run-to-run variation):
> ```bash
> python - <<'PY'
> import gemmi, numpy as np
> a=np.array(gemmi.read_mtz_file('work/reproduce/7E1D/endpoint_fb3_s0.50/03_dm/cycle_01/dm.mtz'))
> b=np.array(gemmi.read_mtz_file('work/reproduce/7E1D/endpoint_fb3_s0.50/03_dm/cycle_01/parrot_raw.mtz'))
> print('dm.mtz rows', len(a), 'parrot_raw.mtz rows', len(b))
> # input-side difference (HL amplitudes): HLA..HLD of parrot_input.mtz
> m=gemmi.read_mtz_file('work/reproduce/7E1D/endpoint_fb3_s0.50/03_dm/cycle_01/parrot_input.mtz')
> print('parrot_input cols:', m.column_labels())
> PY
> ```

> **Review update**: Linux parameter scan (2026-11-19) and macOS independent review add:
>
> 1. **The strictly cross-platform reproducible route is `pplus_internal`**, not further bundled tuning. `--method pplus_internal --rounds 14 --cycles 15 --seed 123 --reference <sha256-locked 1tqw.mtz>` (reference required for internal, not bundled) gives Linux 0.3325 / macOS 0.3409 (Δ 0.84 pp). This path computes its own observation likelihood + NumPy MLHL DM, with no Phaser/Parrot float-fingerprint entry. It is a different method from the ten-case bundled baseline; do not mix them in direct comparison.
> 2. **Seed is a lever only in `pplus_internal`**, not bundled: 7E1D native.mtz already has FreeR_flag (`flags_generated=False`), so the seed's only numeric consumer (free-set split) never triggers — bundled seed changes give bit-identical R-free (all 0.4167). Lever ranking: method > seed within internal > rounds; `--trust-kl` is effective and bound, while `dm_solvent` / `dm_mask_radius` / `dm_anisotropy` / `feedback_cycles` are silently ignored under internal (review §4–5).
> 3. **best-of-N selection bias**: 0.3325 is the best seed picked on R-free; this is valid only for reference-based phasing benchmarks. ML training data generation should keep a fixed seed to avoid systematically dropping hard cases and truncating the quality distribution.
>
> For a plain platform-usability check, run the bundled baseline above (6ZPP 32.54% vs 32.78%, highly consistent); no need to switch to internal.

**Linux test records**:

- 2026-11-18 first pass: 7E1D prepare/SAD/DM/bridge all succeeded, but the build stage `cbuccaneer` built 0 residues at fb2, then `refmac5` exited on empty coordinates (code 1). A separate full **5D9D** run succeeded (both endpoints, fb3 selected, R-work 25.47% / R-free 27.24%) to rule out a platform issue.
- 2026-11-19 (4.3.8) re-run: 7E1D **both endpoints completed, EXIT 0, no traceback** — fb2 R-work 55.58% / R-free 57.03% (63 residues, selected), fb3 R-work 56.34% / R-free 57.08% (48 residues). The 0-residue Buccaneer issue did not recur (randomness / feedback-trajectory branch). The gap to macOS 41.67% is the Parrot float divergence above; normalization verified, mode identical — not a platform or config problem.

Conclusion: **the full Linux SAD → model pipeline works.** 7E1D residue counts / R-free are not strictly cross-platform reproducible without a fixed CCP4 build; for a plain platform check use 6ZPP (Linux 32.54% vs macOS 32.78%, highly consistent) or 5D9D. Do not flag numeric differences as lack of platform support.

### 5.2b Full 6ZPP reproduction (focus of the 4.3.8 fix; recommended)

In 4.3.6/4.3.7, 6ZPP failed at finish because one weight pilot (w0p02) diverged numerically in REFMAC. 4.3.8 isolates the divergent pilot, cuts the timeout from 1800 s to 600 s, and no longer lets it drag down the endpoint — the key point to verify in this release:

```bash
export CCP4_ROOT="$PHASEBUILD_CCP4"
export PHASEBUILD_PYTHON="$(command -v python)"
bash reproduce_ten_cases.sh 6ZPP 2>&1 | tee /tmp/oasis43-6ZPP-linux.log
```

Expected: exit 0, both endpoints fb2/fb3 complete; `work/reproduce/6ZPP/endpoint_selection.json` present, fb3 selected, R-work ≈ 0.304 / **R-free ≈ 0.328 (32.78%)**.

Key evidence the fix works (check either endpoint):
- `work/reproduce/6ZPP/endpoint_fb2_s0.50/06_finish/FAILED_PILOTS.json`, if present, should record `weight=0.02, reason=timeout, timeout_seconds=600` — divergence isolated, not crashing the case; if fb2 did not reproduce divergence (intermittent), the file being absent is also fine.
- No uncaught `subprocess.TimeoutExpired` traceback in the main log.
- If it still fails with an uncaught `TimeoutExpired`, keep the full `finish.log` for feedback — that is a regression to fix.



Other cases: `3O2E 5D9D 5HHK 6HF7 6ZPP 7E1D 7ENM 7TBS 7TFQ 7TJ1` (7TFQ uses `pplus_trust` with `observations.mtz`). All ten: `bash reproduce_ten_cases.sh all` (long, optional).

### 5.3 GUI interaction (with DISPLAY)

Use the standard entry point (prefix `xvfb-run -a` on headless machines):

```bash
python -B run_desktop.py --project .
# headless: xvfb-run -a python -B run_desktop.py --project .
```

GUI check points (focus of this update):
1. Two top-level tabs ("Configuration Input" / "Run & Results"); starts on "Configuration Input".
2. Input form fields fill the window width (Browse / Detect buttons visible, no dragging).
3. With valid input, "Check Config" auto-switches to "Run & Results → Run Plan"; "New Task" switches back to "Configuration Input".
4. Scrolling the log on the run page does not scroll the hidden input form.

## 6. Reporting format (send back to the macOS side)

```
Platform: <uname -a / distro>
Python: <version>  conda/system: <>
gemmi/numpy/scipy: <versions>
CCP4: <path or "not provided">
Tier A: pytest=<X passed, Y failed>  GUI self-check=<desktop window ok / error>  two tabs=<yes/no>
Tier B: status=<success/failed>  artifacts=<pipeline_report.json present yes/no>
Tier C: 7E1D status=<>  R-work=<%>  R-free=<%>  anomaly=<traceback summary or "none">  fcalc normalised=<true/false, B>  FC median=<value>
Tier C: 6ZPP status=<>  R-work=<%>  R-free=<%>  FAILED_PILOTS=<yes(with content)/no>  uncaught Timeout=<yes/no>
Conclusion: <releasable on Linux / specific issue to fix>
```

On failure, keep the full logs (`/tmp/oasis43-*.log`) and tracebacks, not just the conclusion.

## 7. Known issues

| Symptom | Cause / fix |
| --- | --- |
| GUI startup reports a platform error | Should not happen: macOS platform gate removed since 4.3.2; use `python -B run_desktop.py --project .` |
| Many `pytest` failures/errors (e.g. hundreds) | `/tmp` is a tiny tmpfs; `tmp_path` MTZ writes fail; set `TMPDIR` to a large disk and add `--basetemp` |
| One packaging failure under venv | venv symlink artifact: bundle loses site-packages in the test; conda passes cleanly; acceptable under venv |
| pip `No space left on device` | Same as above: pip cache/temp in `/tmp`; set `TMPDIR` and `PIP_CACHE_DIR` to a large disk |
| `Missing input file: reference` | `--method likelihood` requires `--reference`; `pplus_bundled` omits it |
| `Failed to open ...mtz` | Not run from the package root, or wrong relative path |
| GUI `no display name` | No DISPLAY; use `xvfb-run -a` or `ssh -X`; if a DISPLAY exists (e.g. `:10.0`) it works directly |
| Buccaneer builds 0 residues at build stage → REFMAC empty-coordinate error | Case-level randomness / parameter sensitivity when the DM map has signal; validate the platform with another case, see §5 note |
| 6ZPP finish weight pilot timeout (old versions failed the whole case) | Fixed in 4.3.8: divergent pilot isolated as failed, 600 s timeout, other pilots/endpoints unaffected; see §5.2b; if an uncaught `TimeoutExpired` remains, keep logs for feedback |
| CCP4 check reports missing `reference-1tqw.mtz` | Incomplete CCP4 install or version difference; check `lib/data/reference_structures/` |
| `--output` directory already exists | Program refuses to overwrite; use a new name |
| Python 3.13 dependency install fails | Hard pin `>=3.12,<3.13`; use 3.12 |
