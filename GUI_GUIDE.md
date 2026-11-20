# GUI User Guide

The integrated build adds "SAD phasing to model". Requires a Python 3.12 conda environment. Inputs: SAD anomalous data, the original native free set, an existing heavy-atom model, sequence, and an unrelated density-modification reference. Raw experimental HL and the density-modified map are stored separately; resume via the CLI `--resume`. See [INTEGRATION.md](docs/INTEGRATION.md).

This is a source build: no bundled legacy portable package and no third-party CCP4 runtime. Bring your own licensed CCP4. Under "Show Advanced Parameters", set the `CCP4 Installation Root` field or click "Detect" to auto-find a standard install (CLI equivalent: `--ccp4-root`; legacy `PHASEBUILD_RUNTIME` / `--runtime` are kept for compatibility only). From the project root, double-click `OASIS4.3.command` (`PhaseBuild.command` is deprecated — do not use), or run `python -B run_desktop.py --project .` with an existing scientific Python.

Platform: the desktop shell is Tkinter; it runs on macOS and Linux (Apple Silicon via Rosetta; Linux needs a display — add `xvfb-run -a` on headless). The CCP4 engine is user-supplied and not subject to this platform limit. See `LINUX_TESTING.md` for Linux tests.

Native desktop window (Tk, shipped with Python): no web server, no listening port, no browser. Data and computation stay local. Closing the window exits; running tasks stop with it; result directories remain.

## Quick start

1. Click "New Task" (or Cmd+N / Ctrl+N) to clear the previous form and suggest a new output directory; optionally set a task name at the top (e.g. "7E1D rerun") — it appears in the task record.
2. On the left form choose a phasing method and fill in the experimental inputs: SAD data, heavy-atom model, sequence, anomalous element, wavelength.
3. Make sure the output directory does not already exist.
4. If needed, check "Show Advanced Parameters" to review the runtime environment, original free-set value and selected MTZ columns; click "Check Config".
5. Review the steps on the "Run Plan" page, then click "Start Run".
6. Follow progress on "Run Log"; view per-round R-work/R-free curves, stages, engines and refinement cycles on "Refinement Results"; view the final model on "Structure 3D". Failed or cancelled outputs are kept and are not a complete endpoint.

## Form sections

Only per-experiment parameters are shown:

- Phasing method: P+ method, parameter scheme, density-modification engine.
- Input files: SAD data, heavy-atom model, sequence, element, wavelength, output directory.
- Sample and density modification: copy number, solvent content.
- Model building and refinement: whether to build a model, auto hand determination, refinement engine, refinement rounds and cycles per round.

The two engine selectors are explicit choices, not automatic detection. The density-modification selector offers CCP4 Parrot (default) and the in-house `handdm`; the refinement selector offers CCP4 REFMAC (default) and the in-house `minirefmac`. Selecting an in-house engine needs its deployment to be reachable (the in-house programs live outside the package); if it cannot be found, the run stops at configuration validation instead of falling back to the external engine. The in-house engines differ from the external ones in documented ways — for example `handdm` cannot run the trace-feedback loop, so that request is refused rather than silently downgraded.

Remaining parameters (native data, observational priors, combined sequence, MTZ column names, free-set value, loop/sampling limits, P+ optimization constants, reference structure, runtime environment, feedback endpoint combination) are hidden by default. Check "Show Advanced Parameters" at the bottom to expand, grouped under "Input Details / Loops & Sampling / P+ Optimization / Run Resources & Endpoint Combination". They are calibration constants or run resources; defaults are normally fine. Hiding affects display only, not meaning or defaults. Detect external engine paths with "Detect"; if detection fails, enter the `--ccp4-root` directory manually.

## Download results

- "Export Results…" (bottom button, or top-right of "Refinement Results"): copies the whole output directory of the current task to `oasis_export_<time>`, including model, phase MTZ, per-round reports, `endpoint_selection.json` and logs; the status bar then shows file count and size.
- "Save Final Model…" (top-right of "Refinement Results"): saves the final model alone as mmCIF/PDB, for other programs or direct download.

## Refinement curves

The top of "Refinement Results" plots R-work / R-free vs round: solid R-work, dashed R-free, y-axis percent, x-axis round (last point = fixed post-refinement); the legend shows the latest values. The curve and table below come from the same per-round report; with a single round the page shows a hint instead of a curve.

## Structure 3D

"Structure 3D" draws the final model in-window, no browser or external viewer: by default the Cα trace coloured by chain; drag to rotate, scroll to zoom, "Reset View" to restore orientation. "Show All Atoms" loads all atoms from the model file. The page title shows file name, chain count and residue count; with no result it shows a placeholder.

## Task record and parallel runs

The first right-side page, "Task Record", lists every task run on this machine: start time, name, status, method, final R-work and output directory. Selecting a row loads that task into the "Run Log / Refinement Results / Structure 3D / Run Plan" pages; the window stays on "Task Record" (no auto-jump, buttons stay enabled); double-click opens the result directory.

Two buttons at the top-right of the record page:

- "Reuse Parameters": loads the selected task's full parameter set (including advanced) back into the left form with a new suggested output directory. Jobs with identical parameters but different input data need no re-entry — change the input files and click "Start Run". Imported existing results have no saved parameters and cannot be reused; the UI explains why.
- "Delete Record": removes the selected row(s) from the task record (with confirmation). Only entries under `.gui-oasis/jobs` are deleted; result directories and files are kept; a running task must be stopped first. Records live in `.gui-oasis/jobs` under the project and survive restart, so you need not track runs yourself. "Import Existing Result…" only registers an existing output directory in the record for review: status is fixed to "Imported (view only)", no re-run and no run log, but curves, tables, export and save-model work as usual.

"Start Run" no longer waits for the previous task: each click creates a task that runs in parallel in its own directory, log and subprocess group. After start, the form switches to a new suggested output directory so you can queue the next job immediately; up to 4 tasks at once, with a wait prompt beyond that. Each task still needs its own fresh output directory; existing results are never reused or overwritten. The window top-right shows "N tasks running"; closing the window confirms stopping all of them at once.

## Native mode

Model building from a phase map needs a phase MTZ, sequence, and a training reference structure from the user-supplied CCP4 runtime directory. Refinement additionally needs a native MTZ and REFMAC. Extending an existing model additionally needs the model and its sequence-assignment report. Column choices apply only in native flows that support them.

The native GUI can select a stage by Rfree within budget; this is an exploratory choice, not an untouched final test.

## Two HL integration modes

"No Initial Model + Raw HL Refinement" needs no seed; "Buccaneer Integration + Raw HL Refinement" requires an existing seed model. Both build with external Buccaneer and refine with REFMAC, ending at the fixed final planned round.

The current integration interface requires unique column names: map F/SIGF/FWT/PHWT/FOM, native F/SIGF/FreeR_flag, raw experimental priors HLA/HLB/HLC/HLD. Cell, space group and indexing must match; raw priors are matched by exact HKL. Missing priors stay missing.

The free-set value in advanced settings must match the input (default 0); original flags are not re-encoded. The same raw experimental prior is reused each round, never replaced by a model posterior. Running needs CCP4 `cbuccaneer`, `refmac5` and training reference resources.

## Result interpretation

Independent R factors are evaluated at fixed coordinates and are not directly interchangeable with REFMAC values. No fabricated R value is shown when refinement/evaluation is not run. After a service restart, unfinished tasks are marked interrupted; there is currently no general GUI resume button. The new no-seed entry has not yet been validated on new real data.
