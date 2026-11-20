# OASIS4.3-Build (PhaseBuild)

OASIS4.3 is a research pipeline for structural biology / crystallography: SAD phasing → density modification → automated model building, end to end. Version `4.3.8.6` (internal PEP440 version `4.3.8+refmac-timeout-isolation`).

> **Limitations (read before use)**
>
> - **R-free is not an end-to-end independent blind-test metric.** SAD and density modification exclude the original free set, but the downstream external map handoff and model refinement are not certified as work-set-only. R-free is an in-protocol reference, not an estimate of generalisation to unknown structures.
> - **"10/10" means**: all ten cases complete the full pipeline under a fixed protocol, with R-free within ±5 percentage points of an independent reference. Some absolute R-free values remain high (e.g. 3O2E 42.64%, 7E1D 42.33%). It does not mean "ten high-quality structures solved".
> - **Defaults ≠ validation protocol.** Validation tuned each case individually (e.g. 7TFQ uses `--method pplus_trust` with observational priors). GUI/CLI defaults differ; applying defaults to a new sample does not reproduce the ten-case performance.
> - **External engines are user-supplied and licensed by the user.** The source distribution contains no Buccaneer / REFMAC / Parrot / CCP4 components and no bundled 3.5 GB runtime. Users obtain CCP4 and comply with its licence, specifying the install root via `--ccp4-root` (or the GUI field). Density modification uses CCP4 Parrot; trace/Fcalc feedback uses REFMAC zero-cycle from the same CCP4 install. The historical validation runtime contained licensed third-party components and was restricted to internal or authorised use; it is not publicly redistributable.
> - This is a **research implementation requiring experimental validation**, not a calibrated commercial structure-determination tool.

Docs: [user guide](docs/OASIS4.3.md) · [GUI guide](GUI_GUIDE.md) · [Linux testing](LINUX_TESTING.md) · [validation](docs/VALIDATION.md)

## What is it

OASIS4.3 chains SAD phasing, density modification (DM) and automated model building into one end-to-end pipeline, with two fixed pairs:

| Phasing method | Density modification | Downstream |
| --- | --- | --- |
| Continuous SAD (`likelihood`) | MLHL-best | Seedless Buccaneer building + REFMAC observed-HL refinement |
| Direct-method P+ | Legacy DM | Same |

- **Input**: SAD MTZ, raw native observations with free set, heavy-atom coordinates, protein sequence, experimental wavelength (theoretical f″ computed from element and wavelength).
- **Output**: per-stage phases / maps / models, per-cycle R-work and R-free, `primary_model.cif`, `pipeline_report.json`, `endpoint_selection.json`.
- **Scope**: research implementation. The numerical core preserves the upstream OASIS Python implementation; although the full pipeline has been run on ten cases, each new sample requires experimental validation. No universal accuracy is claimed for arbitrary new samples.

## Install

### (a) Source repository (bring your own dependencies)

Requires Python 3.12 and a pinned scientific stack: gemmi 0.7.5, numpy 1.26.4, scipy 1.13.1.

```bash
conda env create -f environment.yml
conda activate oasis43-unified
```

The full pipeline also needs **user-supplied external engines**: CCP4 `phaser`, `cparrot`, `cbuccaneer`, `refmac5` and DM reference resources. For phasing + density modification only, add `--no-build-model`; no external building engine is then needed. The source distribution includes no CCP4 runtime.

### (b) Bring your own external engines (user-supplied, licensed CCP4)

The source distribution ships no 3.5 GB third-party runtime. The CCP4 engines needed for the full pipeline must be supplied and licensed by the user; set the CCP4 root with `--ccp4-root` (the GUI field offers auto-detect). Density modification uses CCP4 Parrot (`cparrot`); trace/Fcalc feedback uses REFMAC (`refmac5`) zero-cycle structure factors, both from the same CCP4 install. The historical dev17 portable runtime (`runtime/`, ~3.5 GB, with `cparrot`/`cbuccaneer`/`refmac5`/`phaser` and `density/`) was a restricted validation asset containing licensed third-party components; per upstream licences it is not publicly redistributable and is limited to internal or authorised channels. See [Linux testing](LINUX_TESTING.md).

### Launch the GUI

- Double-click `OASIS4.3.command` in the project root, or `OASIS4.3.app` on the desktop;
- or run `python -B run_desktop.py --project .` in the conda environment above.

The window is a native Tk desktop app: no web server, no network access, data stays local. Built on Tkinter; runs on macOS and Linux (Linux needs a display; headless: `xvfb-run -a python -B run_desktop.py --project .`; Apple Silicon needs Rosetta). CCP4 engines are user-supplied. Linux test steps: `LINUX_TESTING.md`.

## Quickstart

Insulin S-SAD example in the repository (element S, wavelength 2.0 Å; six sulfur sites from `sg.pdb`):

```bash
python -u run_sad_to_model.py \
  --method likelihood \
  --sad validation_ten_cases/inputs/ins/ins_wave2.0_reso2.1A_sad.mtz \
  --sequence validation_ten_cases/inputs/ins/3INS.pir \
  --heavy validation_ten_cases/inputs/ins/sg.pdb \
  --reference <your CCP4>/lib/data/reference_structures/reference-1tqw.mtz \
  --elements S --wavelength 2.0 --auto-fpp \
  --input-columns --no-build-model \
  --output work/ins_ssad_quickstart
```

This is **S-SAD**: `--elements S --wavelength 2.0 --auto-fpp` lets Gemmi's Cromer–Liberman algorithm compute theoretical f″ offline. `--no-build-model` stops after density modification, so no external building engine is needed. To continue into model building, drop the flag and set `--ccp4-root` to your CCP4 install (the GUI can auto-detect). The output directory must not already exist.

## Results / Reading output

The `--output` directory (`00_sad_to_model/` in the ten cases) is organised by stage; with multiple endpoints each endpoint gets a subdirectory:

| Directory | Contents |
| --- | --- |
| `01_inputs/` | SAD input verification, original free set, expanded modelling sequence |
| `02_sad/` | Phasing result `experimental.mtz` (PHIB/FOM/HL) and pure-observation likelihood `observations.mtz` (fixed refinement prior) |
| `03_dm/` | Density-modified phases / FOM (`map_phases.mtz`) and per-cycle mean-FOM trace |
| `04_bridge/` | Map handoff `building.mtz`: native amplitudes + `FWT = F × final DM FOM` + PHWT |
| `05_build/` | Per-cycle Buccaneer building, REFMAC refinement, refined maps and logs |
| `06_finish/` | Single-copy cleanup / work-set weight calibration / post-refinement; final model copy also here |

Key files:

- `primary_model.cif`: final model from the fixed last cycle (at the output root, with a copy in `06_finish/`). **Not** selected by lowest R-free.
- `pipeline_report.json`: auditable pipeline record. `status: "completed"` means all planned stages finished and file hashes verified; otherwise check `error`. With building off, `primary_model` is `null` and `dm_phases` is the final phase file — still a successful completion, not a building failure.
- R-work / R-free: per cycle in `pipeline_report.json` `primary_metric` (and per-endpoint records) and the GUI refinement curves, as fractions or percentages. R-work fits the working set; R-free is an in-protocol reference only. Neither is comparable across methods or software, nor an end-to-end blind score.
- `endpoint_selection.json`: multi-endpoint grid record. The `rule` states *lowest final work-set R across the predeclared endpoint grid; R-free never read* — endpoints are selected by R-work only, never by R-free. Lists every endpoint's `r_work`/`r_free`, the chosen endpoint and the selection rule.

Field details and CLI resume: [user guide](docs/OASIS4.3.md) and [GUI guide](GUI_GUIDE.md).

## Ten cases

Validated on ten cases under a fixed protocol: 3O2E, 5D9D, 5HHK, 6HF7, 6ZPP, 7E1D, 7ENM, 7TBS, 7TFQ, 7TJ1. See the [user guide](docs/OASIS4.3.md) and deployment notes for per-case settings and results.

## Cite

Method paper: Ding W, Zhang T, He Y, Wang J, Wu L, Han P, Zheng C, Gu Y, Zeng L, Hao Q, Fan H. *IPCAS: a direct-method-based pipeline from phasing to model building and refinement for macromolecular structure determination*. **Journal of Applied Crystallography** 2020;53:253–261. DOI: [10.1107/S1600576719015115](https://doi.org/10.1107/S1600576719015115).

Software: OASIS4.3-Build (PhaseBuild), version `4.3.8.6`, 2026-11-20. Full citation: [CITATION.cff](CITATION.cff). Maintainer: Ding Wei, Institute of Physics, Chinese Academy of Sciences. Licence LGPL-2.1-or-later ([LICENSE](LICENSE)); third-party packages, engines and reference resources follow their own licences — see [PROVENANCE.md](docs/PROVENANCE.md). Source: https://github.com/WeiDing535/OASIS4.3-Build
