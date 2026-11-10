# OASIS4.3-Build

NSF-D2 experimental standalone CLI: Rfree37.42% ->31.28%, still above30%. Options remain off by default. [Study report](docs/NSF_D2_IMPROVEMENT_20261026.md).

当前新增 dev5 自带引擎模式：无需另装套件，SEC17 完整隔离复验 Rfree 27.42%。详见 [独立版验证记录](docs/SELF_CONTAINED_PPLUS_20261026.md)。仅供已有授权组件的本机评估。

[中文](README.md) | English

Research software with two SAD methods, paired density modification, and optional external Buccaneer/REFMAC construction and refinement. The development repository retains the PhaseBuild 0.8.3 components and an unchanged upstream OASIS numerical-source snapshot.

## Integrated development workflow

The focused desktop window offers continuous SAD + MLHL-best or P+ + legacy DM. Model building is enabled by default; disabling it produces a successful phases-only endpoint without checking or invoking Buccaneer/REFMAC. Building starts without a protein model seed and uses a separate, fixed observation-likelihood HL prior for every refinement round. Continuous SAD's unchanged phase-result HL includes the final triplet/partial-structure field and remains the DM input. Since 4.3.0.dev2, `observations.mtz` excludes that field for refinement. P+ observation HL values are retained unchanged. See the [current usage guide](docs/OASIS4.3.md), [integration notes](docs/INTEGRATION.md), and [version policy](docs/VERSIONING.md).

## Linux-first complete bundle

The next complete runtime target is Linux x86_64, with Python and scientific engines included. The Linux launcher, ELF/dependency audit and tar.gz assembler are implemented; an actual Linux runtime bundle still requires native runtime inputs and clean-system validation. See [Linux preparation](docs/LINUX.md). The existing Intel Mac local-evaluation bundle is retained. No public complete binary release has been made.

## Run from source

Use an Python 3.12 environment with Gemmi, NumPy and SciPy. External construction/refinement additionally requires a user-provided CCP4 installation and reference resources. No Python or CCP4 runtime is bundled.

Double-click `OASIS4.3.command` for the native desktop window, or run `python -B run_desktop.py --project .` from the repository root in the Python 3.12 [conda environment](environment.yml). The window is drawn with Tk, which is part of the Python runtime: no browser, no local web server and no port forwarding. No dependencies are installed automatically.

Everything the scientific workflow needs is shipped or installed by the environment above. Optional helpers that are not part of the core — the legacy Qt front end (`oasis_py.gui`) and plotting — are installed by the user with `python -m pip install ".[gui]"` or conda; see the [usage guide](docs/OASIS4.3.md). Closing the window stops the running job. Always use a new output directory.

## Scope and validation

The focused product exposes the two SAD/DM routes and optional building only. Other native building, extension, tracing and experimental workflows remain in the development repository for research and regression testing.

The recovery baseline passed 385 cases; the first recovery changes passed 388, and version 4.3.0.dev2 passed 394, with one reference-data test deselected. Integration contracts include substituted external engines. The separate real-engine smoke results and coverage limits are recorded in the [development recovery report](development/RECOVERY_REPORT_20261023.md). The insulin diagnostic now produces a partial 46/51-residue continuous-SAD model with fixed final Rwork/Rfree 0.2769/0.3036 using independently refined heavy sites; P+ remains unsuccessful. This is not general performance validation or a complete structure. See the [insulin diagnosis](docs/INSULIN_DIAGNOSIS_20261023.md). A new complete focused binary bundle remains unvalidated.

Run tests using an existing pytest installation: `python -B -m pytest -q -p no:cacheprovider tests vendor/oasis43/tests`.

Maintained by DingWei (丁玮), Institute of Physics, Chinese Academy of Sciences. LGPL-2.1-or-later; see LICENSE. Third-party packages, engines and reference resources retain their own terms. See [provenance](docs/PROVENANCE.md).
