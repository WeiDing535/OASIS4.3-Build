# Linux x86_64 first

Linux x86_64 is the primary target for the next complete runtime package. Users of that package should only extract it and launch the application; they should not install CCP4, Python or scientific Python packages separately.

## Current evidence

The source regression suite already runs successfully on Ubuntu 22.04 through GitHub Actions. The new Linux packager and dependency-audit code are tested separately. A complete Linux binary bundle has not yet been built or validated: Linux-native engine and Python runtime inputs will be provided from the user's Linux machine or existing installation archive. The earlier Intel Mac package remains available but its Mach-O binaries cannot be used in Linux.

A source-tree run on a host without CCP4 cannot use Phaser, Buccaneer or REFMAC. Density modification is unaffected because it is in-tree (oasis-dm v2) and needs only the bundled Python. If the in-tree engine cannot start, the CLI prints the reason to stderr, the window shows a red warning, the result page marks the engine as a fallback, and the run record keeps `dm_engine_requested`/`dm_engine_fallback`.

The 2026-10-22 Linux-preparation changes passed 197 regression cases locally on macOS; that local result is not a Linux runtime validation. The GitHub workflow also checks the Linux application-only staging launcher.

## Build operator workflow

These are build-host steps, not end-user installation requirements. Work from staging copies; do not modify the original installation. Obtain a relocatable Linux x86_64 Python prefix containing Gemmi, NumPy and SciPy, and a Linux engine prefix containing Buccaneer, REFMAC, required libraries, dictionaries, startup definitions, license materials and reference resources. Their precise sourcing and relocation must be checked against the supplied Linux installation.

Audit the staged prefixes (requires existing Linux `readelf` and `ldd`):

```bash
python3 -B packaging/audit_linux_runtime.py --root /staging/python --kind python
python3 -B packaging/audit_linux_runtime.py --root /staging/runtime --kind engines
python3 -B packaging/build_standalone.py \
  --python-source /staging/python \
  --runtime-source /staging/runtime \
  --output /output/OASIS4.3-Build-linux-x86_64 \
  --tar /output/OASIS4.3-Build-linux-x86_64.tar.gz
```

The audit refuses missing libraries, non-system libraries loaded from outside the prefix, absolute or escaping RPATHs, broken symlinks and non-x86_64 ELF files. It records glibc version requirements from symbol metadata; this does not by itself certify every Linux distribution. Build on the oldest selected supported glibc baseline and test the resulting bundle on clean target systems. No exact minimum distribution version is promised until native binaries are examined and tested.

The assembler checks ELF architecture, manifest hashes, bundled Python/CLI startup and the doctor result. Linux uses `.tar.gz` to preserve executable modes and relative symlinks. Destinations must be new paths. No build tool installs dependencies automatically. Audit the final relocated and extracted copy, including every ELF binary/extension, before delivery.

## End-user launch (after a verified bundle is delivered)

```bash
tar -xzf OASIS4.3-Build-linux-x86_64.tar.gz
cd OASIS4.3-Build-linux-x86_64
./OASIS4.3-Build.sh
```

The window is a native Tk desktop window; there is no web server, no port and no browser. On a server without a graphical session, use the command line instead — it takes the same options as the form:

```bash
./phasebuild sad-to-model --help
./phasebuild sad-to-model --sad SAD.mtz --native NATIVE.mtz --sequence SEQ.fasta \
   --heavy HEAVY.pdb --wavelength 0.97918 --elements Se --output NEW_DIR
```

For a window over SSH, forward X11 (`ssh -X user@server`) and run `./phasebuild gui`; no public listener is opened by the program.

## Release gate

Verify clean-environment GUI startup, all six GUI mode plans, engine loading, synthetic construction/refinement calls, original-installation isolation, relocation and artifact checksums. Real-data structure-solving validation is separate. Public redistribution conditions for engines, libraries and assets must also be resolved before a public binary release. This source update does not create or authorize that release.
