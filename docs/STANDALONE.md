# Intel Mac complete local package

The user-facing packaging goal is one application folder containing the Python environment, scientific Python packages, external Buccaneer/REFMAC executables, dynamic libraries, monomer dictionary and reference resources. Users supply experimental inputs and do not install CCP4 separately.

`packaging/build_standalone.py` combines already provisioned runtime directories with current application sources. It verifies source and copied file hashes, rejects escaping links and existing destinations, runs the bundled doctor and checks CLI entry points. It does not install dependencies or modify original installations.

The current build is Intel macOS only and is a local-evaluation package. It is not a publicly redistributable binary release: exact external executable, library, data and weight distribution conditions remain unresolved. Source publication and binary-runtime publication are separate decisions. The bundled license texts and notices remain intact.

Local builds can be made using an existing verified portable runtime:

```bash
python3 -B packaging/build_standalone.py \
  --python-source /path/to/verified/python \
  --runtime-source /path/to/verified/runtime \
  --output /path/to/new/OASIS4.3-Build \
  --zip /path/to/new/OASIS4.3-Build.zip
```

Double-click `OASIS4.3-Build.command` in the extracted directory. Keep the entire directory layout together and use a writable location. All six GUI modes use the bundled runtime; the GUI retains PhaseBuild as its internal engine name. Development benchmark scripts may additionally require their own datasets and explicit paths; target datasets are not bundled. The package has no Apple signing/notarization and has not been validated on a fresh second machine.

Acceptance checks include manifest verification, clean-environment startup, denial of reads from the original CCP4 installation, engine dynamic loading, synthetic-data calls, and GUI runtime selection. Startup and synthetic data checks do not establish scientific success on real experimental data.

## Local verification on 2026-10-22

189 source regression cases passed. With reads from the original CCP4 installation and previous project copies denied by an OS sandbox, the bundled doctor, native trace, standalone R factors, amplitude-only and HL-aware refinement adapters, GUI runtime selection, one real REFMAC refinement cycle and one real Buccaneer construction cycle passed on synthetic data. These are runtime/integration checks, not real-data structure-solving validation.
