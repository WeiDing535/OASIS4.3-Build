#!/usr/bin/env python3
"""Build a local relocatable macOS scientific Python, without editing CCP4.

This assembles already installed binaries; it is not a source build and does not
make a redistribution-rights determination.  Run with any Python >= 3.9 and
macOS developer command-line tools.  The destination must not already exist.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time


SYSTEM_PREFIXES = ("/usr/lib/", "/System/Library/")
MACHO_MAGICS = {b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe", b"\xfe\xed\xfa\xcf",
                b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe", b"\xbe\xba\xfe\xca"}
EXCLUDED_STDLIB_DIRS = {"site-packages", "lib-dynload", "test", "idlelib", "tkinter",
                         "turtledemo", "ensurepip", "config-3.9-darwin", "__pycache__"}


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def command(args, **kwargs):
    result = subprocess.run([str(x) for x in args], text=True, capture_output=True,
                            **kwargs)
    if result.returncode:
        raise RuntimeError("Command failed (%d): %r\nstdout:\n%s\nstderr:\n%s" %
                           (result.returncode, result.args, result.stdout, result.stderr))
    return result


def is_macho(path):
    if not path.is_file():
        return False
    with path.open("rb") as f:
        return f.read(4) in MACHO_MAGICS


def ignore_python(directory, names):
    return [n for n in names if n == "__pycache__" or n in {"tests", "test"}
            or n.endswith((".pyc", ".pyo"))]


def dylib_commands(path):
    """Parse load commands, distinguishing LC_ID_DYLIB from dependencies."""
    out = command(["/usr/bin/otool", "-l", path]).stdout
    kind = None
    result = {"dependencies": [], "rpaths": [], "id": None}
    for line in out.splitlines():
        fields = line.strip().split()
        if len(fields) > 1 and fields[0] == "cmd":
            kind = fields[1]
        elif fields and fields[0] == "name" and " (offset " in line:
            name = line.strip()[5:].rsplit(" (offset ", 1)[0]
            if kind == "LC_ID_DYLIB":
                result["id"] = name
            elif kind in {"LC_LOAD_DYLIB", "LC_LOAD_WEAK_DYLIB", "LC_REEXPORT_DYLIB",
                          "LC_LAZY_LOAD_DYLIB", "LC_LOAD_UPWARD_DYLIB"}:
                result["dependencies"].append(name)
        elif fields and fields[0] == "path" and kind == "LC_RPATH":
            result["rpaths"].append(line.strip()[5:].rsplit(" (offset ", 1)[0])
    return result


PROBE = r'''
import bz2, ctypes, decimal, hashlib, importlib, json, lzma, os, pathlib, sqlite3, ssl, sys
import numpy as np
import scipy
from scipy import fft, linalg, ndimage, optimize, spatial, special
import gemmi
a = np.array([[3., 1.], [1., 2.]])
assert np.allclose(np.linalg.solve(a, [1., 0.]), [0.4, -0.2])
assert np.allclose(linalg.eigvalsh(a), [1.381966011250105, 3.618033988749895])
assert abs(spatial.cKDTree([[0., 0., 0.]]).query([1., 0., 0.])[0] - 1.) < 1.e-12
assert np.allclose(ndimage.map_coordinates(np.arange(9.).reshape(3, 3), [[1.], [1.]]), [4.])
assert optimize.minimize(lambda x: (x[0] - 2.)**2, [0.]).success
assert abs(special.expit(0.) - 0.5) < 1.e-12
assert np.allclose(fft.irfft(fft.rfft(np.arange(8.))), np.arange(8.))
grid = gemmi.FloatGrid(4, 4, 4)
grid.unit_cell = gemmi.UnitCell(10, 10, 10, 90, 90, 90)
root = pathlib.Path(os.environ['PHASEBUILD_PYTHON_ROOT']).resolve()
def inside(path):
    try:
        pathlib.Path(path).resolve().relative_to(root)
        return True
    except ValueError:
        return False
modules = {name: str(pathlib.Path(mod.__file__).resolve())
           for name, mod in list(sys.modules.items()) if getattr(mod, '__file__', None)}
outside = {name: p for name, p in modules.items() if not inside(p)}
assert not outside, outside
lib = ctypes.CDLL(None)
lib._dyld_image_count.restype = ctypes.c_uint32
lib._dyld_get_image_name.argtypes = [ctypes.c_uint32]
lib._dyld_get_image_name.restype = ctypes.c_char_p
loaded = [lib._dyld_get_image_name(i).decode() for i in range(lib._dyld_image_count())]
outside_libraries = [p for p in loaded if not inside(p)
    and not p.startswith(('/usr/lib/', '/System/Library/', '/System/Volumes/Preboot/Cryptexes/OS/'))]
assert not outside_libraries, outside_libraries
result = {'executable': sys.executable, 'prefix': sys.prefix, 'base_prefix': sys.base_prefix,
          'sys_path': sys.path, 'python_version': sys.version,
          'packages': {name: {'version': getattr(mod, '__version__', None), 'file': mod.__file__}
                      for name, mod in [('numpy', np), ('scipy', scipy), ('gemmi', gemmi)]},
          'module_paths': modules, 'loaded_images': loaded,
          'outside_module_paths': outside, 'outside_non_system_images': outside_libraries,
          'functional_checks': ['numpy_solve', 'scipy_eigvalsh', 'ckdtree', 'map_coordinates',
                                'minimize', 'special_expit', 'fft_roundtrip', 'gemmi_grid'],
          'no_ccp4_runtime_dependency': not any('/Applications/ccp4' in p for p in list(modules.values())+loaded)}
print(json.dumps(result, indent=2, sort_keys=True))
'''


def build(args):
    source = args.framework_version.resolve()
    architectures = command(["/usr/bin/lipo", "-archs", source / "Python"]).stdout.strip().split()
    if architectures != ["x86_64"]:
        raise SystemExit("This builder currently validates x86_64 only: " + repr(architectures))
    ccp4_root = source.parents[3]
    dest = args.output.resolve()
    if dest.exists():
        raise SystemExit("Destination already exists; preserve it and choose another --output: " + str(dest))
    dest.parent.mkdir(parents=True, exist_ok=True)
    stage = dest.parent / (".python-build-" + str(os.getpid()))
    stage.mkdir()
    log = []
    manifest = {"format_version": 1, "created_unix": time.time(),
                "purpose": "local development runtime; distribution clearance is separate",
                "platform": "macos_x86_64", "source_framework_version": str(source),
                "architectures": architectures,
                "entrypoint": "bin/python3", "files": [], "patches": [],
                "third_party_licenses": [], "limitations": [
                    "Copied from an existing installation, not reproducibly rebuilt from source.",
                    "Same-host environment-isolation and physical-relocation tests; no clean-machine test.",
                    "macOS x86_64 only; no native Apple Silicon, Linux or Windows claim.",
                    "No pip, development headers, test suites, Tk, IDLE, or CCP4-specific Python modules.",
                    "External library redistribution and source obligations require separate review."]}
    mapping = {}
    source_hash = {}

    def remember(src, dst):
        src = Path(src).resolve()
        dst = Path(dst)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        dst.chmod(dst.stat().st_mode | 0o200)
        mapping[src] = dst
        source_hash[dst] = {"source": str(src), "source_sha256": sha256(src)}

    version = stage / "Frameworks/Python.framework/Versions/3.9"
    remember(source / "Python", version / "Python")
    remember(source / "bin/python3.9", version / "bin/python3.9")
    (version / "Resources").mkdir()
    remember(source / "Resources/Info.plist", version / "Resources/Info.plist")
    # Framework builds launch the application executable from the bin stub.
    # It is part of CPython, even for command-line operation.
    shutil.copytree(source / "Resources/Python.app", version / "Resources/Python.app",
                    ignore=ignore_python)
    stdlib = source / "lib/python3.9"
    outlib = version / "lib/python3.9"
    outlib.mkdir(parents=True)
    for child in stdlib.iterdir():
        if child.name in EXCLUDED_STDLIB_DIRS or child.suffix in {".pyc", ".pyo"}:
            continue
        if child.is_dir():
            shutil.copytree(child, outlib / child.name, ignore=ignore_python)
        elif child.is_file():
            remember(child, outlib / child.name)
    dynload = outlib / "lib-dynload"
    dynload.mkdir()
    for child in sorted((stdlib / "lib-dynload").glob("*.cpython-39-darwin.so")):
        if child.name.startswith(("_test", "_tkinter", "xxlimited")):
            continue
        remember(child, dynload / child.name)
    site = outlib / "site-packages"
    site.mkdir()
    for package in ("numpy", "scipy", "gemmi"):
        shutil.copytree(stdlib / "site-packages" / package, site / package, ignore=ignore_python)
        for metadata in (stdlib / "site-packages").glob(package + "-*"):
            if metadata.is_dir() and metadata.name.endswith((".dist-info", ".egg-info")):
                shutil.copytree(metadata, site / metadata.name, ignore=ignore_python)
    # Map all copied package/stdlib Mach-O files back to their original location.
    for copied in stage.rglob("*"):
        if copied.is_file() and copied not in source_hash:
            rel = copied.relative_to(version)
            original = source / rel
            source_hash[copied] = {"source": str(original), "source_sha256": sha256(original)}
            if is_macho(copied):
                mapping[original.resolve()] = copied

    def resolve(dep, owner, rpaths):
        if dep.startswith("@loader_path/"):
            candidate = owner.parent / dep[len("@loader_path/"):]
        elif dep.startswith("@executable_path/"):
            candidate = source / "bin" / dep[len("@executable_path/"):]
        elif dep.startswith("@rpath/"):
            candidates = []
            for rp in rpaths:
                rp = rp.replace("@loader_path", str(owner.parent))
                rp = rp.replace("@executable_path", str(source / "bin"))
                candidates.append(Path(rp) / dep[len("@rpath/"):])
            candidates += [ccp4_root / "lib" / Path(dep).name,
                           ccp4_root / "opt/lib" / Path(dep).name]
            candidate = next((p for p in candidates if p.exists()), None)
            if candidate is None:
                raise RuntimeError("Unresolved @rpath dependency %s of %s" % (dep, owner))
        else:
            candidate = Path(dep)
        if not candidate.exists():
            raise RuntimeError("Missing dependency %s of %s -> %s" % (dep, owner, candidate))
        return candidate.resolve()

    queue = list(mapping)
    plans = {}
    for original in queue:
        copied = mapping[original]
        loads = dylib_commands(original)
        changes = []
        for dep in loads["dependencies"]:
            if dep.startswith(SYSTEM_PREFIXES):
                continue
            dep_source = resolve(dep, original, loads["rpaths"])
            if dep_source not in mapping:
                dep_dest = stage / "lib" / dep_source.name
                if dep_dest.exists() and sha256(dep_dest) != sha256(dep_source):
                    raise RuntimeError("Conflicting library basenames: " + str(dep_source))
                remember(dep_source, dep_dest)
                queue.append(dep_source)
            relative = os.path.relpath(mapping[dep_source], copied.parent)
            changes.append((dep, "@loader_path/" + relative))
        plans[original] = (loads, changes)
    # Resolve the entire closure before changing any source-derived binary.
    for original, (loads, changes) in plans.items():
        copied = mapping[original]
        cmd = ["/usr/bin/install_name_tool"]
        for old, new in changes:
            cmd += ["-change", old, new]
        if loads["id"] is not None:
            cmd += ["-id", "@loader_path/" + copied.name]
        for rp in loads["rpaths"]:
            cmd += ["-delete_rpath", rp]
        if len(cmd) > 1:
            log.append(command(cmd + [copied]).stderr)
        # New ad-hoc signatures describe modified bytes; they confer no publisher trust.
        signed = command(["/usr/bin/codesign", "--force", "--sign", "-", copied])
        log.append(signed.stderr)
        manifest["patches"].append({"file": str(copied.relative_to(stage)),
                                    "dependency_changes": changes,
                                    "removed_rpaths": loads["rpaths"], "ad_hoc_signed": True})
        final_loads = dylib_commands(copied)
        for dep in final_loads["dependencies"]:
            if dep.startswith(SYSTEM_PREFIXES):
                continue
            if not dep.startswith("@loader_path/"):
                raise RuntimeError("Non-relative dependency remains: " + dep)
            resolved = (copied.parent / dep[len("@loader_path/"):]).resolve()
            resolved.relative_to(stage.resolve())
            if not resolved.exists():
                raise RuntimeError("Missing bundled dependency: " + str(resolved))

    launcher = stage / "bin/python3"
    launcher.parent.mkdir()
    launcher.write_text('''#!/bin/sh
set -eu
_phasebuild_python_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd -P)
unset PYTHONHOME PYTHONSTARTUP PYTHONUSERBASE DYLD_LIBRARY_PATH DYLD_FALLBACK_LIBRARY_PATH DYLD_INSERT_LIBRARIES
export PYTHONHOME="$_phasebuild_python_root/Frameworks/Python.framework/Versions/3.9"
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1
exec "$PYTHONHOME/bin/python3.9" "$@"
''')
    launcher.chmod(0o755)
    (stage / "bin/python3.9").symlink_to("python3")
    licenses = stage / "licenses"
    licenses.mkdir()
    license_sources = {
        "Python-LICENSE.txt": stdlib / "LICENSE.txt",
        "NumPy-LICENSE.txt": stdlib / "site-packages/numpy/LICENSE.txt",
        "SciPy-LICENSE.txt": stdlib / "site-packages/scipy/LICENSE.txt",
        "Gemmi-LICENSE.txt": stdlib / "site-packages/gemmi-0.7.5.dist-info/licenses/LICENSE.txt",
    }
    for name, src in license_sources.items():
        remember(src, licenses / name)
        manifest["third_party_licenses"].append({"file": "licenses/" + name, "source": str(src)})
    (licenses / "DEPENDENCY_REVIEW.md").write_text(
        "# Local scientific Python runtime\n\n"
        "Python, NumPy, SciPy and Gemmi license texts are included. Recursive binary "
        "dependencies, source versions and hashes are in python_manifest.json. "
        "This directory is a local development bundle, not cleared for public redistribution. "
        "The surrounding distribution audit must also cover the bundled BLAS/LAPACK, "
        "GCC runtime libraries, OpenSSL and other standard-library dependencies, "
        "including any applicable source/relinking obligations. No licenses are inferred "
        "from filenames. Apple system libraries are dynamically linked, not copied.\n")
    (stage / "build.log").write_text("".join(log))
    (stage / "probe.py").write_text(PROBE)
    stage.rename(dest)
    env = {"PATH": "/usr/bin:/bin:/usr/sbin:/sbin", "LANG": "C", "HOME": str(dest / ".isolated-home"),
           "PHASEBUILD_PYTHON_ROOT": str(dest), "DYLD_PRINT_LIBRARIES": "1"}
    for label, root in [("direct", dest)]:
        result = subprocess.run([str(root / "bin/python3"), str(root / "probe.py")],
                                cwd="/private/tmp", env=env, capture_output=True, text=True)
        (dest / ("probe_" + label + ".json")).write_text(result.stdout)
        (dest / ("probe_" + label + ".dyld.log")).write_text(result.stderr)
        result.check_returncode()
        manifest[label + "_probe"] = json.loads(result.stdout)
    # Physical relocation, including a space in the directory name, not a symlink test.
    with tempfile.TemporaryDirectory(prefix="phasebuild-python-relocation-") as temporary:
        relocated = Path(temporary) / "renamed scientific python"
        shutil.copytree(dest, relocated, symlinks=True)
        env["PHASEBUILD_PYTHON_ROOT"] = str(relocated)
        env["HOME"] = str(relocated / ".isolated-home")
        result = subprocess.run([str(relocated / "bin/python3"), str(relocated / "probe.py")],
                                cwd="/private/tmp", env=env, capture_output=True, text=True)
        (dest / "probe_relocated.json").write_text(result.stdout)
        (dest / "probe_relocated.dyld.log").write_text(result.stderr)
        result.check_returncode()
        manifest["relocated_probe"] = json.loads(result.stdout)
    for copied in sorted(dest.rglob("*")):
        if not copied.is_file() or copied.is_symlink():
            continue
        old_key = stage / copied.relative_to(dest)
        record = {"file": str(copied.relative_to(dest)),
                  "relative_path": str(copied.relative_to(dest)), "size": copied.stat().st_size,
                  "sha256": sha256(copied), **source_hash.get(old_key, {})}
        manifest["files"].append(record)
    manifest["total_bytes"] = sum(record["size"] for record in manifest["files"])
    manifest["symlinks"] = [{"relative_path": str(p.relative_to(dest)), "target": os.readlink(p)}
                            for p in sorted(dest.rglob("*")) if p.is_symlink()]
    manifest["macho_count"] = len(plans)
    manifest["status"] = "local_runtime_verified"
    (dest / "python_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(dest), "entrypoint": str(dest / "bin/python3"),
                      "status": manifest["status"], "total_bytes": manifest["total_bytes"],
                      "macho_count": manifest["macho_count"],
                      "direct_import_paths": manifest["direct_probe"]["packages"],
                      "relocated_no_ccp4": manifest["relocated_probe"]["no_ccp4_runtime_dependency"]},
                     indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--framework-version", type=Path,
                        default=Path("/Applications/ccp4-9/Frameworks/Python.framework/Versions/3.9"))
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).resolve().parents[1] / "dist/local_macos_x86_64/python")
    build(parser.parse_args())
