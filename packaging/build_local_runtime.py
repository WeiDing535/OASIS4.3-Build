"""Build a private macOS runtime from the user's existing CCP4 installation.

This does not create or authorize a public redistribution. It preserves the
installed binaries and complete lib/data tree byte-for-byte, records provenance
and copied license materials, and only performs non-structural startup probes.
"""
import argparse
import concurrent.futures
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import tempfile
import time


REPOSITORY = Path(__file__).resolve().parents[1]
SYSTEM_PREFIXES = ('/usr/lib/', '/System/Library/')
STARTUP_DEFINITIONS = ('environ.def', 'default.def')


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def dependency_names(path, run=subprocess.run):
    listing = run(['/usr/bin/otool', '-L', str(path)], capture_output=True, text=True, check=True)
    identity = run(['/usr/bin/otool', '-D', str(path)], capture_output=True, text=True, check=True)
    own_ids = {line.strip() for line in identity.stdout.splitlines()[1:] if line.strip()}
    names = []
    for line in listing.stdout.splitlines():
        if not line.startswith('\t'):
            continue
        name = line.strip().split(' (compatibility version', 1)[0]
        if name not in own_ids:
            names.append(name)
    return names, {'otool_L': listing.stdout, 'otool_D': identity.stdout,
                   'dylib_install_ids_excluded': sorted(own_ids)}


def resolve_dependency(name, loader, source):
    if name.startswith(SYSTEM_PREFIXES):
        return None
    if name.startswith('@loader_path/'):
        target = loader.parent / name[len('@loader_path/'):]
    elif name.startswith('@executable_path/'):
        target = source / 'bin' / name[len('@executable_path/'):]
    elif name.startswith('@rpath/'):
        # Do not silently pick a host library and then leave an unsafe rpath in
        # the bundle. The installed target binaries use loader-relative edges.
        raise RuntimeError('An @rpath dependency requires an explicit relocatability review: ' + name)
    elif name.startswith('/'):
        raise RuntimeError('Non-system absolute dependency requires explicit relocation: ' + name)
    else:
        raise RuntimeError('Unsupported dependency load path: ' + name)
    target = Path(os.path.abspath(target))
    try:
        target.relative_to(source)
    except ValueError:
        raise RuntimeError('Dependency escapes the installed source prefix: ' + str(target))
    if not target.is_file():
        raise RuntimeError('Missing installed dependency: ' + str(target))
    return target


def runtime_graph(source):
    pending = [source / 'bin/cbuccaneer', source / 'bin/refmac5']
    seen, graph = set(), []
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        if not path.is_file():
            raise RuntimeError('Missing installed executable/library: ' + str(path))
        names, evidence = dependency_names(path)
        edges = []
        for name in names:
            resolved = resolve_dependency(name, path, source)
            edge = {'load_path': name, 'system': resolved is None}
            if resolved is not None:
                edge['relative_target'] = resolved.relative_to(source).as_posix()
                if resolved not in seen:
                    pending.append(resolved)
            edges.append(edge)
        architecture = subprocess.run(['/usr/bin/file', str(path)], capture_output=True,
                                      text=True, check=True).stdout.strip()
        if 'x86_64' not in architecture:
            raise RuntimeError('Expected installed macOS x86_64 binary: ' + architecture)
        graph.append({'relative_path': path.relative_to(source).as_posix(),
                      'dependencies': edges, 'load_command_evidence': evidence,
                      'architecture_evidence': architecture})
    return sorted(graph, key=lambda item: item['relative_path'])


def tree_files(folder):
    # Refuse directory symlinks rather than silently omitting data or following
    # a loop. The audited installed lib/data and license trees have none.
    result = []
    for parent, directories, files in os.walk(folder, followlinks=False):
        for name in directories:
            path = Path(parent) / name
            if path.is_symlink():
                raise RuntimeError('Directory symlink needs explicit portable staging: ' + str(path))
        result.extend(Path(parent) / name for name in sorted(files))
    return sorted(result)


def copy_verified(source, destination, role, output_root, license_group):
    before = source.stat()
    source_hash = sha256(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Dereference file symlinks so every copied load-path name is self-contained.
    shutil.copy2(source, destination)
    copied_hash = sha256(destination)
    after = source.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError('Installed source changed while copying: ' + str(source))
    if source_hash != copied_hash:
        raise RuntimeError('Copied bytes differ from installed source: ' + str(source))
    relative_path = destination.relative_to(output_root).as_posix()
    return {'path': relative_path, 'relative_path': relative_path, 'sha256': copied_hash,
            'bytes': destination.stat().st_size, 'role': role,
            'source_path': str(source), 'source_resolved_path': str(source.resolve()),
            'source_sha256': source_hash, 'source_file_symlink': os.readlink(source) if source.is_symlink() else None,
            'copy_strategy': 'Unmodified exact bytes; file symlink dereferenced if present',
            'mode': oct(stat.S_IMODE(destination.stat().st_mode)),
            'license_evidence_group': license_group, 'public_redistribution_cleared': False}


def startup_definition_plan(source, output):
    """CCP4FYP reads these logical-name definitions through CINCL at startup."""
    plan = []
    for name in STARTUP_DEFINITIONS:
        original = source / 'include' / name
        if not original.is_file():
            raise RuntimeError('Missing installed runtime startup definition: ' + str(original))
        plan.append((original, output / 'include' / name,
                     'runtime_startup_logical_name_definition', 'installed_component_licenses'))
    return plan


def probe_runtime_environment(output, scratch):
    """Use the application's real environment constructor, not a test facsimile."""
    path = REPOSITORY / 'src/phasebuild/runtime.py'
    spec = importlib.util.spec_from_file_location('phasebuild_runtime_probe_environment', path)
    runtime = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runtime)
    environment = runtime.runtime_environment(output,
        inherited={'PATH': '/usr/bin:/bin', 'LANG': 'C', 'LC_ALL': 'C'}, scratch=scratch)
    return environment, {'path': str(path), 'sha256': sha256(path)}


NOTICES = '''# Third-party notices: local evaluation runtime

This is a LOCAL-EVALUATION bundle assembled from the user's existing installed
CCP4 copy for this user's evaluation. It is not an authorized public release,
does not grant new rights to any component, and must not be described as a
publicly redistributable PhaseBuild binary package.

The original installation is not modified. Binary/library/data bytes are copied
without alteration. The complete installed lib/data tree is retained, including
the complete monomer chemical dictionary, fixed 1TQW reference and neural-network
weights. No dictionary, scattering factors or scientific model is substituted.

## Components and attribution

- Buccaneer: Kevin Cowtan and contributors; installed cbuccaneer 1.6.12 is the
  tracing engine. Historical Coot fragments with LGPLv3+ headers do not establish
  a redistribution license for this complete installed executable.
- REFMAC5: Garib Murshudov, Alexei Vagin and contributors; installed 5.8.0431.
  No matching complete source/build/redistribution grant was established here.
- Clipper, MMDB2, CCP4 C/Fortran IO, FFTW2, GNU Fortran/GCC/quadmath runtimes and
  their authors supply the libraries recorded in runtime_manifest.json.
- The original monomer dictionary, reference structures, pretrained weights and
  other installed data retain their original attribution and conditions.

## License materials and unresolved release conditions

The installed share/licenses tree is copied intact under share/licenses.
The installed conditions.txt and academic licence PDF/RTF are copied under
share/licenses/installed_ccp4_agreement. These documents are evidence, not a new
license grant. Local documentation includes LGPL/GPL component-specific files;
in particular MMDB2 is listed under GPLv3 and FFTW2 under GPLv2-or-later. Do not
replace the complete dependency analysis with the license of one source file.

The official CCP4 Academic Software Licence (UKRI, July 2021), clauses2.2.1 and
2.2.2, restrict third-party distribution/sharing of CCP4 Applications and copying
or modification except as specified; GNU and third-party components have their
own terms under clauses2.5/2.6. The installed older agreement likewise contains
application restrictions. This local evaluation copy is NOT a determination of
which separate component license or permission would authorize a public package.

Outstanding before any public redistribution: establish the precise cbuccaneer
and REFMAC executable permissions; obtain matching corresponding sources and
build recipes where required; verify Clipper/MMDB/CCP4 IO/FFTW/GCC runtime
licenses and exceptions against the exact versions; review all data/weights;
provide the corresponding sources, notices and installation information required
by the eventual distribution mode. No external party was contacted here.

## Official references

- https://www.ccp4.ac.uk/ccp4-licence-and-conditions-2/
- https://www.ccp4.ac.uk/licensing/CCP4_Academic_Software_Licence_UKRI_July_2021.pdf
- https://www.ccp4.ac.uk/licensing/GPL_and_LGPL_software.html
- https://fftw.org/fftw2_doc/fftw_8.html
- https://www.ccp4.ac.uk/html/cbuccaneer.html
- https://www.ccp4.ac.uk/html/refmac5/description.html

The inspected macOS files are x86_64. Apple system libraries/Accelerate are
provided by macOS and are not copied. Other platforms require their own builds.
Startup checks establish loading/banners only; structural quality requires the
separate fixed-protocol validation performed by the application project.
'''


def startup_probe(output, name, probe_parent=None):
    probe = (output / 'share/startup_checks' if probe_parent is None else probe_parent) / name
    probe.mkdir(parents=True, exist_ok=False)
    with tempfile.TemporaryDirectory(prefix='phasebuild_runtime_startup_') as temporary:
        home = Path(temporary)
        environment, environment_source = probe_runtime_environment(output, home / 'scratch')
        environment['DYLD_PRINT_LIBRARIES'] = '1'  # Loading evidence only.
        command = [str(output / 'bin' / name), '-help' if name == 'cbuccaneer' else '-i']
        run = subprocess.run(command, input='END\n', capture_output=True, text=True,
                             cwd=str(home), env=environment, timeout=30)
    (probe / 'stdout.txt').write_text(run.stdout)
    (probe / 'stderr.txt').write_text(run.stderr)
    combined = run.stdout + run.stderr
    dynamic = [line for line in run.stderr.splitlines() if '.dylib' in line or 'dyld' in line]
    original_loaded = [line for line in dynamic if '/Applications/ccp4-9/' in line]
    loader_error = any(value in combined for value in
        ['Library not loaded:', 'Symbol not found:', 'Incompatible library version', 'wrong architecture'])
    initialization_error = any(value in combined for value in
        ['Cannot open environ.def', 'Cannot open default.def'])
    definitions = {name: (output / 'include' / name).is_file() for name in STARTUP_DEFINITIONS}
    banner = bool(re.search('buccaneer' if name == 'cbuccaneer' else 'refmac', combined, re.I))
    result = {'program': name, 'command': ['${RUNTIME_ROOT}/bin/' + name, command[-1]],
        'environment_policy': 'Application runtime_environment from PATH/LANG/LC_ALL only; private scratch supplied through its public scratch parameter',
        'environment_constructor': environment_source,
        'environment': environment,
        'exit_code': run.returncode, 'banner_present': banner,
        'original_installation_libraries_loaded': original_loaded,
        'loader_failure': loader_error, 'structural_experiment': False,
        'runtime_startup_definition_files_present': definitions,
        'runtime_initialization_failure': initialization_error,
        'expected_no_input_behavior': 'Nonzero usage/missing-input exits can be normal; no reflection/model input was provided',
        'passed': banner and not original_loaded and not loader_error and not initialization_error
                  and all(definitions.values()) and run.returncode >= 0}
    (probe / 'result.json').write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=Path('/Applications/ccp4-9'))
    parser.add_argument('--output', type=Path, default=REPOSITORY / 'dist/local_macos_x86_64/runtime')
    parser.add_argument('--workers', type=int, default=8)
    parser.add_argument('--skip-startup-check', action='store_true')
    arguments = parser.parse_args()
    source, output = arguments.source.resolve(), arguments.output.resolve()
    if output.exists():
        raise RuntimeError('Output already exists; never overwrite an existing evaluation runtime: ' + str(output))
    try:
        output.relative_to(source)
    except ValueError:
        pass
    else:
        raise RuntimeError('Output must not be inside the original installation')
    if arguments.workers < 1 or arguments.workers > 32:
        raise ValueError('Workers must be1..32')
    started = time.time()
    graph = runtime_graph(source)
    plan = []
    for item in graph:
        relative = Path(item['relative_path'])
        plan.append((source / relative, output / relative, 'executable' if relative.parts[0] == 'bin' else 'dynamic_library', 'installed_component_licenses'))
    plan.extend(startup_definition_plan(source, output))
    for path in tree_files(source / 'lib/data'):
        plan.append((path, output / path.relative_to(source), 'scientific_data_complete_installed_tree', 'installed_data_licenses'))
    for path in tree_files(source / 'share/licenses'):
        plan.append((path, output / path.relative_to(source), 'license_material', 'license_text'))
    for name in ['conditions.txt', 'academic_software_licence.pdf', 'academic_software_licence.rtf']:
        path = source / name
        if path.is_file():
            plan.append((path, output / 'share/licenses/installed_ccp4_agreement' / name, 'installed_ccp4_agreement', 'license_text'))
    for name in ['SOURCE_LICENSE_AUDIT.md', 'source_license_audit.json']:
        path = REPOSITORY / 'research/standalone_20261019' / name
        if path.is_file():
            plan.append((path, output / 'share/provenance' / name, 'local_license_review_evidence', 'not_a_license_grant'))
    assert len({str(item[1]) for item in plan}) == len(plan)
    output.mkdir(parents=True)
    manifest = {'schema_version': 1, 'runtime_kind': 'local-evaluation',
        'platform': 'macos-x86_64', 'local_evaluation_only': True,
        'public_redistribution_authorized': False, 'corresponding_sources_complete': False,
        'source_installation': str(source), 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'builder': {'path': str(Path(__file__).resolve()), 'sha256': sha256(__file__)},
        'layout': {'bin': 'bin', 'libraries': 'lib', 'clibd': 'lib/data', 'monomer_library': 'lib/data/monomers',
                   'runtime_startup_definitions': 'include'},
        'environment_templates': {'CCP4': '${RUNTIME_ROOT}', 'CLIBD': '${RUNTIME_ROOT}/lib/data',
                                  'CLIBD_MON': '${RUNTIME_ROOT}/lib/data/monomers/',
                                  'CINCL': '${RUNTIME_ROOT}/include',
                                  'CCP4_SCR': '${WRITABLE_JOB_SCRATCH}'},
        'data_policy': 'Complete installed lib/data and complete monomer dictionary, unmodified bytes; no scientific data pruning',
        'no_original_installation_modified': True, 'binary_load_commands_modified': False,
        'dependency_graph': graph, 'files': [], 'startup_checks': [], 'build_complete': False,
        'missing_corresponding_source_and_permissions': [
            'Exact installed cbuccaneer1.6.12 source/build/redistribution terms not established',
            'Exact installed REFMAC5.8.0431 source/build/redistribution terms not established',
            'Exact library corresponding sources and GCC runtime license/exception materials not supplied',
            'Public redistribution rights for complete data/weights/assets not established'],
        'license_materials': {'installed_directory': 'share/licenses',
            'installed_agreements': 'share/licenses/installed_ccp4_agreement',
            'local_review': 'share/provenance', 'notice': 'THIRD_PARTY_NOTICES.md'},
        'inventory_excludes': ['runtime_manifest.json', 'runtime_manifest.sha256'],
        'inventory_exclusion_reason': 'Self-referential manifest/checksum bookkeeping; every other bundled regular file is inventoried'}
    try:
        print('Copying and verifying %d installed files' % len(plan), flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=arguments.workers) as executor:
            futures = [executor.submit(copy_verified, original, destination, role, output, group)
                       for original, destination, role, group in plan]
            for count, done in enumerate(concurrent.futures.as_completed(futures), 1):
                manifest['files'].append(done.result())
                if count % 4000 == 0:
                    print('Verified %d/%d files' % (count, len(plan)), flush=True)
        (output / 'THIRD_PARTY_NOTICES.md').write_text(NOTICES)
        (output / 'README_LOCAL_EVALUATION.txt').write_text(
            'Private local evaluation runtime; not a public redistribution.\n'
            'Resolve this directory dynamically; bin/ and lib/ are relocatable together.\n'
            'Set CCP4=<runtime>, CLIBD=<runtime>/lib/data, CLIBD_MON=<runtime>/lib/data/monomers/.\n'
            'Set CINCL=<runtime>/include for CCP4FYP logical-name startup definitions.\n'
            'Set CCP4_SCR to an existing writable job-specific scratch directory.\n'
            'Complete installed monomer dictionary and data were copied unchanged.\n'
            'See runtime_manifest.json and THIRD_PARTY_NOTICES.md before use or release.\n')
        if not arguments.skip_startup_check:
            for name in ['cbuccaneer', 'refmac5']:
                result = startup_probe(output, name)
                manifest['startup_checks'].append(result)
                print('Startup %s: %s' % (name, result['passed']), flush=True)
                if not result['passed']:
                    raise RuntimeError('Bundled startup check failed: ' + name)
        copied_paths = {item['path'] for item in manifest['files']}
        for path in sorted(output.rglob('*')):
            if path.is_file() and path.relative_to(output).as_posix() not in copied_paths:
                relative_path = path.relative_to(output).as_posix()
                manifest['files'].append({'path': relative_path, 'relative_path': relative_path,
                    'sha256': sha256(path), 'bytes': path.stat().st_size,
                    'role': 'generated_documentation_or_startup_evidence', 'source_path': None,
                    'generated_by_builder_sha256': manifest['builder']['sha256'],
                    'public_redistribution_cleared': False})
        manifest['files'].sort(key=lambda item: item['path'])
        # Verify closure paths exist under this root; no system library is copied.
        assert all((output / entry['relative_path']).is_file() and all(
            edge['system'] or (output / edge['relative_target']).is_file()
            for edge in entry['dependencies']) for entry in graph)
        monomer_relative = {str(path.relative_to(source)) for path in tree_files(source / 'lib/data/monomers')}
        assert monomer_relative == {item['path'] for item in manifest['files'] if item['path'].startswith('lib/data/monomers/')}
        manifest.update(build_complete=True, copied_file_count=len(plan), inventory_file_count=len(manifest['files']),
            monomer_file_count=len(monomer_relative),
            total_bytes=sum(item['bytes'] for item in manifest['files']), seconds=time.time() - started)
    except Exception as error:
        manifest.update(build_error=str(error), seconds=time.time() - started)
        raise
    finally:
        target = output / 'runtime_manifest.json'
        target.write_text(json.dumps(manifest, indent=2) + '\n')
        (output / 'runtime_manifest.sha256').write_text(sha256(target) + '  runtime_manifest.json\n')
    print(json.dumps({key: manifest[key] for key in ['build_complete', 'copied_file_count',
        'inventory_file_count', 'monomer_file_count', 'total_bytes', 'seconds']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
