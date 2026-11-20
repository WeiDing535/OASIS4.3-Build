#!/usr/bin/env python3
"""Build the OASIS4.3 source release archive from a Git commit.

The archive is produced by ``git archive``, so it contains exactly the files
tracked by that commit: nothing untracked, nothing ignored, no local run
output, no Git history. The file list and the tree object id are recorded in a
sidecar manifest so the published archive can be checked against the commit
that produced it.
"""
import argparse
import hashlib
import json
import subprocess
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREFIX = 'OASIS4.3'


def git(*args, text=True):
    out = subprocess.run(['git', '-C', str(ROOT), *args],
                         check=True, capture_output=True)
    return out.stdout.decode() if text else out.stdout


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--ref', default='HEAD', help='Commit-ish to archive')
    parser.add_argument('--version', required=True, help='External release number, e.g. 4.3.8.7')
    args = parser.parse_args(argv)
    target = args.output.resolve()
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists():
        target.unlink()
    tree = git('rev-parse', args.ref + '^{tree}').strip()
    commit = git('rev-parse', args.ref).strip()
    with tarfile.open(target, 'w:gz', format=tarfile.GNU_FORMAT) as archive:
        archive.add(str(ROOT), arcname=PREFIX, filter=lambda info: (
            info if info.isdir() else None))
        # git archive writes the authoritative tree; stream it into the file.
    with target.open('wb') as handle:
        subprocess.run(['git', '-C', str(ROOT), 'archive', '--format=tar.gz',
                        '--prefix=' + PREFIX + '/', args.ref],
                       check=True, stdout=handle)
    with tarfile.open(target, 'r:gz') as archive:
        members = archive.getmembers()
        files = [m for m in members if m.isfile()]
        names = sorted(m.name[len(PREFIX) + 1:] for m in files)
        digest = hashlib.sha256()
        entries = []
        for member in sorted(files, key=lambda m: m.name):
            payload = archive.extractfile(member).read()
            digest.update(member.name.encode())
            digest.update(b'\0')
            digest.update(hashlib.sha256(payload).digest())
            entries.append({'path': member.name[len(PREFIX) + 1:],
                            'bytes': len(payload),
                            'sha256': hashlib.sha256(payload).hexdigest()})
    summary = {'release': args.version, 'prefix': PREFIX, 'ref': args.ref,
               'commit': commit, 'tree': tree, 'files': len(entries),
               'tree_sha256': digest.hexdigest(),
               'bytes': target.stat().st_size, 'names': names, 'entries': entries}
    (target.parent / (target.name + '.manifest.json')).write_text(
        json.dumps(summary, indent=2) + '\n')
    print(json.dumps({k: v for k, v in summary.items()
                      if k not in ('entries', 'names')}, indent=2))
    return summary


if __name__ == '__main__':
    main()
