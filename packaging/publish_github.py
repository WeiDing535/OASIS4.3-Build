#!/usr/bin/env python3
"""Push a local commit to GitHub through the Git Data API, incrementally.

The remote branch holds an independent root commit, so a normal ``git push``
is rejected as non-fast-forward. This script creates the remote commit from
the local tree instead: it uploads only the blobs that changed relative to
``--base-tree``, builds a new tree on top of the remote tree, commits it and
moves the branch. The remote tree object id is then compared against the local
tree object id, so the content identity is proven rather than assumed.

Usage::

    python packaging/publish_github.py \
        --repo WeiDing535/OASIS4.3-Build --branch main \
        --message-file /tmp/commit_msg.txt --verify-tree <local tree sha>
"""
import argparse
import base64
import json
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
API = 'https://api.github.com'


def git(*args):
    out = subprocess.run(['git', '-C', str(ROOT), *args], check=True,
                         capture_output=True)
    return out.stdout.decode().strip()


def token():
    out = subprocess.run(['gh', 'auth', 'token'], check=True, capture_output=True)
    return out.stdout.decode().strip()


class Client:
    def __init__(self, repo):
        self.repo = repo
        self.token = token()

    def call(self, method, path, payload=None):
        url = API + '/repos/' + self.repo + path
        data = json.dumps(payload).encode() if payload is not None else None
        request = urllib.request.Request(url, data=data, method=method)
        request.add_header('Authorization', 'Bearer ' + self.token)
        request.add_header('Accept', 'application/vnd.github+json')
        request.add_header('X-GitHub-Api-Version', '2022-11-28')
        if data is not None:
            request.add_header('Content-Type', 'application/json')
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                body = response.read().decode()
                return json.loads(body) if body else {}
        except urllib.error.HTTPError as error:
            detail = error.read().decode()
            raise SystemExit('%s %s -> %s: %s' % (method, path, error.code, detail[:600]))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', required=True)
    parser.add_argument('--branch', default='main')
    parser.add_argument('--base-tree', required=True,
                        help='Tree of the commit the remote branch currently holds')
    parser.add_argument('--message-file', required=True)
    parser.add_argument('--verify-tree', required=True,
                        help='Local tree object id the remote tree must equal')
    parser.add_argument('--dry-run', action='store_true')
    args = parser.parse_args(argv)

    client = Client(args.repo)
    head = client.call('GET', '/git/ref/heads/' + args.branch)
    remote_commit = head['object']['sha']
    remote = client.call('GET', '/git/commits/' + remote_commit)
    remote_tree = remote['tree']['sha']
    print('remote commit %s tree %s' % (remote_commit, remote_tree))
    if remote_tree != args.base_tree:
        raise SystemExit('remote tree %s does not match the expected base %s'
                         % (remote_tree, args.base_tree))

    # Files that differ between the base commit and the local HEAD.
    changed = [line.split('\t', 1)[1]
               for line in git('diff', '--name-status', '1727c75', 'HEAD').splitlines()
               if line.strip()]
    deleted = [line.split('\t', 1)[1]
               for line in git('diff', '--name-status', '1727c75', 'HEAD').splitlines()
               if line.startswith('D\t')]
    print('%d changed paths (%d deleted)' % (len(changed), len(deleted)))
    if args.dry_run:
        for path in changed:
            print('  ', path)
        return

    entries = []
    for index, path in enumerate(changed, 1):
        if path in deleted:
            entries.append({'path': path, 'mode': '100644', 'type': 'blob', 'sha': None})
            continue
        payload = (ROOT / path).read_bytes()
        blob = client.call('POST', '/git/blobs',
                           {'content': base64.b64encode(payload).decode(),
                            'encoding': 'base64'})
        mode = '100755' if (ROOT / path).stat().st_mode & 0o111 else '100644'
        entries.append({'path': path, 'mode': mode, 'type': 'blob', 'sha': blob['sha']})
        print('  [%d/%d] %s' % (index, len(changed), path), flush=True)

    tree = client.call('POST', '/git/trees',
                       {'base_tree': remote_tree, 'tree': entries})
    print('new tree %s' % tree['sha'])
    if tree['sha'] != args.verify_tree:
        raise SystemExit('remote tree %s != local tree %s; the push is not '
                         'content-identical' % (tree['sha'], args.verify_tree))

    message = Path(args.message_file).read_text().strip()
    commit = client.call('POST', '/git/commits',
                         {'message': message, 'tree': tree['sha'],
                          'parents': [remote_commit]})
    print('new commit %s' % commit['sha'])
    client.call('PATCH', '/git/refs/heads/' + args.branch,
                {'sha': commit['sha'], 'force': False})
    print('branch %s moved to %s' % (args.branch, commit['sha']))
    print('verified: remote tree == local tree == %s' % args.verify_tree)


if __name__ == '__main__':
    sys.exit(main())
