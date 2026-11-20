"""Locate and validate the user-supplied CCP4 installation.

This module is the single source of truth for external scientific engines.
No third-party binaries or data are shipped with the source distribution; users
point OASIS at their own licensed CCP4 installation. Density modification uses
CCP4 Parrot and the trace/Fcalc feedback uses REFMAC zero-cycle structure factors
from the same CCP4 installation.
"""
from __future__ import annotations
import os
from pathlib import Path

# Binaries provided by a standard CCP4 installation.
CCP4_PROGRAMS = ('phaser', 'cparrot', 'cbuccaneer', 'refmac5')
# Data files a standard CCP4 installation provides (generic unrelated reference
# histogram used only for density-scaling histograms, never a target structure).
CCP4_DATA = ('lib/data/reference_structures/reference-1tqw.mtz',
             'lib/data/reference_structures/reference-1tqw.pdb',
             'lib/data/monomers')


def _fmt_missing(items, label):
    return label + ' missing required item: ' + ', '.join(items) + '. Please confirm this points to a complete installation directory.'


def ccp4_candidates(value=None):
    """Yield candidate CCP4 roots in priority order.

    An explicit path is authoritative: if it is incomplete we report the exact
    missing items instead of silently falling back to another installation.
    This is what makes user-pinned engine paths reproducible and auditable.
    """
    if value:
        yield Path(value).expanduser()
        return
    for env in ('PHASEBUILD_CCP4', 'CCP4'):
        if os.environ.get(env):
            yield Path(os.environ[env]).expanduser()
    for base in ('/Applications/ccp4-9', '/Applications/ccp4', '/opt/ccp4',
                 '/usr/local/ccp4'):
        yield Path(base)


def ccp4_directory(value=None, required=CCP4_PROGRAMS, *, data=True, inspect=True):
    """Resolve and validate a CCP4 installation root.

    ``value`` is an explicit path (CLI ``--ccp4-root`` / env). Returns the
    resolved root as str. Raises ValueError with a clear Chinese message when
    the location lacks required binaries or data. ``inspect=False`` performs
    only path resolution (used for previews).
    """
    tried = []
    for raw in ccp4_candidates(value):
        root = raw.resolve()
        if root in tried:
            continue
        tried.append(root)
        if not root.is_dir():
            continue
        if not inspect:
            return str(root)
        missing = [name for name in required
                   if not (root / 'bin' / name).is_file()
                   or not os.access(root / 'bin' / name, os.X_OK)]
        if data:
            missing += [rel for rel in CCP4_DATA if not (root / rel).exists()]
        if not missing:
            return str(root)
        # Keep scanning other candidates but remember why this one failed.
        tried[-1] = root
    hint = value or os.environ.get('PHASEBUILD_CCP4') or os.environ.get('CCP4')
    if hint:
        root = Path(hint).expanduser().resolve()
        missing = [name for name in required
                   if not (root / 'bin' / name).is_file()
                   or not os.access(root / 'bin' / name, os.X_OK)]
        if data:
            missing += [rel for rel in CCP4_DATA if not (root / rel).exists()]
        raise ValueError(_fmt_missing(missing, 'CCP4 directory ' + str(root)))
    raise ValueError('could not find a complete CCP4 installation directory; please use --ccp4-root or the environment variable '
                     'PHASEBUILD_CCP4 to specify the path (the directory must contain bin/ and lib/data/). '
                     'tried: ' + ', '.join(str(p) for p in tried[:6]))


def detect_ccp4(value=None):
    """Best-effort detect a complete CCP4 installation for UI auto-fill.

    Returns ``(root, message)`` and never raises. ``root`` is empty when no
    complete installation was found; the message explains what was checked.
    """
    try:
        root = ccp4_directory(value, inspect=True)
        return root, 'detected a complete CCP4 installation: ' + root
    except ValueError as error:
        return '', str(error)
