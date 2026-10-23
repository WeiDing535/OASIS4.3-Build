# SPDX-License-Identifier: LGPL-2.1-or-later
# Copyright (C) 2026 DingWei (丁玮) and contributors.
"""Default reference provisioned from upstream data, never from a target model.

The public repository deliberately excludes CCP4-processed reference MTZ data.
Only explicitly requested provisioning performs network access; calculations
and ordinary module imports remain offline.
"""
from pathlib import Path
import hashlib
import json
import os
import shutil
import tempfile
import urllib.request
import numpy as np
import gemmi

RESOURCE_DIR=Path(__file__).resolve().parent/'resources'
F_LABEL='FC_BEST.F_phi.F'
PHI_LABEL='FC_BEST.F_phi.phi'


def reference_info():
    return json.loads((RESOURCE_DIR/'reference-1tqw.json').read_text(encoding='utf-8'))


def reference_cache():
    override=os.environ.get('OASIS_REFERENCE_CACHE')
    return Path(override).expanduser() if override else Path.home()/'.cache'/'oasis-sad'/'reference'


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def bundled_reference():
    """Compatibility name: resolve validated packaged data or the default cache."""
    for path in [RESOURCE_DIR/'reference-1tqw.mtz',reference_cache()/'reference-1tqw.mtz']:
        if path.is_file():
            if digest(path)!=reference_info()['mtz_sha256']:
                raise ValueError('Default reference checksum mismatch; run oasis-py reference-fetch --force')
            return path
    raise ValueError('Default reference is not installed. Run oasis-py reference-fetch (internet needed once), or supply --reference-mtz. Offline setup: oasis-py reference-fetch --source /path/to/reference-1tqw.mtz')


def _extract_reference(source,destination):
    m=gemmi.read_mtz_file(str(source))
    values=np.column_stack([m.make_miller_array(),m.column_with_label(F_LABEL).array,m.column_with_label(PHI_LABEL).array]).astype('float32')
    out=gemmi.Mtz(with_base=True);out.spacegroup=m.spacegroup;out.cell=m.cell
    out.add_dataset('generic_1tqw');out.title='OASIS generic reference: original CCP4 1TQW map coefficients, unchanged'
    out.add_column(F_LABEL,'F');out.add_column(PHI_LABEL,'P')
    out.set_data(values);out.update_reso();out.write_to_file(str(destination))


def fetch_reference(source=None,force=False):
    """Provision exact pinned reference data, atomically; no software is downloaded."""
    target=reference_cache()/'reference-1tqw.mtz';info=reference_info()
    if target.is_file() and not force:
        if digest(target)==info['mtz_sha256']:return target
        raise ValueError('Cached reference is corrupt; use reference-fetch --force to replace it')
    target.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.reference-',dir=target.parent) as folder:
        original=Path(folder)/'upstream.mtz';prepared=Path(folder)/'reference-1tqw.mtz'
        if source is not None:
            source=Path(source).expanduser().resolve()
            if not source.is_file():raise ValueError(f'Reference source not found: {source}')
            if source.stat().st_size>32*1024*1024:raise ValueError('Reference source exceeds the expected size limit')
            shutil.copyfile(source,original)
        else:
            url=info['public_data_url']
            if not url.startswith('https://'):raise ValueError('Reference download requires HTTPS')
            request=urllib.request.Request(url,headers={'User-Agent':'OASIS4.3-GUI-0.1 reference provisioner'})
            total=0
            with urllib.request.urlopen(request,timeout=60) as response,original.open('wb') as output:
                if not response.url.startswith('https://'):raise ValueError('Reference redirect is not HTTPS')
                for block in iter(lambda:response.read(1024*1024),b''):
                    total+=len(block)
                    if total>32*1024*1024:raise ValueError('Reference download exceeds the expected size limit')
                    output.write(block)
        actual=digest(original)
        if actual==info['mtz_sha256']:shutil.copyfile(original,prepared)
        elif actual==info['source_sha256']:_extract_reference(original,prepared)
        else:raise ValueError('Upstream reference SHA256 mismatch; cache was not changed')
        if digest(prepared)!=info['mtz_sha256']:
            raise ValueError('Extracted reference checksum mismatch; check the documented Gemmi version')
        read_reference(prepared)
        os.replace(prepared,target)
    return target


def read_reference(path=None):
    path=Path(path) if path is not None else bundled_reference()
    m=gemmi.read_mtz_file(str(path))
    if m.spacegroup is None or not m.cell.is_crystal():raise ValueError('Reference MTZ needs a valid crystal cell and space group')
    columns={c.label:c.type for c in m.columns}
    if columns.get(F_LABEL)!='F' or columns.get(PHI_LABEL)!='P':
        raise ValueError(f'Reference MTZ requires amplitude {F_LABEL} (F) and phase {PHI_LABEL} (P)')
    m.ensure_asu();f=np.asarray(m.column_with_label(F_LABEL).array);phi=np.asarray(m.column_with_label(PHI_LABEL).array)
    if np.count_nonzero(np.isfinite(f)&(f>0)&np.isfinite(phi))<100:
        raise ValueError('Reference MTZ has fewer than 100 finite positive map coefficients')
    return m
