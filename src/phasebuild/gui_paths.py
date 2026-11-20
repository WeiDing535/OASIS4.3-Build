"""Resolve paths pasted from Finder or a file URL, with explicit navigation errors."""
from pathlib import Path
from urllib.parse import unquote,urlparse


def browse_path(raw, allow_parent=False):
    text=str(raw).strip()
    if len(text)>=2 and text[0]==text[-1] and text[0] in ('"',"'"):
        text=text[1:-1]
    if text.startswith('file:'):
        parsed=urlparse(text)
        if parsed.netloc not in ('','localhost'):
            raise ValueError('only local file paths are supported')
        text=unquote(parsed.path)
    if not text:
        raise ValueError('please provide the full path of a file or directory')
    path=Path(text).expanduser().resolve()
    requested=str(path)
    if not path.exists() and not allow_parent:
        raise ValueError('path does not exist, please check the spelling: '+str(path))
    while not path.exists() and path!=path.parent:
        path=path.parent
    selected=str(path) if path.is_file() else None
    directory=path.parent if selected else path
    if not directory.is_dir():
        raise ValueError('not a browsable directory')
    return directory,selected,requested!=str(path)
