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
            raise ValueError('仅支持本机 file 路径')
        text=unquote(parsed.path)
    if not text:
        raise ValueError('请输入文件或目录的完整路径')
    path=Path(text).expanduser().resolve()
    requested=str(path)
    if not path.exists() and not allow_parent:
        raise ValueError('路径不存在，请检查拼写：'+str(path))
    while not path.exists() and path!=path.parent:
        path=path.parent
    selected=str(path) if path.is_file() else None
    directory=path.parent if selected else path
    if not directory.is_dir():
        raise ValueError('该路径不是可浏览的目录')
    return directory,selected,requested!=str(path)
