"""Wrap an OASIS4.3 tree in a double-clickable macOS ``.app``.

The shipped launcher is ``OASIS4.3.command``: Finder shows it as a Terminal
document and double-clicking it leaves a Terminal window behind.  This writes a
minimal ``OASIS4.3.app`` next to it whose executable runs that same script with
no Terminal, and generates a ``.icns`` from the drawing code below.

Nothing is installed and no runtime is copied: the app only remembers where the
tree lives and execs ``<tree>/OASIS4.3.command``, so the icon and the script can
never drift apart.
"""
import argparse
import os
from pathlib import Path
import plistlib
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
NAME = 'OASIS4.3'
ICON_SIZES = [16, 32, 64, 128, 256, 512, 1024]
FONTS = ['/System/Library/Fonts/Supplemental/Arial Bold.ttf',
         '/System/Library/Fonts/Supplemental/Arial.ttf',
         '/System/Library/Fonts/Helvetica.ttc',
         '/Library/Fonts/Arial Bold.ttf']

EXECUTABLE = r"""#!/bin/bash
# Run this tree's OASIS4.3.command without a Terminal window.
set -u
here=$(cd -- "$(dirname -- "$0")" && pwd)
log="$HOME/Library/Logs/OASIS4.3-app.log"
mkdir -p "$(dirname "$log")"
tree=""
for candidate in "${OASIS43_TREE:-}" __BAKED__; do
  [ -n "$candidate" ] && [ -f "$candidate/OASIS4.3.command" ] && tree="$candidate" && break
done
if [ -z "$tree" ] && [ -f "$here/../../../OASIS4.3.command" ]; then
  tree=$(cd -- "$here/../../.." && pwd)
fi
if [ -z "$tree" ]; then
  osascript -e 'display alert "OASIS4.3 启动失败" message "找不到 OASIS4.3.command。请把 OASIS4.3.app 放回 OASIS4.3-Build 目录，或重新运行 packaging/make_macos_app.py。'
  exit 1
fi
cd -- "$tree" || exit 1
exec /bin/bash "$tree/OASIS4.3.command" "$@" </dev/null 2>>"$log"
"""


def _font(size):
    from PIL import ImageFont
    for path in FONTS:
        if Path(path).is_file():
            try:
                return ImageFont.truetype(path, size)
            except OSError:
                continue
    return ImageFont.load_default()


def draw_icon(path, size=1024):
    """Draw the app icon: a navy plate with the name and the version."""
    from PIL import Image, ImageDraw
    image = Image.new('RGBA', (size, size), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    for y in range(size):
        t = y / (size - 1)
        draw.line([(0, y), (size, y)], fill=(int(9 + 17 * t), int(35 + 25 * t), int(63 + 33 * t), 255))
    mask = Image.new('L', (size, size), 0)
    ImageDraw.Draw(mask).rounded_rectangle([0, 0, size - 1, size - 1], radius=int(size * 0.22), fill=255)
    image.putalpha(mask)
    # A thin diffraction ring behind the lettering reads as "phasing" at 32 px.
    ring = ImageDraw.Draw(image)
    inset = int(size * 0.13)
    ring.ellipse([inset, inset, size - inset, size - inset], outline=(90, 200, 235, 255),
                 width=max(2, int(size * 0.012)))
    title = _font(int(size * 0.17))
    version = _font(int(size * 0.12))
    box = draw.textbbox((0, 0), 'OASIS', font=title)
    draw.text(((size - box[2] + box[0]) / 2, size * 0.32), 'OASIS', font=title, fill=(255, 255, 255, 255))
    box = draw.textbbox((0, 0), '4.3', font=version)
    draw.text(((size - box[2] + box[0]) / 2, size * 0.58), '4.3', font=version, fill=(120, 214, 244, 255))
    image.save(path)
    return path


def make_icns(resources, size=1024):
    """Render the icon at every size macOS asks for and compile a .icns."""
    if not shutil.which('iconutil'):
        return None
    master = resources / 'icon_master.png'
    try:
        draw_icon(master, size)
    except ImportError:
        return None
    iconset = resources / 'OASIS4.3.iconset'
    if iconset.exists():
        shutil.rmtree(iconset)
    iconset.mkdir()
    from PIL import Image
    source = Image.open(master)
    for pixels in ICON_SIZES:
        for scale, suffix in ((1, ''), (2, '@2x')):
            edge = pixels * scale
            if edge > size:
                continue
            name = 'icon_%dx%d%s.png' % (pixels, pixels, suffix)
            source.resize((edge, edge), Image.LANCZOS).save(iconset / name)
    target = resources / 'OASIS4.3.icns'
    done = subprocess.run(['iconutil', '-c', 'icns', str(iconset), '-o', str(target)],
                          capture_output=True, text=True)
    shutil.rmtree(iconset)
    master.unlink()
    if done.returncode != 0:
        raise RuntimeError('iconutil 失败：' + (done.stderr or '').strip())
    return target


def build(tree=ROOT, target=None, icon=True):
    """Write ``<tree>/OASIS4.3.app``; return the bundle path."""
    tree = Path(tree).resolve()
    script = tree / 'OASIS4.3.command'
    if not script.is_file():
        raise FileNotFoundError('缺少启动脚本：%s' % script)
    target = Path(target) if target else tree / (NAME + '.app')
    target = target.resolve()
    if target.exists():
        shutil.rmtree(target)
    macos = target / 'Contents/MacOS'
    resources = target / 'Contents/Resources'
    macos.mkdir(parents=True)
    resources.mkdir(parents=True)
    baked = '"%s"' % str(tree).replace('"', '\\"')
    executable = target / 'Contents/MacOS' / NAME
    executable.write_text(EXECUTABLE.replace('__BAKED__', baked))
    executable.chmod(0o755)
    (target / 'Contents/PkgInfo').write_text('APPL????')
    icns = make_icns(resources) if icon else None
    info = dict(CFBundleName=NAME, CFBundleDisplayName=NAME, CFBundleExecutable=NAME,
                CFBundleIdentifier='cn.ac.iphy.oasis43', CFBundlePackageType='APPL',
                CFBundleVersion='4.3', CFBundleShortVersionString='4.3',
                CFBundleInfoDictionaryVersion='6.0', NSHighResolutionCapable=True,
                LSMinimumSystemVersion='11.0', NSHumanReadableCopyright='Institute of Physics, CAS')
    if icns:
        info['CFBundleIconFile'] = 'OASIS4.3'
    (target / 'Contents/Info.plist').write_bytes(plistlib.dumps(info))
    return target


def main(argv=None):
    parser = argparse.ArgumentParser(description='生成可双击的 OASIS4.3.app')
    parser.add_argument('--tree', default=str(ROOT), help='含 OASIS4.3.command 的目录')
    parser.add_argument('--target', help='生成的 .app 路径（默认 <tree>/OASIS4.3.app）')
    parser.add_argument('--no-icon', action='store_true', help='不生成图标')
    args = parser.parse_args(argv)
    app = build(args.tree, args.target, icon=not args.no_icon)
    print(app)
    return 0


if __name__ == '__main__':
    sys.exit(main())
