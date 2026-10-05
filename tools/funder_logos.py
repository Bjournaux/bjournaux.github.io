"""
Make the one-color versions of the footer logos.

    pip install pillow
    python tools/funder_logos.py

For every logo in images/funders/color/ (PNG with a transparent background),
writes images/funders/mono/ with the same name: the logo in plain white on a
transparent background. White parts of a logo (e.g. the letters in the NASA or
NSF logos) become see-through, so they stay readable.

The footer shows color/ or mono/ depending on `footer-logos` in _config.yaml.
If a mono/ file is missing, the footer turns the color/ file white with CSS.
"""
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parent.parent / "images" / "funders"
WHITE = 0.92  # pixels lighter than this count as white (see-through)
SOFT = 0.12   # soft edge between ink and white


def mono(src, dst):
    im = Image.open(src).convert("RGBA")
    r, g, b, a = (band.load() for band in im.split())
    out = Image.new("RGBA", im.size, (255, 255, 255, 0))
    px = out.load()
    for y in range(im.height):
        for x in range(im.width):
            lum = (0.299 * r[x, y] + 0.587 * g[x, y] + 0.114 * b[x, y]) / 255
            ink = min(1.0, max(0.0, (WHITE - lum) / SOFT))
            px[x, y] = (255, 255, 255, round(a[x, y] * ink))
    out.save(dst, optimize=True)


def main():
    (ROOT / "mono").mkdir(exist_ok=True)
    for src in sorted((ROOT / "color").glob("*.png")):
        mono(src, ROOT / "mono" / src.name)
        print("wrote", (ROOT / "mono" / src.name).relative_to(ROOT.parent.parent))


if __name__ == "__main__":
    main()
