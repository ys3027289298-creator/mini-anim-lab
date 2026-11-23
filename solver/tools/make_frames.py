#!/usr/bin/env python3
"""Turn the solver's PPM frames into review-friendly PNG montages and GIFs.

Outputs under ARTIFACT_DIR/frames/:
  montage_hi.png   fixed vs legacy side by side at the high raster tier
  fixed_lo.gif     full fixed-order loop, low raster tier
  legacy_lo.gif    full legacy-order loop, low raster tier

Usage: python3 -I make_frames.py ARTIFACT_DIR
"""
import site
import sys
from pathlib import Path

# The Makefile runs this interpreter with -I (isolated: no cwd, no user site),
# because the frames it reads are data. Pillow only lives in the user site, so
# add that one directory back explicitly rather than dropping -I.
sys.path.append(site.getusersitepackages())

from PIL import Image, ImageDraw  # noqa: E402  (import after path fix)

KEY_FRAMES = (0, 10, 20, 30, 39)


def load_ppm(path):
    return Image.open(path).convert("RGB")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1])
    frames = root / "frames"

    # --- high-tier montage: rows = policy, columns = key frames ----------
    samples = []
    for f in KEY_FRAMES:
        samples.append(load_ppm(frames / "hi" / f"fixed_f{f:02d}.ppm"))
        samples.append(load_ppm(frames / "hi" / f"legacy_f{f:02d}.ppm"))
    fw, fh = samples[0].size
    label_h, pad = 26, 8
    montage = Image.new(
        "RGB",
        (5 * fw + 6 * pad, 2 * (fh + label_h) + 3 * pad),
        (245, 245, 245),
    )
    draw = ImageDraw.Draw(montage)
    for row, policy in enumerate(("fixed", "legacy")):
        y0 = pad + row * (fh + label_h + pad)
        draw.text((pad, y0), f"{policy} order", fill=(0, 0, 0))
        for col, f in enumerate(KEY_FRAMES):
            img = load_ppm(frames / "hi" / f"{policy}_f{f:02d}.ppm")
            x0 = pad + col * (fw + pad)
            montage.paste(img, (x0, y0 + label_h))
            draw.text((x0 + 6, y0 + 4), f"frame {f}", fill=(0, 0, 0))
    montage.save(frames / "montage_hi.png")

    # --- low-tier animated loops -----------------------------------------
    for policy in ("fixed", "legacy"):
        imgs = [
            load_ppm(frames / "lo" / f"{policy}_f{f:02d}.ppm")
            for f in range(40)
        ]
        imgs[0].save(
            frames / f"{policy}_lo.gif",
            save_all=True,
            append_images=imgs[1:],
            duration=50,
            loop=0,
        )
    print(f"wrote {frames / 'montage_hi.png'}")
    print(f"wrote {frames / 'fixed_lo.gif'}")
    print(f"wrote {frames / 'legacy_lo.gif'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
