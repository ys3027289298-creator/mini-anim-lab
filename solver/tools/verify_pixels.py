#!/usr/bin/env python3
"""Independent cross-check of the C solver's clipping-pixel statistics.

Three independent counts are compared per frame and tier:
  1. hand/arm clipping pixels recomputed from pose_*.csv by this script's
     own rasterizer (no code shared with solver.c);
  2. the penetration_*.csv numbers emitted by the C program;
  3. red pixels counted directly inside the rendered PPM frames (arm only,
     since hand and arm clipping share one red in the visualization).

Usage: python3 -I verify_pixels.py ARTIFACT_DIR
Exits 0 only if every count matches and fixed-order counts are all zero.
"""
import csv
import math
import sys
from pathlib import Path

BODY_X, BODY_Y, BODY_R = 0.6, 1.5, 0.9
R_HAND, R_ELBOW = 0.12, 0.10
VIEW_X0, VIEW_Y0, VIEW_W, VIEW_H = -1.5, -0.5, 4.5, 4.0
SCALES = {"lo": 40, "hi": 120}


def read_poses(path):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append({
                "frame": int(row["frame"]),
                "elbow": (float(row["elbow_x"]), float(row["elbow_y"])),
                "hand": (float(row["hand_x"]), float(row["hand_y"])),
            })
    return rows


def count_pixels(poses, scale):
    """Re-rasterize from poses; returns {frame: (hand_pen, arm_pen)}."""
    w, h = int(VIEW_W * scale), int(VIEW_H * scale)
    out = {}
    for pose in poses:
        hx, hy = pose["hand"]
        ex, ey = pose["elbow"]
        hand_pen = arm_pen = 0
        for py in range(h):
            wy = VIEW_Y0 + (py + 0.5) / scale
            for px in range(w):
                wx = VIEW_X0 + (px + 0.5) / scale
                body = math.hypot(wx - BODY_X, wy - BODY_Y) <= BODY_R
                if not body:
                    continue
                hand = math.hypot(wx - hx, wy - hy) <= R_HAND
                elbow = math.hypot(wx - ex, wy - ey) <= R_ELBOW
                if hand:
                    hand_pen += 1
                if hand or elbow:
                    arm_pen += 1
        out[pose["frame"]] = (hand_pen, arm_pen)
    return out


def read_pen_csv(path):
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            f = int(row["frame"])
            out[f] = {
                "lo": (int(row["hand_pen_lo"]), int(row["arm_pen_lo"])),
                "hi": (int(row["hand_pen_hi"]), int(row["arm_pen_hi"])),
            }
    return out


def read_ppm_red(path):
    """Count exact red pixels (230,0,0) in a binary P6 PPM."""
    data = Path(path).read_bytes()
    assert data[:2] == b"P6", f"not a P6 PPM: {path}"
    idx, tokens = 2, []
    while len(tokens) < 3:
        j = idx
        while data[j] in b" \t\r\n":
            j += 1
        idx = j
        while data[idx] not in b" \t\r\n":
            idx += 1
        tokens.append(int(data[j:idx]))
    w, h, _maxv = tokens
    idx += 1  # single whitespace after maxval
    rgb = data[idx:idx + w * h * 3]
    assert len(rgb) == w * h * 3, f"truncated PPM: {path}"
    red = 0
    for i in range(0, len(rgb), 3):
        if rgb[i] == 230 and rgb[i + 1] == 0 and rgb[i + 2] == 0:
            red += 1
    return red


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1])
    failures = []

    for policy in ("fixed", "legacy"):
        poses = read_poses(root / f"pose_{policy}.csv")
        pen_csv = read_pen_csv(root / f"penetration_{policy}.csv")
        for tier, scale in SCALES.items():
            mine = count_pixels(poses, scale)
            ppm_dir = root / "frames" / tier
            ppm_frames = {
                int(p.stem.rsplit("_f", 1)[1]): p
                for p in ppm_dir.glob(f"{policy}_f*.ppm")
            }
            for pose in poses:
                f = pose["frame"]
                hand_mine, arm_mine = mine[f]
                hand_csv, arm_csv = pen_csv[f][tier]
                if (hand_mine, arm_mine) != (hand_csv, arm_csv):
                    failures.append(
                        f"{policy} f{f} {tier}: python={hand_mine},{arm_mine} "
                        f"csv={hand_csv},{arm_csv}")
            # PPMs cover all frames at the low tier, key frames at the high tier
            for f, ppm_path in sorted(ppm_frames.items()):
                arm_csv = pen_csv[f][tier][1]
                red_ppm = read_ppm_red(ppm_path)
                if red_ppm != arm_csv:
                    failures.append(
                        f"{policy} f{f} {tier}: ppm_red={red_ppm} csv_arm={arm_csv}")

    fixed_pen = read_pen_csv(root / "penetration_fixed.csv")
    fixed_total = sum(v[t][i] for v in fixed_pen.values()
                      for t in ("lo", "hi") for i in (0, 1))
    legacy_pen = read_pen_csv(root / "penetration_legacy.csv")
    legacy_total = sum(v[t][i] for v in legacy_pen.values()
                       for t in ("lo", "hi") for i in (0, 1))

    print(f"frames per policy : {len(fixed_pen)}")
    print(f"fixed  total clipping pixels (both tiers, hand+arm): {fixed_total}")
    print(f"legacy total clipping pixels (both tiers, hand+arm): {legacy_total}")
    for tier in ("lo", "hi"):
        fh = sum(v[tier][0] for v in fixed_pen.values())
        fa = sum(v[tier][1] for v in fixed_pen.values())
        lh = sum(v[tier][0] for v in legacy_pen.values())
        la = sum(v[tier][1] for v in legacy_pen.values())
        print(f"  {tier:>2}: fixed hand={fh:6d} arm={fa:6d} | "
              f"legacy hand={lh:6d} arm={la:6d}")

    if failures:
        print("MISMATCH:")
        for msg in failures[:20]:
            print("  " + msg)
        return 1
    assert fixed_total == 0, "fixed order must have zero clipping pixels"
    assert legacy_total > 0, "legacy order must reproduce clipping"
    print("verify=PASS (python re-raster == C csv == PPM red counts; "
          "fixed clipping is zero)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
