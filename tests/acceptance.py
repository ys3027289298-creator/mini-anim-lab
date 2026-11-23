#!/usr/bin/env python3
"""
acceptance.py - unified acceptance driver for the new pose export format.

For each of skeleton/ik/curve/path it verifies:

  1. Old format untouched: running the binary with no arguments still prints
     exactly tests/golden/<module>.out, byte for byte.
  2. New format complete: --pose-export <dir> produces pose.json + pose.ppm,
     and pose.json contains all four required deliverables:
       pose values, image hash, two timing tiers, reproduction command.
  3. Image hash truthful: the sha256 recorded in pose.json equals the actual
     sha256 of pose.ppm (recomputed here, outside the C code).
  4. Deterministic: a second export renders a byte-identical pose.ppm and
     identical pose values.
  5. Reproducible: the exact "reproduce" command recorded in pose.json can be
     re-run from the repository root and regenerates the same image.

A module missing any of the four data items counts as incomplete and the
script exits non-zero.

Usage:
  python3 -I tests/acceptance.py                 # full acceptance run
  python3 -I tests/acceptance.py --update-golden # refresh tests/golden/*.out
"""
import argparse
import filecmp
import hashlib
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MODULES = ["skeleton", "ik", "curve", "path"]
GOLDEN_DIR = ROOT / "tests" / "golden"

# pose.json must contain these four deliverables, in this structure.
REQUIRED = {
    "module": str,
    "format_version": int,
    "pose": dict,
    "image": dict,
    "timing": dict,
    "reproduce": list,
}
REQUIRED_IMAGE_KEYS = {"file", "format", "width", "height", "sha256"}
REQUIRED_TIMING_KEYS = {"compute_pose_us", "compute_repeats", "render_image_us"}


class Failure(Exception):
    pass


def run(cmd, cwd):
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def check_old_format(module):
    """No-argument invocation must match the golden stdout byte for byte."""
    rc, out, err = run([f"./{module}/{module}"], cwd=ROOT)
    if rc != 0:
        raise Failure(f"old-format run exited {rc}: {err.strip()}")
    golden = (GOLDEN_DIR / f"{module}.out").read_bytes()
    actual = out.encode()
    if actual != golden:
        raise Failure(
            f"old output differs from golden "
            f"({len(actual)} vs {len(golden)} bytes)"
        )
    return f"stdout {len(actual)} B identical to tests/golden/{module}.out"


def load_export(module, outdir):
    rc, out, err = run(
        [f"./{module}/{module}", "--pose-export", str(outdir)], cwd=ROOT
    )
    if rc != 0:
        raise Failure(f"pose export exited {rc}: {err.strip()}")
    if out or err:
        raise Failure("export mode must stay silent on stdout/stderr")
    json_path = outdir / "pose.json"
    ppm_path = outdir / "pose.ppm"
    if not json_path.is_file() or not ppm_path.is_file():
        raise Failure("pose.json or pose.ppm missing")
    data = json.loads(json_path.read_text())
    return data, json_path, ppm_path


def validate_schema(module, data, ppm_path):
    for key, typ in REQUIRED.items():
        if key not in data:
            raise Failure(f"pose.json missing top-level key: {key}")
        if not isinstance(data[key], typ):
            raise Failure(f"pose.json key {key} has wrong type")
    if data["module"] != module:
        raise Failure(f"module field is {data['module']!r}, expected {module!r}")
    if not data["pose"]:
        raise Failure("deliverable missing: pose values (pose object empty)")
    img = data["image"]
    missing = REQUIRED_IMAGE_KEYS - img.keys()
    if missing:
        raise Failure(f"deliverable missing: image block lacks {sorted(missing)}")
    if len(img["sha256"]) != 64:
        raise Failure("image sha256 is not 64 hex chars")
    # verify the PPM container itself
    head = ppm_path.read_bytes()[:16]
    if not head.startswith(b"P6\n"):
        raise Failure("pose.ppm is not a P6 PPM")
    actual_hash = hashlib.sha256(ppm_path.read_bytes()).hexdigest()
    if actual_hash != img["sha256"]:
        raise Failure(
            f"recorded hash {img['sha256'][:12]}... != actual {actual_hash[:12]}..."
        )
    timing = data["timing"]
    missing_t = REQUIRED_TIMING_KEYS - timing.keys()
    if missing_t:
        raise Failure(f"deliverable missing: timing lacks {sorted(missing_t)}")
    # two timing tiers must both be present and non-negative integers
    for k in ("compute_pose_us", "render_image_us"):
        if not isinstance(timing[k], int) or timing[k] < 0:
            raise Failure(f"timing tier {k} must be a non-negative integer")
    if not data["reproduce"] or not all(
        isinstance(x, str) and x for x in data["reproduce"]
    ):
        raise Failure("deliverable missing: reproduction command empty")
    return actual_hash


def check_determinism(module, first_ppm, first_json):
    with tempfile.TemporaryDirectory() as td:
        second = Path(td) / module
        data2, _, ppm2 = load_export(module, second)
        if not filecmp.cmp(first_ppm, ppm2, shallow=False):
            raise Failure("pose.ppm differs between two runs (not deterministic)")
        # pose values must be identical even though timing naturally varies
        p1 = json.loads(first_json.read_text())["pose"]
        if data2["pose"] != p1:
            raise Failure("pose values differ between two runs")
        if data2["image"]["sha256"] != hashlib.sha256(ppm2.read_bytes()).hexdigest():
            raise Failure("second run hash not internally consistent")


def check_reproduce(module, data):
    with tempfile.TemporaryDirectory() as td:
        repro_dir = Path(td) / "repro"
        cmd = list(data["reproduce"])
        # the recorded command points at the original outdir; redirect the
        # re-run to a fresh dir while keeping every other argument identical
        cmd[-1] = str(repro_dir)
        rc, out, err = run(cmd, cwd=ROOT)
        if rc != 0:
            raise Failure(f"reproduce command exited {rc}: {err.strip()}")
        repro_ppm = repro_dir / "pose.ppm"
        if not repro_ppm.is_file():
            raise Failure("reproduce command produced no pose.ppm")
        h = hashlib.sha256(repro_ppm.read_bytes()).hexdigest()
        if h != data["image"]["sha256"]:
            raise Failure("reproduced image hash differs from recorded hash")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--update-golden", action="store_true")
    args = ap.parse_args()

    if args.update_golden:
        for m in MODULES:
            rc, out, err = run([f"./{m}/{m}"], cwd=ROOT)
            if rc != 0:
                print(f"[{m}] cannot update golden, run failed: {err}")
                return 1
            (GOLDEN_DIR / f"{m}.out").write_text(out)
            print(f"[{m}] golden updated ({len(out)} bytes)")
        return 0

    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    report = {}
    all_ok = True

    for m in MODULES:
        result = {"ok": False, "checks": [], "deliverables": {}}
        try:
            result["checks"].append(("old format byte-identical",
                                     check_old_format(m)))
            with tempfile.TemporaryDirectory() as td:
                outdir = Path(td) / m
                data, json_path, ppm_path = load_export(m, outdir)
                h = validate_schema(m, data, ppm_path)
                result["checks"].append(
                    ("four deliverables present + hash verified",
                     f"pose={len(json.dumps(data['pose']))} B, "
                     f"sha256={h[:16]}..."))
                check_determinism(m, ppm_path, json_path)
                result["checks"].append(
                    ("deterministic across reruns",
                     "identical pose.ppm and pose values"))
                check_reproduce(m, data)
                result["checks"].append(
                    ("reproduce command reruns cleanly",
                     " ".join(data["reproduce"])))
                result["deliverables"] = {
                    "pose_values": data["pose"],
                    "image_sha256": h,
                    "timing": data["timing"],
                    "reproduce": data["reproduce"],
                }
            result["ok"] = True
        except Failure as exc:
            result["error"] = str(exc)
            all_ok = False
        report[m] = result

    # ---------------- unified acceptance report ----------------
    print("=" * 72)
    print("mini-anim-lab pose export acceptance report")
    print("=" * 72)
    for m in MODULES:
        r = report[m]
        status = "PASS" if r["ok"] else "FAIL"
        print(f"\n[{m}] {status}")
        if not r["ok"]:
            print(f"  ERROR: {r['error']}")
            continue
        for name, detail in r["checks"]:
            print(f"  - {name}: {detail}")
        d = r["deliverables"]
        print(f"  deliverables (all four required):")
        print(f"    1) pose values   : {json.dumps(d['pose_values'])[:80]}...")
        print(f"    2) image hash    : {d['image_sha256']}")
        print(f"    3) timing (us)   : compute={d['timing']['compute_pose_us']}"
              f" (x{d['timing']['compute_repeats']}), "
              f"render={d['timing']['render_image_us']}")
        print(f"    4) reproduce     : {' '.join(d['reproduce'])}")

    print("\n" + "=" * 72)
    passed = sum(1 for m in MODULES if report[m]["ok"])
    print(f"summary: {passed}/{len(MODULES)} modules accepted "
          f"({'ALL OK' if all_ok else 'INCOMPLETE - four datasets required'})")
    print("=" * 72)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
