#!/usr/bin/env python3
"""Unified self-check runner for the five mini-anim-lab modules.

Each module is built and run five ways:

  normal   the module's own self-check (every check it already had)
  pose     the same source plus -DMAL_EVIDENCE_POSE, printing mal_pose lines
  time     two builds, -O0 and -O2, run the same workload for milliseconds
  image    -DMAL_EVIDENCE_IMAGE, rasterizes a pose and prints a raster CRC
  replay   rebuild+rerun command carrying the source checksum

Failures are reported per module with a non-zero exit code, and the five
pieces of evidence are written to artifacts/selfcheck.json (plus a readable
text summary and the raster PNGs).
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

MODULES = ["skeleton", "ik", "curve", "constraint", "path"]

# Every check the module had before the unified runner existed. The runner
# fails a module if any of these names disappear or flips to 0.
BASELINE = {
    "skeleton": [
        "add_root", "add_child", "world_chain", "rotate_chain",
        "segment_dir", "reject_bad_parent", "capacity",
    ],
    "ik": [
        "reachable", "joint_on_circle", "end_matches", "pole_up", "pole_down",
        "unreachable", "too_close", "zero_len",
    ],
    "curve": [
        "add_first", "add_second", "clamp_left", "clamp_right",
        "mid_monotonic", "endpoints_exact", "flat_curve", "flat_deriv_zero",
        "reject_unsorted", "capacity",
    ],
    "constraint": [
        "init_len", "stays_in_stretch", "pinned_anchor", "falls_down",
        "more_iters_better", "capacity_ok",
    ],
    "path": [
        "add0", "add1", "add2", "total_len", "sample_mid", "sample_corner",
        "sample_end", "sample_out", "dup_rejected", "closest_dist",
        "closest_s", "curvature_straight", "curvature_turn",
    ],
}

CHECK_RE = re.compile(r"^([A-Za-z0-9_]+)=([01])$")
POSE_RE = re.compile(r"^mal_pose k=(\S+) v=(\S+)$")
POSE_HASH_RE = re.compile(r"^mal_pose_hash v=([0-9a-f]+) n=(\d+)$")
TIME_RE = re.compile(r"^mal_time tier=(\S+) ms=([0-9.]+) iterations=(\d+)$")
PNG_RE = re.compile(r"^mal_png label=(\S+) w=(\d+) h=(\d+) crc=([0-9a-f]+)$")


class ModuleError(Exception):
    pass


def run(cmd: list[str], cwd: Path) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_output(text: str) -> dict:
    checks: dict[str, str] = {}
    order: list[str] = []
    pose: dict[str, float] = {}
    pose_hash = None
    pose_n = None
    times: dict[str, dict] = {}
    png = None
    all_ok = None
    for line in text.splitlines():
        if line.startswith("all_ok="):
            all_ok = line.split("=", 1)[1] == "1"
            continue
        m = CHECK_RE.match(line)
        if m:
            name, value = m.group(1), m.group(2)
            if name not in checks:
                order.append(name)
            checks[name] = value
            continue
        m = POSE_RE.match(line)
        if m:
            pose[m.group(1)] = float(m.group(2))
            continue
        m = POSE_HASH_RE.match(line)
        if m:
            pose_hash = m.group(1)
            pose_n = int(m.group(2))
            continue
        m = TIME_RE.match(line)
        if m:
            times[m.group(1)] = {
                "ms": float(m.group(2)),
                "iterations": int(m.group(3)),
            }
            continue
        m = PNG_RE.match(line)
        if m:
            png = {
                "label": m.group(1),
                "w": int(m.group(2)),
                "h": int(m.group(3)),
                "crc": m.group(4),
            }
            continue
    return {
        "checks": checks,
        "order": order,
        "pose": pose,
        "pose_hash": pose_hash,
        "pose_n": pose_n,
        "times": times,
        "png": png,
        "all_ok": all_ok,
    }


def compile_module(root: Path, module: str, cc: str, opt: str, defs: list[str],
                   out: Path) -> Path:
    src = root / module / "src" / f"{module}.c"
    include = [str(root / module / "src"), str(root / "common")]
    cmd = [cc, opt, "-std=c11", "-Wall", "-Wextra"]
    for d in defs:
        cmd.append(f"-D{d}")
    for inc in include:
        cmd.extend(["-I", inc])
    cmd.extend([str(src), "-o", str(out), "-lm"])
    res = run(cmd, root)
    if res.returncode != 0:
        raise ModuleError(f"compile failed ({opt}): {res.stderr.strip()[:400]}")
    return out


def load_baseline(module: str) -> list[str]:
    return BASELINE[module]


def validate_checks(module: str, parsed: dict) -> list[str]:
    errors = []
    expected = load_baseline(module)
    got = parsed["checks"]
    for name in expected:
        if name not in got:
            errors.append(f"missing check {name!r}")
        elif got[name] != "1":
            errors.append(f"check {name!r} failed")
    extra = [n for n in parsed["order"] if n not in expected]
    if extra:
        errors.append(f"unexpected checks {extra!r}")
    if len(parsed["order"]) != len(set(parsed["order"])):
        errors.append("duplicate check names")
    if not parsed["all_ok"]:
        errors.append("all_ok is not 1")
    return errors


def check_module(root: Path, module: str, cc: str, build_dir: Path) -> dict:
    start = time.monotonic()
    src = root / module / "src" / f"{module}.c"
    sha = sha256_file(src)
    base = build_dir / module
    if base.exists():
        shutil.rmtree(base)
    base.mkdir(parents=True)

    record: dict = {
        "module": module,
        "source": str(src.relative_to(root)),
        "sha256": sha,
        "checks": {},
        "pose": None,
        "timing_ms": {"slow": None, "fast": None},
        "image_crc": None,
        "image_size": None,
        "replay": None,
    }
    errors: list[str] = []
    transcript: dict[str, dict] = {}

    try:
        # 1. normal build: the module's own self-check
        exe = compile_module(root, module, cc, "-O2", [], base / f"{module}-normal")
        res = run([str(exe)], root)
        parsed = parse_output(res.stdout)
        transcript["normal"] = parsed
        if res.returncode != 0:
            errors.append(f"normal self-check exited {res.returncode}")
        errors.extend(validate_checks(module, parsed))
        record["checks"] = dict(parsed["checks"])
        if errors:
            raise ModuleError("; ".join(errors))

        # 2. pose evidence
        exe = compile_module(root, module, cc, "-O2", ["MAL_EVIDENCE_POSE"],
                             base / f"{module}-pose")
        res = run([str(exe)], root)
        parsed = parse_output(res.stdout)
        transcript["pose"] = parsed
        if res.returncode != 0:
            raise ModuleError(f"pose build exited {res.returncode}")
        errors.extend(validate_checks(module, parsed))
        if not parsed["pose"] or not parsed["pose_hash"]:
            raise ModuleError("pose evidence missing")
        record["pose"] = {
            "hash": parsed["pose_hash"],
            "count": parsed["pose_n"],
            "values": parsed["pose"],
        }
        if errors:
            raise ModuleError("; ".join(errors))

        # 3. timing: two optimisation tiers, same source and iteration count.
        # Each tier runs TIMING_REPS times; the median absorbs scheduler noise.
        TIMING_REPS = 3
        timings = {}
        iters = None
        for tier, opt in (("slow", "-O0"), ("fast", "-O2")):
            exe = compile_module(
                root, module, cc, opt,
                ["MAL_EVIDENCE_TIMING", f"MAL_TIME_TIER=\"{tier}\""],
                base / f"{module}-{tier}",
            )
            samples = []
            for _ in range(TIMING_REPS):
                res = run([str(exe)], root)
                parsed = parse_output(res.stdout)
                transcript[f"{tier}-{len(samples)}"] = parsed
                if res.returncode != 0:
                    raise ModuleError(f"{tier} timing run exited {res.returncode}")
                errors.extend(validate_checks(module, parsed))
                if tier not in parsed["times"]:
                    raise ModuleError(f"{tier} timing line missing")
                sample = parsed["times"][tier]
                if sample["iterations"] <= 0 or sample["ms"] <= 0:
                    raise ModuleError(f"{tier} timing is not a positive measurement")
                if iters is not None and sample["iterations"] != iters:
                    raise ModuleError("timing tiers ran different iteration counts")
                iters = sample["iterations"]
                samples.append(sample["ms"])
            if errors:
                raise ModuleError("; ".join(errors))
            samples.sort()
            timings[tier] = {
                "ms": samples[len(samples) // 2],
                "samples": samples,
            }
        record["timing_ms"] = {
            "slow": round(timings["slow"]["ms"], 3),
            "fast": round(timings["fast"]["ms"], 3),
            "iterations": iters,
            "reps": TIMING_REPS,
            "samples_ms": {"slow": [round(x, 3) for x in timings["slow"]["samples"]],
                           "fast": [round(x, 3) for x in timings["fast"]["samples"]]},
            "speedup": round(timings["slow"]["ms"] / timings["fast"]["ms"], 3),
        }

        # 4. image evidence
        exe = compile_module(root, module, cc, "-O2", ["MAL_EVIDENCE_IMAGE"],
                             base / f"{module}-image")
        res = run([str(exe)], root)
        parsed = parse_output(res.stdout)
        transcript["image"] = parsed
        if res.returncode != 0:
            raise ModuleError(f"image build exited {res.returncode}")
        errors.extend(validate_checks(module, parsed))
        if not parsed["png"]:
            raise ModuleError("image evidence missing")
        record["image_crc"] = parsed["png"]["crc"]
        record["image_size"] = [parsed["png"]["w"], parsed["png"]["h"]]

        # 5. replay: a command that rebuilds and reruns the check, which we
        # execute here so the recorded evidence is known to work verbatim
        record["replay"] = (
            f"{cc} -O2 -std=c11 -Wall -Wextra -I {module}/src -I common "
            f"{module}/src/{module}.c -o {build_dir.relative_to(root)}/{module}-replay -lm "
            f"&& ./{build_dir.relative_to(root)}/{module}-replay  # {module} sha256={sha}"
        )
        res = run(["bash", "-c", record["replay"]], root)
        parsed = parse_output(res.stdout)
        if res.returncode != 0 or not parsed["all_ok"]:
            raise ModuleError("replay command did not reproduce the self-check")
        errors.extend(validate_checks(module, parsed))
    except ModuleError as e:
        errors.append(str(e))

    # a stage that fails mid-validation can report the same problem twice
    errors = list(dict.fromkeys(errors))
    record["errors"] = errors
    record["ok"] = not errors
    record["elapsed_s"] = round(time.monotonic() - start, 3)
    record["transcript"] = transcript
    return record


def represent(value: float) -> str:
    return f"{value:.9g}"


def write_text_report(records: list[dict], path: Path) -> None:
    lines = []
    total = 0
    passed = 0
    for rec in records:
        status = "PASS" if rec["ok"] else "FAIL"
        n = len(rec["checks"])
        n_ok = sum(1 for v in rec["checks"].values() if v == "1")
        total += n
        passed += n_ok
        lines.append(f"[{status}] {rec['module']}: checks {n_ok}/{n}  hash={rec['sha256'][:12]}")
        if rec["pose"]:
            lines.append(f"    pose   hash={rec['pose']['hash']} n={rec['pose']['count']}")
        else:
            lines.append("    pose   MISSING")
        t = rec["timing_ms"]
        if t.get("fast") is not None:
            lines.append(f"    time   slow={t['slow']:.3f}ms fast={t['fast']:.3f}ms "
                         f"x{t['speedup']:.2f} ({t['iterations']} iters)")
        else:
            lines.append("    time   MISSING")
        if rec["image_crc"]:
            lines.append(f"    image  crc={rec['image_crc']}")
        else:
            lines.append("    image  MISSING")
        if rec["replay"]:
            lines.append(f"    replay {rec['replay']}")
        else:
            lines.append("    replay MISSING")
        for err in rec["errors"]:
            lines.append(f"    ERROR  {err}")
        lines.append("")
    lines.append(f"TOTAL checks {passed}/{total}")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cc", default=os.environ.get("CC", "cc"))
    ap.add_argument("--build-dir", default="artifacts/build")
    ap.add_argument("--json", default="artifacts/selfcheck.json")
    ap.add_argument("--report", default="artifacts/selfcheck.txt")
    ap.add_argument("--module", action="append", choices=MODULES,
                    help="limit to one or more modules (default: all)")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    root = Path(__file__).resolve().parent.parent
    build_dir = (root / args.build_dir).resolve()
    build_dir.mkdir(parents=True, exist_ok=True)
    modules = args.module or MODULES

    records = []
    for module in modules:
        rec = check_module(root, module, args.cc, build_dir)
        records.append(rec)
        if not args.quiet:
            status = "PASS" if rec["ok"] else "FAIL"
            n_ok = sum(1 for v in rec["checks"].values() if v == "1")
            print(f"[{status}] {module}: {n_ok}/{len(rec['checks'])} checks", flush=True)
            for err in rec["errors"]:
                print(f"         {err}", flush=True)

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "cc": args.cc,
        "modules": records,
    }
    json_path = root / args.json
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    report_path = root / args.report
    write_text_report(records, report_path)

    ok = all(rec["ok"] for rec in records)
    total_checks = sum(len(r["checks"]) for r in records)
    if not args.quiet:
        print(f"\njson:   {json_path.relative_to(root)}")
        print(f"report: {report_path.relative_to(root)}")
        print(f"result: {'PASS' if ok else 'FAIL'} ({total_checks} checks across "
              f"{len(records)} modules)")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
