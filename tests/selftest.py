#!/usr/bin/env python3
"""selftest.py - evidence-driven self-test bench for mini-anim-lab.

For every module it:
  1. builds clean in an isolated worktree copy and runs the module self-test,
  2. records pass count and sha256 of every visible artifact,
  3. injects a battery of defects and asks whether the self-test catches them,
  4. reports capture rate; anything under 80% lists the missed defects,
  5. prints exact reproduction commands.

Standard library only. Run with:  python3 -I tests/selftest.py
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODULES = ["skeleton", "ik", "constraint", "path", "curve"]
CAPTURE_FLOOR = 0.80

# A defect is (id, description, list of (old, new) exact source replacements).
# Capture is detected via: build failure / nonzero test exit / artifact hash
# drift / missing artifact. The healthy build must always pass first.
DEFECTS = {
    "skeleton": [
        ("rot_sign_x", "rotate() flips rotation x sign (mirrored skeleton)",
         [("*ox = x * c - y * s;", "*ox = x * c + y * s;")]),
        ("ang_sum_sign", "world angle uses subtraction instead of sum",
         [("s->world[i].ang = s->world[p].ang + s->local[i].ang;",
           "s->world[i].ang = s->world[p].ang - s->local[i].ang;")]),
        ("seg_wrong_axis", "bone segment points along y instead of x",
         [("rotate(1.0f, 0.0f, s->world[i].ang, &rx, &ry);",
           "rotate(0.0f, 1.0f, s->world[i].ang, &rx, &ry);")]),
        ("parent_offset_dropped", "child world position ignores parent origin",
         [("s->world[i].x = s->world[p].x + rx;", "s->world[i].x = rx;"),
          ("s->world[i].y = s->world[p].y + ry;", "s->world[i].y = ry;")]),
        ("raster_gap", "canvas line width shrunk (image-only defect, hash catch)",
         [("for (int i = 0; i <= n; i++) {", "for (int i = 0; i <= n; i += 2) {")]),
        ("csv_precision_cut", "pose CSV loses two decimal digits (artifact drift)",
         [("%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\\n",
           "%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\\n")]),
    ],
    "ik": [
        ("law_of_cosines", "cosine law numerator sign flipped (wrong joint)",
         [("float cos_a = (l1 * l1 + dist * dist - l2 * l2)",
           "float cos_a = (l1 * l1 - dist * dist - l2 * l2)")]),
        ("pole_side_flip", "pole vector bend direction inverted",
         [("float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? 1.0f : -1.0f;",
           "float side = (dir.x * pole.y - dir.y * pole.x) >= 0 ? -1.0f : 1.0f;")]),
        ("reach_limit_loose", "unreachable targets accepted (guard widened)",
         [("if (dist > l1 + l2 + 1e-5f) return r;",
           "if (dist > l1 + l2 + 5.0f) return r;")]),
        ("zero_len_guard_gone", "non-positive bone length guard removed",
         [("if (l1 <= 0 || l2 <= 0) return r;", "if (0) return r;")]),
        ("raster_gap", "rasterized pose skips pixels (image-only, hash catch)",
         [("for (int i = 0; i <= n; i++) {", "for (int i = 0; i <= n; i += 2) {")]),
        ("csv_angle_dropped", "pose CSV stops writing joint angles",
         [("fprintf(f, \"up,%d,%.6f,%.6f,%.6f,%.6f\\n\", r.ok, r.joint.x, r.joint.y, r.a1, r.a2);",
           "fprintf(f, \"up,%d,%.6f,%.6f,0.0,0.0\\n\", r.ok, r.joint.x, r.joint.y);")]),
    ],
    "constraint": [
        ("verlet_gravity_doubled", "gravity term doubled in integration",
         [("p->x += vx + gx * dt * dt;", "p->x += vx + 2.0f * gx * dt * dt;"),
          ("p->y += vy + gy * dt * dt;", "p->y += vy + 2.0f * gy * dt * dt;")]),
        ("constraint_relaxed", "distance constraint only half-satisfied",
         [("float diff = (d - c->rest[i]) / d;",
           "float diff = 0.5f * (d - c->rest[i]) / d;")]),
        ("pin_ignored", "pinned particles integrate like free ones",
         [("if (p->pinned) continue;", "if (0) continue;")]),
        ("damping_killed", "velocity damping set to 1 (no energy loss)",
         [("c->damping = 0.99f;", "c->damping = 1.0f;")]),
        ("raster_gap", "chain raster skips pixels (image-only, hash catch)",
         [("for (int i = 0; i <= n; i++) {", "for (int i = 0; i <= n; i += 2) {")]),
        ("csv_error_zeroed", "reported segment error always written as 0",
         [("i, c->p[i].x, c->p[i].y, c->p[i].pinned, c->p[i].inv_mass, worst);",
           "i, c->p[i].x, c->p[i].y, c->p[i].pinned, c->p[i].inv_mass, 0.0f);")]),
    ],
    "path": [
        ("arclength_overcount", "cumulative arc length adds each segment twice",
         [("pa->cum[pa->n] = pa->cum[pa->n - 1] + d;",
           "pa->cum[pa->n] = pa->cum[pa->n - 1] + 2.0f * d;")]),
        ("sample_lerp_broken", "constant-speed sample ignores interpolation t",
         [("out->x = pa->p[i - 1].x + (pa->p[i].x - pa->p[i - 1].x) * t;",
           "out->x = pa->p[i - 1].x;"),
          ("out->y = pa->p[i - 1].y + (pa->p[i].y - pa->p[i - 1].y) * t;",
           "out->y = pa->p[i - 1].y;")]),
        ("closest_unclamped", "closest-point t not clamped to [0,1]",
         [("if (t < 0.0f) t = 0.0f;\n        if (t > 1.0f) t = 1.0f;",
           "/* clamp removed */;")]),
        ("curvature_factor", "curvature formula factor doubled",
         [("return 2.0f * cross / (la * lb * lc);",
           "return 4.0f * cross / (la * lb * lc);")]),
        ("length_map_scale", "arc-length map shading scale changed (image hash)",
         [("unsigned char v = (unsigned char)(as * 36.0f);",
           "unsigned char v = (unsigned char)(as * 24.0f);")]),
        ("csv_sample_stride", "CSV emits half the constant-speed samples",
         [("for (int k = 0; k <= 14; k++) {", "for (int k = 0; k <= 14; k += 2) {")]),
    ],
    "curve": [
        ("hermite_h00_broken", "Hermite h00 basis polynomial corrupted",
         [("float h00 = 2 * u3 - 3 * u2 + 1;", "float h00 = 2 * u3 - 3 * u2;")]),
        ("clamp_disabled_left", "left clamp returns interior value",
         [("if (t <= c->keys[0].t) return c->keys[0].v;",
           "if (0) return c->keys[0].v;")]),
        ("deriv_step_wrong", "finite-difference derivative step 100x too big",
         [("float h = 1e-3f;", "float h = 0.1f;")]),
        ("order_check_inverted", "unsorted key insertion accepted",
         [("if (c->count > 0 && t <= c->keys[c->count - 1].t) return -1;",
           "if (c->count > 0 && t < c->keys[c->count - 1].t) return -1;")]),
        ("raster_gap", "curve raster skips pixels (image-only, hash catch)",
         [("for (int k = 0; k < 32; k++) {", "for (int k = 0; k < 32; k += 2) {")]),
        ("csv_range_shrunk", "CSV sampling loop covers less of the clamp region",
         [("for (int k = -4; k <= 28; k++) {", "for (int k = 0; k <= 28; k++) {")]),
    ],
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def run(cmd, cwd, timeout=120):
    proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                          timeout=timeout)
    return proc.returncode, proc.stdout + proc.stderr


def parse_checks(output):
    names = re.findall(r"^check (\w+)=(\d)$", output, re.M)
    passed = sum(1 for _, v in names if v == "1")
    return passed, len(names), dict(names)


def artifact_hashes(art_dir):
    out = {}
    if not os.path.isdir(art_dir):
        return out
    for name in sorted(os.listdir(art_dir)):
        p = os.path.join(art_dir, name)
        if os.path.isfile(p):
            out[name] = sha256_file(p)
    return out


def apply_defects(src_path, replacements):
    text = open(src_path).read()
    for old, new in replacements:
        if old not in text:
            return False
        text = text.replace(old, new, 1)
    open(src_path, "w").write(text)
    return True


def module_bench(module, work_root, baseline_only=False):
    src_rel = os.path.join(module, "src", module + ".c")
    work = os.path.join(work_root, module)
    shutil.copytree(os.path.join(ROOT, module), os.path.join(work, module))
    src = os.path.join(work, module, "src", module + ".c")
    art_dir = os.path.join(work, "artifacts")
    env = dict(os.environ, ARTIFACT_DIR=art_dir)

    def make_target(target):
        proc = subprocess.run(["make", target], cwd=os.path.join(work, module),
                              capture_output=True, text=True, env=env, timeout=120)
        return proc.returncode, proc.stdout + proc.stderr

    # clean baseline build + run
    make_target("clean")
    rc, out = make_target("test")
    base_passed, base_total, base_detail = parse_checks(out)
    base_hash = artifact_hashes(art_dir)
    baseline = {
        "module": module, "rc": rc, "passed": base_passed, "total": base_total,
        "hashes": base_hash, "detail": base_detail, "artifacts": sorted(base_hash),
    }
    if baseline_only:
        return baseline, []

    results = []
    for did, desc, repls in DEFECTS[module]:
        # fresh copy per defect so replacements always apply to pristine source
        dwork = os.path.join(work_root, module + "_" + did)
        shutil.copytree(os.path.join(ROOT, module), os.path.join(dwork, module))
        dsrc = os.path.join(dwork, module, "src", module + ".c")
        dart = os.path.join(dwork, "artifacts")
        denv = dict(os.environ, ARTIFACT_DIR=dart)
        applied = apply_defects(dsrc, repls)
        if not applied:
            results.append({"id": did, "desc": desc, "applied": False,
                            "caught": None, "how": "REPLACEMENT_FAILED"})
            continue

        def dmake(target):
            proc = subprocess.run(["make", target], cwd=os.path.join(dwork, module),
                                  capture_output=True, text=True, env=denv, timeout=120)
            return proc.returncode, proc.stdout + proc.stderr

        dmake("clean")
        brc, bout = dmake("test")
        dpassed, dtotal, ddetail = parse_checks(bout)
        dhash = artifact_hashes(dart)

        if brc != 0 and dtotal == 0:
            caught, how = True, "build_failure"
        elif brc != 0:
            caught, how = True, "check_failure(%d/%d)" % (dpassed, dtotal)
        elif dhash != base_hash:
            drifted = sorted(set(base_hash) | set(dhash)
                             if base_hash.keys() ^ dhash.keys()
                             else [k for k in base_hash if base_hash[k] != dhash.get(k)])
            caught, how = True, "artifact_hash_drift(" + ",".join(drifted) + ")"
        elif set(dhash) != set(base_hash):
            caught, how = True, "artifact_missing_or_extra"
        else:
            caught, how = False, "NOT_CAUGHT"
        results.append({"id": did, "desc": desc, "applied": True,
                        "caught": caught, "how": how,
                        "passed": dpassed, "total": dtotal})
    return baseline, results


def main():
    work_root = tempfile.mkdtemp(prefix="selftest_")
    print("work root:", work_root)
    print("=" * 78)
    baselines = {}
    all_results = {}
    for m in MODULES:
        baseline, results = module_bench(m, work_root)
        baselines[m] = baseline
        all_results[m] = results

    # ---- evidence block 1: pass counts ----
    print("\n[1] BASELINE PASS COUNTS (healthy source, isolated builds)")
    print("-" * 78)
    grand_passed = grand_total = 0
    for m in MODULES:
        b = baselines[m]
        grand_passed += b["passed"]; grand_total += b["total"]
        flag = "OK " if b["rc"] == 0 and b["passed"] == b["total"] else "BAD"
        print("%-11s %s  checks %2d/%2d  artifacts=%d"
              % (m, flag, b["passed"], b["total"], len(b["hashes"])))
    print("-" * 78)
    print("TOTAL checks passed: %d/%d" % (grand_passed, grand_total))

    # ---- evidence block 2: artifact hashes ----
    print("\n[2] VISIBLE ARTIFACT HASHES (sha256, re-opened from disk)")
    print("-" * 78)
    for m in MODULES:
        print(m + ":")
        for name, h in baselines[m]["hashes"].items():
            print("  %-24s %s" % (name, h))

    # ---- evidence block 3: capture rate table ----
    print("\n[3] DEFECT INJECTION CAPTURE RATE (%d defects/module)"
          % len(DEFECTS[MODULES[0]]))
    print("-" * 78)
    total_def = total_caught = 0
    weak = []
    for m in MODULES:
        res = all_results[m]
        applied = [r for r in res if r["applied"]]
        caught = [r for r in applied if r["caught"]]
        rate = len(caught) / len(applied) if applied else 0.0
        total_def += len(applied); total_caught += len(caught)
        tag = "PASS" if rate >= CAPTURE_FLOOR else "BELOW_80%"
        print("%-11s %d/%d caught = %5.1f%%  %s"
              % (m, len(caught), len(applied), 100.0 * rate, tag))
        for r in applied:
            mark = "+" if r["caught"] else "MISS"
            print("    [%s] %-20s %s" % (mark, r["id"], r["how"]))
        if rate < CAPTURE_FLOOR:
            weak.append((m, [r["id"] for r in applied if not r["caught"]]))
    overall = total_caught / total_def if total_def else 0.0
    print("-" * 78)
    print("OVERALL capture: %d/%d = %.1f%% (floor %.0f%%)"
          % (total_caught, total_def, 100.0 * overall, 100.0 * CAPTURE_FLOOR))

    # ---- evidence block 4: missed defects below 80% ----
    print("\n[4] MODULES BELOW 80% CAPTURE RATE - MISSED DEFECTS")
    print("-" * 78)
    if not weak:
        print("none - every module meets the 80% floor")
    else:
        for m, ids in weak:
            print("%-11s missed: %s" % (m, ", ".join(ids)))

    # ---- evidence block 5: reproduction commands ----
    print("\n[5] REPRODUCTION COMMANDS")
    print("-" * 78)
    print("# full bench (isolated temp worktrees, no repo pollution):")
    print("python3 -I tests/selftest.py")
    print("# single module healthy self-test with visible artifacts:")
    for m in MODULES:
        print("make -C %s clean test   # artifacts -> artifacts/run_%s/" % (m, m))
    print("# reproduce one injected defect by hand, e.g. skeleton rot_sign_x:")
    print("cp -r skeleton /tmp/sk && sed -i 's/\\*ox = x \\* c - y \\* s;/*ox = x * c + y * s;/' \\\n"
          "    /tmp/sk/src/skeleton.c && make -C /tmp/sk clean test")

    shutil.rmtree(work_root, ignore_errors=True)
    ok = (grand_passed == grand_total and overall >= CAPTURE_FLOOR
          and all(b["rc"] == 0 for b in baselines.values()) and not weak)
    print("\nRESULT:", "GREEN" if ok else "RED")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
