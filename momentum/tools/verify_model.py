#!/usr/bin/env python3
"""Independently re-derive the momentum model from its spec and check the CSVs.

The C module writes runs.csv / joints.csv / contacts.csv / path.csv. This
script carries its own copy of the constants and formulas (it shares no code
with the C) and re-computes every written number, then checks the invariants the
analysis relies on:

  * quota_j is exactly E * (1/m_j) / sum(1/m_k) and sums to E
  * residual = S - E = min(S, threshold)          (the requested residual)
  * grounded view: only leg-chain joints can deliver, grounded <= min(demand,
    available slip); per-contact grounded sums to the run's grounded total
  * the pelvis path ends at rest on the stop line

Usage: python3 -I verify_model.py ARTIFACT_DIR
Exit code 0 only if every check passes.
"""
import csv
import sys
from pathlib import Path

GRAV = 9.807
FRICTION = 0.62
STOP_DIST = 2.0
RUNUP_DIST = 5.0
CSV_EPS = 1e-6      # CSVs carry six decimals
SHARE_EPS = 1e-3    # share columns carry four decimals

JOINTS = [
    ("pelvis", "pelvis"), ("spine", "none"), ("head", "none"),
    ("upperarm_l", "none"), ("forearm_l", "none"), ("hand_l", "none"),
    ("upperarm_r", "none"), ("forearm_r", "none"), ("hand_r", "none"),
    ("thigh_l", "left"), ("shin_l", "left"), ("foot_l", "left"),
    ("thigh_r", "right"), ("shin_r", "right"), ("foot_r", "right"),
]

MASS = {
    "A": [12.0, 13.8, 5.0, 1.8, 1.1, 0.45, 1.8, 1.1, 0.45, 9.5, 5.0, 1.75, 9.5, 5.0, 1.75],
    "B": [13.0, 17.5, 6.5, 3.5, 1.8, 1.1, 3.5, 1.8, 1.1, 6.0, 3.0, 1.1, 6.0, 3.0, 1.1],
}

CONTACTS = [
    ("heel_l", 0, 0.00, RUNUP_DIST + 0.01, +0.13),
    ("ball_l", 0, 0.03, RUNUP_DIST + 0.19, +0.13),
    ("heel_r", 1, 0.06, RUNUP_DIST - 0.19, -0.13),
    ("ball_r", 1, 0.09, RUNUP_DIST - 0.01, -0.13),
]

RUNS = [5.0, 5.2, 5.6, 6.0, 6.4]

failures = []
checks = 0


def check(name, ok, detail=""):
    global checks
    checks += 1
    if not ok:
        failures.append(f"{name}: {detail}")
    print(f"{name}={1 if ok else 0}")


def near(a, b, eps=CSV_EPS):
    return abs(a - b) <= eps


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def simulate(v):
    a_req = v * v / (2.0 * STOP_DIST)
    stop_time = v / a_req
    threshold = 0.5 * STOP_DIST
    a_slip = max(0.0, a_req - FRICTION * GRAV)
    slip = [0.5 * a_slip * max(0.0, stop_time - d) ** 2 for _, _, d, _, _ in CONTACTS]
    total = sum(slip)
    triggered = total > threshold
    excess = total - threshold if triggered else 0.0
    return dict(a_req=a_req, stop_time=stop_time, threshold=threshold,
                a_slip=a_slip, slip=slip, total=total,
                triggered=triggered, excess=excess)


def weights(name):
    inv = [1.0 / m for m in MASS[name]]
    s = sum(inv)
    return [x / s for x in inv]


def compensate(name, rr):
    w = weights(name)
    quota = [rr["excess"] * x for x in w]
    avail = [0.0, 0.0]
    for i, (_, leg, _, _, _) in enumerate(CONTACTS):
        avail[leg] += rr["slip"][i]
    demand = [quota[0] / 2.0, quota[0] / 2.0]
    for j, (_, chain) in enumerate(JOINTS):
        if chain == "left":
            demand[0] += quota[j]
        elif chain == "right":
            demand[1] += quota[j]
    grounded_leg = [min(demand[l], avail[l]) for l in (0, 1)]
    ratio = [grounded_leg[l] / demand[l] if demand[l] else 0.0 for l in (0, 1)]
    grounded = [0.0] * len(JOINTS)
    for j, (_, chain) in enumerate(JOINTS):
        if chain == "left":
            grounded[j] = quota[j] * ratio[0]
        elif chain == "right":
            grounded[j] = quota[j] * ratio[1]
    grounded[0] = 0.5 * quota[0] * (ratio[0] + ratio[1])
    return dict(weights=w, quota=quota, quota_sum=sum(quota),
                available=avail, demand=demand, grounded_leg=grounded_leg,
                grounded=grounded, grounded_total=sum(grounded_leg),
                residual=rr["total"] - sum(quota),
                residual_grounded=rr["total"] - sum(grounded_leg),
                leg_share=(sum(demand) / rr["excess"]) if rr["excess"] else 0.0)


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1])
    run_rows = rows(root / "runs.csv")
    joint_rows = rows(root / "joints.csv")
    contact_rows = rows(root / "contacts.csv")
    path_rows = rows(root / "path.csv")

    check("ten_trials", len(run_rows) == 10, f"got {len(run_rows)}")
    check("per_contact_rows", len(contact_rows) == 10 * len(CONTACTS), str(len(contact_rows)))

    for r in run_rows:
        name, run = r["skeleton"], int(r["run"])
        v = RUNS[run - 1]
        rr = simulate(v)
        comp = compensate(name, rr)
        tag = f"run_{name}{run}"
        check(f"{tag}_speed", near(float(r["speed"]), v))
        check(f"{tag}_a_req", near(float(r["a_req"]), rr["a_req"]))
        check(f"{tag}_stop_time", near(float(r["stop_time"]), rr["stop_time"]))
        check(f"{tag}_threshold", near(float(r["threshold"]), rr["threshold"]))
        for i, key in enumerate(("slip_heel_l", "slip_ball_l", "slip_heel_r", "slip_ball_r")):
            check(f"{tag}_{key}", near(float(r[key]), rr["slip"][i]))
        check(f"{tag}_total_slide", near(float(r["total_slide"]), rr["total"]))
        check(f"{tag}_trigger_flag", int(r["triggered"]) == int(rr["triggered"]))
        check(f"{tag}_excess", near(float(r["excess"]), rr["excess"]))
        check(f"{tag}_quota_sum", near(float(r["quota_sum"]), rr["excess"]))
        check(f"{tag}_residual", near(float(r["residual"]), comp["residual"]))
        check(f"{tag}_residual_ratio",
              near(float(r["residual_ratio"]),
                   comp["residual"] / rr["total"] if rr["total"] else 0.0))
        check(f"{tag}_leg_quota_share", near(float(r["leg_quota_share"]), comp["leg_share"]))
        check(f"{tag}_grounded", near(float(r["grounded"]), comp["grounded_total"]))
        check(f"{tag}_residual_grounded",
              near(float(r["residual_grounded"]), comp["residual_grounded"]))
        check(f"{tag}_residual_is_min", near(comp["residual"], min(rr["total"], rr["threshold"])))
        if rr["triggered"]:
            check(f"{tag}_residual_below_total", comp["residual"] < rr["total"])
            check(f"{tag}_grounded_positive", comp["grounded_total"] > 0.0)
            check(f"{tag}_grounded_le_excess", comp["grounded_total"] <= rr["excess"] + CSV_EPS)
        else:
            check(f"{tag}_off_no_comp", comp["grounded_total"] == 0.0)

    for s in ("A", "B"):
        w = weights(s)
        quota = [0.0] * len(JOINTS)
        grounded = [0.0] * len(JOINTS)
        for run in range(1, len(RUNS) + 1):
            rr = simulate(RUNS[run - 1])
            comp = compensate(s, rr)
            quota = [a + b for a, b in zip(quota, comp["quota"])]
            grounded = [a + b for a, b in zip(grounded, comp["grounded"])]
        qsum, gsum = sum(quota), sum(grounded)
        by_joint = {r["joint"]: r for r in joint_rows if r["skeleton"] == s}
        check(f"joints_{s}_count", len(by_joint) == len(JOINTS), str(len(by_joint)))
        for j, (jn, chain) in enumerate(JOINTS):
            r = by_joint[jn]
            check(f"joint_{s}_{jn}_mass", near(float(r["mass"]), MASS[s][j]))
            check(f"joint_{s}_{jn}_weight", near(float(r["inverse_mass_weight"]), w[j]))
            check(f"joint_{s}_{jn}_chain", r["chain"] == chain)
            check(f"joint_{s}_{jn}_quota", near(float(r["quota_total"]), quota[j]))
            check(f"joint_{s}_{jn}_quota_share",
                  near(float(r["quota_share"]), 100.0 * quota[j] / qsum if qsum else 0.0,
                       SHARE_EPS))
            check(f"joint_{s}_{jn}_grounded", near(float(r["grounded_total"]), grounded[j]))
            check(f"joint_{s}_{jn}_grounded_share",
                  near(float(r["grounded_share"]), 100.0 * grounded[j] / gsum if gsum else 0.0,
                       SHARE_EPS))
            if chain == "none":
                check(f"joint_{s}_{jn}_no_grounding", grounded[j] == 0.0)

    for s in ("A", "B"):
        for run in range(1, len(RUNS) + 1):
            rr = simulate(RUNS[run - 1])
            comp = compensate(s, rr)
            sel = [r for r in contact_rows if r["skeleton"] == s and int(r["run"]) == run]
            tag = f"contacts_{s}{run}"
            check(f"{tag}_count", len(sel) == len(CONTACTS), str(len(sel)))
            check(f"{tag}_slip_sum", near(sum(float(r["slip"]) for r in sel), rr["total"]))
            check(f"{tag}_grounded_sum",
                  near(sum(float(r["grounded"]) for r in sel), comp["grounded_total"], 5 * CSV_EPS))
            check(f"{tag}_residual_sum",
                  near(sum(float(r["residual"]) for r in sel), comp["residual_grounded"], 5 * CSV_EPS))
            for i, (cname, leg, delay, x, y) in enumerate(CONTACTS):
                r = sel[i]
                ct = f"contact_{s}{run}_{cname}"
                check(f"{ct}_delay", near(float(r["delay"]), delay))
                check(f"{ct}_leg", r["leg"] == ("left" if leg == 0 else "right"))
                check(f"{ct}_plant", near(float(r["x_plant"]), x) and near(float(r["y"]), y))
                check(f"{ct}_slip", near(float(r["slip"]), rr["slip"][i]))
                check(f"{ct}_grounded_le_slip",
                      float(r["grounded"]) <= float(r["slip"]) + CSV_EPS)
                check(f"{ct}_residual",
                      near(float(r["residual"]),
                           float(r["slip"]) - float(r["grounded"]), 5 * CSV_EPS))

    check("path_rows", len(path_rows) == 5 * 121, str(len(path_rows)))
    for run in range(1, 6):
        ps = [r for r in path_rows if int(r["run"]) == run]
        tag = f"path_run{run}"
        check(f"{tag}_count", len(ps) == 121, str(len(ps)))
        check(f"{tag}_start", near(float(ps[0]["x"]), 0.0) and near(float(ps[0]["v"]), 0.0))
        check(f"{tag}_accel_end",
              near(float(ps[60]["x"]), RUNUP_DIST) and near(float(ps[60]["v"]), RUNS[run - 1]))
        check(f"{tag}_stop",
              near(float(ps[-1]["x"]), RUNUP_DIST + STOP_DIST) and near(float(ps[-1]["v"]), 0.0))
        check(f"{tag}_monotonic",
              all(float(ps[i]["x"]) >= float(ps[i - 1]["x"]) for i in range(1, len(ps))))
        check(f"{tag}_phase", int(ps[60]["phase"]) == 0 and int(ps[-1]["phase"]) == 1)

    for s in ("A", "B"):
        mine = [r for r in run_rows if r["skeleton"] == s]
        # The rule removes the whole excess, so a triggered run's residual is
        # flat at the threshold; untriggered runs keep their whole slide.
        check(f"residual_is_flat_at_threshold_{s}",
              all(near(float(r["residual"]), STOP_DIST / 2.0)
                  for r in mine if int(r["triggered"]) == 1))
        check(f"residual_never_exceeds_threshold_{s}",
              all(float(r["residual"]) <= STOP_DIST / 2.0 + CSV_EPS for r in mine))
        grounded = [float(r["grounded"]) for r in mine if int(r["triggered"]) == 1]
        check(f"grounded_grows_with_speed_{s}",
              all(b > a for a, b in zip(grounded, grounded[1:])), str(grounded))
        res_g = [float(r["residual_grounded"]) for r in mine]
        check(f"residual_grounded_grows_with_speed_{s}",
              all(b > a for a, b in zip(res_g, res_g[1:])), str(res_g))
    ra = [r for r in run_rows if r["skeleton"] == "A"]
    rb = [r for r in run_rows if r["skeleton"] == "B"]
    # B has heavier legs, so more of its correction can reach the contacts.
    check("B_grounds_more_than_A",
          all(float(b["grounded"]) >= float(a["grounded"]) for a, b in zip(ra, rb))
          and any(float(b["grounded"]) > float(a["grounded"]) for a, b in zip(ra, rb))
          and all(float(b["residual_grounded"]) <= float(a["residual_grounded"])
                  for a, b in zip(ra, rb)))
    # Same speeds, same slide; only the mass distribution differs.
    check("slide_independent_of_mass",
          all(near(float(a["total_slide"]), float(b["total_slide"])) for a, b in zip(ra, rb)))

    print(f"checks={checks}")
    if failures:
        print(f"FAILED={len(failures)}")
        for f in failures[:40]:
            print("  " + f)
        return 1
    print("all_ok=1")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
