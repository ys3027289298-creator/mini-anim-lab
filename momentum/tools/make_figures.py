#!/usr/bin/env python3
"""Render the momentum module's figures from its CSV artefacts.

  slide_vs_speed.png - S per trial. Slide is mass-independent, so the two
                       skeletons coincide; drawn once with that stated.
  residual.png       - residual after compensation: the rule view (flat at the
                       threshold once triggered) and the grounded sensitivity
                       view, where the two skeletons separate.
  joint_quota.png    - per joint: inverse-mass quota vs correction that can act
                       on the ground, for both skeletons.
  track.png          - the playground lane, the four contact footprints and
                       their slip, plus the pelvis speed profile.

Palette: reference categorical slots 1-3 (blue, orange, aqua), validated for
both modes; a sequential light step for the de-emphasised quota. Text uses ink
tokens only; series colour appears on marks and line keys, never on labels.

Usage: python3 -I make_figures.py ARTIFACT_DIR
"""
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.patches import FancyArrowPatch, Rectangle

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
QUOTA_FILL = "#d8d7d0"      # de-emphasis gray for the allocated quota
SKEL_COLOR = {"A": BLUE, "B": ORANGE}
FONT = "Noto Sans CJK SC"


def setup_style():
    try:
        fm.fontManager.addfont("/home/ys/.fonts/NotoSansCJKsc-Regular.otf")
        plt.rcParams["font.family"] = FONT
    except Exception:
        pass
    plt.rcParams.update({
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "axes.edgecolor": BASELINE,
        "axes.labelcolor": INK_2,
        "text.color": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 1.0,
        "grid.linestyle": "-",
        "axes.axisbelow": True,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.size": 10,
        "axes.titlesize": 12,
        "figure.dpi": 160,
    })


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def f(r, k):
    return float(r[k])


def save(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {path}")


def fig_slide(root, runs):
    fig, ax = plt.subplots(figsize=(7.0, 3.9))
    threshold = f(runs[0], "threshold")
    mine = sorted((r for r in runs if r["skeleton"] == "A"), key=lambda r: f(r, "speed"))
    xs = [f(r, "speed") for r in mine]
    ys = [f(r, "total_slide") for r in mine]
    ax.plot(xs, ys, color=BLUE, linewidth=2, marker="o", markersize=7,
            markerfacecolor=BLUE, markeredgecolor=SURFACE, markeredgewidth=2)
    ax.axhline(threshold, color=INK_2, linewidth=1.2)
    ax.annotate(f"触发阈值 D/2 = {threshold:.2f} m",
                xy=(xs[0], threshold), xytext=(xs[0] + 0.05, threshold - 0.42),
                color=INK_2, fontsize=9,
                arrowprops=dict(arrowstyle="-", color=BASELINE, linewidth=1))
    ax.annotate(f"两套骨架逐点重合\n（滑动与质量分布无关）",
                xy=(xs[3], ys[3]), xytext=(xs[1] + 0.02, ys[4] + 0.18),
                color=INK_2, fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=BASELINE, linewidth=1))
    for i, (x, y) in enumerate(zip(xs, ys)):
        if i in (0, len(xs) - 1):
            ax.annotate(f"{y:.2f}", (x, y), textcoords="offset points",
                        xytext=(0, 10), ha="center", color=INK_2, fontsize=8.5)
    ax.set_xlabel("入区速度 (m/s)")
    ax.set_ylabel("四个接触点滑动总量 S (m)")
    ax.set_title("急停滑动总量随入区速度增长")
    ax.set_ylim(0, max(ys) * 1.2)
    ax.set_xlim(xs[0] - 0.15, xs[-1] + 0.15)
    ax.grid(axis="x", visible=False)
    save(fig, root / "slide_vs_speed.png")


def fig_residual(root, runs):
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    threshold = f(runs[0], "threshold")
    rule = sorted((r for r in runs if r["skeleton"] == "A"), key=lambda r: f(r, "speed"))
    xs = [f(r, "speed") for r in rule]
    ys = [f(r, "residual") for r in rule]
    ax.plot(xs, ys, color=INK, linewidth=2, marker="o", markersize=7,
            markerfacecolor=INK, markeredgecolor=SURFACE, markeredgewidth=2,
            label="规则口径残留 = min(S, 阈值)")
    for s in ("A", "B"):
        mine = sorted((r for r in runs if r["skeleton"] == s), key=lambda r: f(r, "speed"))
        ax.plot([f(r, "speed") for r in mine], [f(r, "residual_grounded") for r in mine],
                color=SKEL_COLOR[s], linewidth=1.6, linestyle=(0, (1, 1.6)),
                marker="s", markersize=6, markerfacecolor=SURFACE,
                markeredgecolor=SKEL_COLOR[s], markeredgewidth=1.6,
                label=f"骨架 {s} 仅腿部可落地（敏感性）")
        last = mine[-1]
        ax.annotate(f"{f(last,'residual_grounded'):.2f}",
                    (f(last, "speed"), f(last, "residual_grounded")),
                    textcoords="offset points", xytext=(8, -2), color=INK_2, fontsize=8.5)
    ax.axhline(threshold, color=INK_2, linewidth=1.2)
    ax.annotate("阈值 1.00 m：触发后规则口径残留被压平",
                xy=(xs[1] + 0.02, ys[1]), xytext=(xs[1] + 0.10, 0.42),
                color=INK_2, fontsize=8.5,
                arrowprops=dict(arrowstyle="-", color=BASELINE, linewidth=1,
                                connectionstyle="arc3,rad=-0.2"))
    ax.set_xlabel("入区速度 (m/s)")
    ax.set_ylabel("补偿后残留滑动 (m)")
    ax.set_title("补偿后残留：超出阈值的部分按规则全部收回")
    ax.set_ylim(0, max(f(r, "residual_grounded") for r in runs) * 1.28)
    ax.set_xlim(xs[0] - 0.15, xs[-1] + 0.3)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")
    ax.grid(axis="x", visible=False)
    save(fig, root / "residual.png")


def fig_joint_quota(root, joints):
    names = [r["joint"] for r in joints if r["skeleton"] == "A"]
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.7), sharex=True)
    for ax, s in zip(axes, ("A", "B")):
        sel = {r["joint"]: r for r in joints if r["skeleton"] == s}
        quota = [f(sel[n], "quota_total") for n in names]
        grounded = [f(sel[n], "grounded_total") for n in names]
        y = list(range(len(names)))
        ax.barh([i + 0.21 for i in y], quota, height=0.38, color=QUOTA_FILL,
                edgecolor="none", label="逆质量配额（E 的全部分配）")
        ax.barh([i - 0.21 for i in y], grounded, height=0.38, color=SKEL_COLOR[s],
                edgecolor="none", label="能落在接触链上的量（敏感性）")
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8.5)
        ax.invert_yaxis()
        ax.set_title(f"骨架 {s}", fontsize=11)
        ax.set_xlabel("五遍累计补偿量 (m)")
        ax.tick_params(labelleft=s == "A")
        ax.grid(axis="y", visible=False)
        pass
    axes[0].legend(frameon=False, fontsize=8.5, loc="lower right")
    fig.suptitle("补偿落在最轻的关节，但只有腿部能把修正送到地面",
                 fontsize=12, y=1.0)
    fig.tight_layout()
    save(fig, root / "joint_quota.png")


def fig_track(root, runs, contacts, path_rows):
    fig, (ax, axv) = plt.subplots(2, 1, figsize=(9.8, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [1.25, 1.0]})
    runup, stop = 5.0, 2.0

    ax.add_patch(Rectangle((0, -1.0), runup + stop, 2.0, facecolor="#f4f4f0",
                           edgecolor="none", zorder=0))
    ax.axvline(runup, color=BASELINE, linewidth=1.0, zorder=1)
    ax.text(runup / 2, 1.18, "加速段（固定 5 m）", ha="center", color=INK_2, fontsize=9)
    ax.text(runup + stop / 2, 1.18, "急停段（2 m）", ha="center", color=INK_2, fontsize=9)
    ax.text(runup, -1.42, "制动线", ha="center", color=MUTED, fontsize=8)

    fastest = max(runs, key=lambda r: f(r, "speed"))
    sel = [r for r in contacts if r["skeleton"] == "A"
           and int(r["run"]) == int(fastest["run"])]
    lane_y = {"left": 0.34, "right": -0.34}
    # label rows so the two contacts of a foot don't collide
    label_y = {"heel_l": 0.62, "ball_l": 0.34, "heel_r": -0.34, "ball_r": -0.62}
    for r in sel:
        x = f(r, "x_plant")
        y = lane_y[r["leg"]]
        slip = f(r, "slip")
        ax.plot([x], [y], marker="o", markersize=7, color=ORANGE,
                markeredgecolor=SURFACE, markeredgewidth=2, zorder=4)
        ax.add_patch(FancyArrowPatch((x, y), (x + slip, y),
                                     arrowstyle="-|>", mutation_scale=9,
                                     color=ORANGE, linewidth=1.4, zorder=4))
        ly = label_y[r["contact"]]
        ax.annotate(f"{r['contact']}  {slip:.2f} m", xy=(x + slip, y),
                    xytext=(x + slip + 0.16, ly), fontsize=8, color=INK_2,
                    va="center", arrowprops=dict(arrowstyle="-", color=BASELINE,
                                                 linewidth=0.8))
    ax.set_ylim(-1.55, 1.4)
    ax.set_yticks([])
    ax.grid(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_visible(spine == "bottom")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.set_title("操场跑道 7 m：5 m 加速接 2 m 急停，四个接触点同向滑移")

    # pelvis speed profile: all trials muted, fastest emphasised
    for r in runs:
        if r["skeleton"] != "A":
            continue
        pr = [p for p in path_rows if int(p["run"]) == int(r["run"])]
        xs = [f(p, "x") for p in pr]
        vs = [f(p, "v") for p in pr]
        fast = float(r["run"]) == float(fastest["run"])
        axv.plot(xs, vs, color=BLUE if fast else MUTED, linewidth=2 if fast else 1.2,
                 alpha=1.0 if fast else 0.55, zorder=3 if fast else 2)
        axv.annotate(f"{f(r,'speed'):.1f}", (xs[60], vs[60]),
                     textcoords="offset points", xytext=(6, 4 if fast else -10),
                     color=INK_2, fontsize=8)
    axv.axvline(runup, color=BASELINE, linewidth=1.0)
    axv.set_ylabel("骨盆速度 (m/s)")
    axv.set_xlabel("沿跑道位置 (m)")
    axv.set_xlim(0, runup + stop)
    axv.set_ylim(0, 7.0)
    axv.grid(axis="x", visible=False)
    axv.set_title("骨盆速度剖面（五遍试验）", fontsize=10, loc="left")
    save(fig, root / "track.png")


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1])
    setup_style()
    runs = rows(root / "runs.csv")
    joints = rows(root / "joints.csv")
    contacts = rows(root / "contacts.csv")
    path_rows = rows(root / "path.csv")
    fig_slide(root, runs)
    fig_residual(root, runs)
    fig_joint_quota(root, joints)
    fig_track(root, runs, contacts, path_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
