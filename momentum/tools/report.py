#!/usr/bin/env python3
"""Build REPORT.md from the momentum CSVs and figures.

Reads runs.csv / joints.csv / contacts.csv / path.csv plus the four figures
under ARTIFACT_DIR and writes ARTIFACT_DIR/REPORT.md: the ten-trial table, the
per-joint ranking, the mass-distribution comparison and the reproduction steps.

Usage: python3 -I report.py ARTIFACT_DIR
"""
import csv
import sys
from pathlib import Path

SKEL = {"A": "A（躯干重、四肢轻）", "B": "B（四肢重、躯干轻）"}


def rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def f(r, k):
    return float(r[k])


def md_table(header, body, align=None):
    align = align or ["---"] * len(header)
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join(align) + " |"]
    out += ["| " + " | ".join(row) + " |" for row in body]
    return "\n".join(out)


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    root = Path(sys.argv[1])
    runs = rows(root / "runs.csv")
    joints = rows(root / "joints.csv")
    contacts = rows(root / "contacts.csv")

    L = []
    L.append("# 急停滑移动量补偿：十组数据与关节补偿排序")
    L.append("")
    L.append("跑步急停时骨盆还在前滑。把四个地面接触点（左右脚跟、左右脚掌）的滑动量")
    L.append("加总为 S；骨盆在急停段前移 D = 2.0 m，当 S 超过 D 的一半（阈值 1.0 m）")
    L.append("触发动量补偿，把超出量 E = S − 阈值按各关节质量的**反比例**分配回去：")
    L.append("")
    L.append("```")
    L.append("quota_j = E · (1/m_j) / Σ_k (1/m_k)      Σ_j quota_j = E")
    L.append("residual = S − Σ_j quota_j = min(S, 阈值)")
    L.append("```")
    L.append("")

    L.append("## 1. 十组数据")
    L.append("")
    L.append("两套骨架总质量都是 70 kg，层级相同，只有质量分布不同；同一组五遍入区速度")
    L.append("（5.0–6.4 m/s）各跑一次，共 2 × 5 = 10 组。")
    L.append("")
    body = []
    for r in runs:
        body.append([
            r["trial_label"], r["speed"], r["total_slide"],
            "是" if int(r["triggered"]) else "否", r["excess"],
            r["residual"], f"{f(r, 'residual_ratio') * 100:.1f}%",
            r["grounded"], r["residual_grounded"],
        ])
    L.append(md_table(
        ["组", "入区速度 m/s", "滑动总量 S m", "触发", "超出量 E m",
         "补偿后残留 m", "残留/S", "可落地量 m", "仅腿部落地时残留 m"],
        body,
        ["---", "---:", "---:", "---", "---:", "---:", "---:", "---:", "---:"]))
    L.append("")
    L.append("要点：")
    L.append("")
    L.append("- **滑动总量与质量分布无关**，两套骨架逐组相同：S 只由速度、摩擦和接触点")
    L.append("  的错开时序决定，质量分布改变的是**谁**来补偿，不是**要补多少**。")
    L.append("- 5.0 与 5.2 m/s 两组未触发（S 低于 1.0 m）；5.6 m/s 起三组全部触发。")
    L.append("- **补偿后残留恒等于阈值 1.0 m**：规则要求把 E 全部分配回去，于是触发后的")
    L.append("  残留被压平在 1.0 m，不再随速度增长。这正是这条规则的口径；最后两列的")
    L.append("  敏感性视图说明它在物理上并不完整（见第 3 节）。")
    L.append("")

    L.append("## 2. 谁承担了最多的补偿量")
    L.append("")
    L.append("按反比例分配，补偿量最大的永远是**最轻的关节**。两套骨架五遍累计：")
    L.append("")
    L.append("（两套骨架的累计总量都是 3.652 m —— Σquota = E，总量只由速度决定，")
    L.append("质量分布改变的是这 3.652 m 怎么切分，不是切多少。）")
    L.append("")
    for s in ("A", "B"):
        js = [j for j in joints if j["skeleton"] == s]
        qsum = sum(f(j, "quota_total") for j in js)
        ranked = sorted(js, key=lambda j: -f(j, "quota_total"))
        L.append(f"**骨架 {SKEL[s]}** — 五遍累计补偿 {qsum:.3f} m")
        L.append("")
        body = []
        for j in ranked[:8]:
            body.append([j["joint"], j["mass"], f"{f(j, 'quota_total'):.3f}",
                         f"{f(j, 'quota_share'):.1f}%",
                         f"{f(j, 'grounded_total'):.3f}",
                         f"{f(j, 'grounded_share'):.1f}%"])
        L.append(md_table(
            ["关节", "质量 kg", "补偿配额 m", "占补偿总量", "可落地 m", "占落地总量"],
            body, ["---", "---:", "---:", "---:", "---:", "---:"]))
        L.append("")
    L.append("结论：")
    L.append("")
    qA = {j["joint"]: f(j, "quota_total") for j in joints if j["skeleton"] == "A"}
    qB = {j["joint"]: f(j, "quota_total") for j in joints if j["skeleton"] == "B"}
    tA = sum(qA.values())
    L.append(f"- 骨架 A 的前两名是 **hand_l / hand_r**（各 {qA['hand_l']:.3f} m，")
    L.append(f"  合计 {200 * qA['hand_l'] / tA:.1f}%），因为手只有 0.45 kg；骨架 B 的手")
    L.append(f"  加重到 1.1 kg，首位变成 **hand_l / hand_r / foot_l / foot_r 四者并列**")
    L.append(f"  （各 {qB['hand_l']:.3f} m，合计 {400 * qB['hand_l'] / tA:.1f}%）。")
    L.append(f"- 两套骨架里 **foot（脚）都是腿链中承担最多的关节**：A 为")
    L.append(f"  {qA['foot_l']:.3f} m（{100 * qA['foot_l'] / tA:.1f}%），B 为 "
             f"{qB['foot_l']:.3f} m（{100 * qB['foot_l'] / tA:.1f}%）。腿链越轻，脚分到的越多。")
    L.append(f"- 质量分布的杠杆很长：手从 0.45 kg 变成 1.1 kg（×2.4），配额从")
    L.append(f"  {100 * qA['hand_l'] / tA:.1f}% 掉到 {100 * qB['hand_l'] / tA:.1f}%，"
             f"同时脚的配额从 {100 * qA['foot_l'] / tA:.1f}% 涨到 "
             f"{100 * qB['foot_l'] / tA:.1f}%。")
    L.append("")

    L.append("## 3. 规则的盲点：配额落在手臂上，但滑移只发生在地面")
    L.append("")
    L.append("按规则分配完，E 是**全部收回**的，所以规则口径的残留等于阈值。但滑动只")
    L.append("发生在四个地面接触点，修正只有沿着**接触链**（腿，加上作为两条腿共同根")
    L.append("节点的骨盆）才能作用到地面上。落在脊柱、头、手臂上的配额会改变身体姿态，")
    L.append("却碰不到滑动的地面。")
    L.append("")
    L.append("把这一层算清楚（每侧腿最多只能收回自己接触点累计的滑移量，骨盆配额按")
    L.append("左右各一半计入），就得到上表最后两列。五遍合计：")
    L.append("")
    body = []
    for s in ("A", "B"):
        mine = [r for r in runs if r["skeleton"] == s]
        S = sum(f(r, "total_slide") for r in mine)
        R = sum(f(r, "residual") for r in mine)
        G = sum(f(r, "grounded") for r in mine)
        RG = sum(f(r, "residual_grounded") for r in mine)
        share = f(mine[-1], "leg_quota_share")
        body.append([f"骨架 {s}", f"{S:.3f}", f"{R:.3f}", f"{G:.3f}", f"{RG:.3f}",
                     f"{share * 100:.1f}%"])
    L.append(md_table(
        ["骨架", "滑动总量 S", "规则口径残留", "可落地补偿", "仅腿部落地时残留",
         "E 中落在腿链的比例"],
        body, ["---", "---:", "---:", "---:", "---:", "---:"]))
    L.append("")
    legA = f([r for r in runs if r["skeleton"] == "A"][-1], "leg_quota_share")
    legB = f([r for r in runs if r["skeleton"] == "B"][-1], "leg_quota_share")
    gA = sum(f(r, "grounded") for r in runs if r["skeleton"] == "A")
    gB = sum(f(r, "grounded") for r in runs if r["skeleton"] == "B")
    rA = sum(f(r, "residual_grounded") for r in runs if r["skeleton"] == "A")
    rB = sum(f(r, "residual_grounded") for r in runs if r["skeleton"] == "B")
    L.append(f"骨架 A 只有 {legA * 100:.1f}% 的 E 落在腿链上（其余落在手臂、头、脊柱和")
    L.append(f"骨盆的另一半），五遍里真正能压住滑动的只有 {gA:.3f} m，残留 {rA:.3f} m；")
    L.append(f"骨架 B 腿更重，落在腿链上的比例升到 {legB * 100:.1f}%，可落地 {gB:.3f} m，")
    L.append(f"残留 {rB:.3f} m。")
    L.append("**质量分布越偏向四肢，脚承担得越多，能真正落地的补偿也越多。**")
    L.append("")
    L.append("这不是要改掉用户给的规则，而是它的直接后果：反比例分配天然偏向最轻的")
    L.append("末端，末端又大多不在接触链上。如果目标是压住地面滑动，配额应当只在接触")
    L.append("链内做反比例分配（脚在腿链里最轻，会立刻成为最大承担者）。")
    L.append("")

    L.append("## 4. 场景与路径")
    L.append("")
    L.append("操场上一段 7 m 跑道：5 m 固定加速段从静止加到入区速度，紧接着 2 m 急停段")
    L.append("刹停。骨盆沿这条路径运动，四个接触点按脚跟先、脚掌后、左脚先、右脚后的")
    L.append("顺序投入使用（0 / 0.03 / 0.06 / 0.09 s），因此同一时刻各点的累计滑移")
    L.append("不同。摩擦系数取 μ = 0.62，刹车需求减速度 a_req = v²/(2D)。")
    L.append("")
    L.append("![跑道与接触点滑移](track.png)")
    L.append("")
    L.append("![滑动总量与阈值](slide_vs_speed.png)")
    L.append("")
    L.append("![补偿后残留](residual.png)")
    L.append("")
    L.append("![各关节补偿量](joint_quota.png)")
    L.append("")

    L.append("## 5. 复现")
    L.append("")
    L.append("```sh")
    L.append("cd momentum")
    L.append("make test      # C 模块自检")
    L.append("make accept    # 重写 ../artifacts/momentum 下的四个 CSV")
    L.append("make verify    # 独立 Python 复算 + 重画四张图 + 重建本报告")
    L.append("```")
    L.append("")
    L.append("`tools/verify_model.py` 不共享 C 的任何代码，按规格重算每个数字并核对")
    L.append("不变量（配额总和 = E、残留 = min(S, 阈值)、只有腿链关节可落地、")
    L.append("逐接触点落地量之和等于该组落地总量、路径终点静止在停止线上）。")
    L.append("")
    L.append("数据文件：`runs.csv`（十组）、`joints.csv`（逐关节）、")
    L.append("`contacts.csv`（逐接触点）、`path.csv`（骨盆路径采样）。")
    L.append("")
    (root / "REPORT.md").write_text("\n".join(L), encoding="utf-8")
    print(f"wrote {root / 'REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
