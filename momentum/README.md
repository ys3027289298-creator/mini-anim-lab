# momentum

急停时骨盆还在前滑：把地面接触点的滑动量加总，超过骨盆位移一半就触发动量
补偿，把超出部分按各关节质量的**反比例**分配回去。`make test` 跑自检（44 项
断言：质量表、滑动模型、触发阈值、反比例配额、接触链可达性、路径与确定性）。

## 规则

四个接触点（左右脚跟、左右脚掌）在急停段各需提供 `a_req = v²/(2D)` 的减速度，
摩擦最多给 `μ·g`，差额即为该接触点执不住、会向前滑的减速度：

```
a_slip  = max(0, a_req − μ·g)
T       = v / a_req
slip_i  = ½ · a_slip · max(0, T − delay_i)²      四个点错开投入使用
S       = Σ slip_i
阈值     = D / 2 = 1.0 m                          D = 2.0 m 为急停段骨盆位移
E       = max(0, S − 阈值)                        触发后才非零
quota_j = E · (1/m_j) / Σ_k (1/m_k)               Σ_j quota_j = E
残留     = S − Σ_j quota_j = min(S, 阈值)
```

## 两套骨架

总质量都是 70 kg、层级相同，只改质量分布：

- **A**（躯干重、四肢轻）：手 0.45 kg，脚 1.75 kg
- **B**（四肢重、躯干轻）：手 1.1 kg，脚 1.1 kg

同五档入区速度各跑一遍，2 × 5 = 10 组。滑动总量与质量分布无关（两套骨架逐组
相同）；质量分布只改变**谁**来补偿。

## 谁承担了最多

反比例分配永远偏向最轻的关节：

| 骨架 | 补偿最多 | 腿链中最多 |
|---|---|---|
| A | hand_l / hand_r 各 23.4% | foot_l / foot_r 各 6.0% |
| B | hand_l / hand_r / foot_l / foot_r 各 13.8% | foot_l / foot_r 各 13.8% |

但滑移只发生在地面，只有接触链（腿 + 作为共同根节点的骨盆）上的修正能作用到
接触点。模块额外报告这一层（`grounded_*` 列，属敏感性分析，不是规则口径的
残留）。细节与数字见 `../artifacts/momentum/REPORT.md`。

## 命令

```sh
make test     # 自检
make accept   # 重写 ../artifacts/momentum 下 runs/joints/contacts/path.csv
make verify   # 独立 Python 复算 + 重画四张图 + 重建 REPORT.md
make clean
```

`tools/verify_model.py` 不共享 C 的任何代码；`make verify` 依赖 matplotlib 与
Noto Sans CJK 字体，`make test` 只依赖 C 工具链。
