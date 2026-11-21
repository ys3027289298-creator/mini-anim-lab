# constraint

Verlet 距离约束链条 + 逐帧守恒校验。

`make test` 跑自检。

## 守恒校验

链条"自己飘起来"时长度往往还正常，所以校验同时看总量与逐约束：

- **逐帧统计**：每帧记录总长度 `L`、动能 `K`、势能 `U`、总能量 `E = K + U`，
  以及每条约束的残差 `|len_i - rest_i|`。
- **与理论值比对**：
  - 长度理论值 `L_ref = sum(rest_i)`，恒定；`|L - L_ref| > len_tol` 记长度违反。
  - 能量理论基线取首个观测帧的 `E`。默认策略 `CONS_ENERGY_NONINCREASE`
    只报能量**增长**超过 `en_tol`（阻尼 ≤ 1 时系统不应增能，这正是"飘起来"的特征）；
    `CONS_ENERGY_ABSOLUTE` 改为报任意方向的偏差。
- **定位到具体约束**：约束 `i` 在 `|len_i - rest_i| > len_tol` 时被点名，
  判据与手工逐约束核对完全一致。能量是全局量：能量违反但没有任何单条约束超限时，
  记为全局违反（`cons_worst_constraint` 返回 -1），不会硬栽给某一段。
- **漏报为零**：自测里把校验器的点名集合与独立的手工逐约束集合逐帧对比，
  含 4800 帧随机差分（含随机脉冲），漏报 0、误报 0。
- **采样不足**：观测帧数 `< min_frames` 且尚未看到违反时判 `undetermined`
  （"无法判定"），绝不返回 `ok`；已观测到的违反仍如实报出，否则会相对手工核对漏报。
- **非法阈值拒绝**：`len_tol`/`en_tol` 为 0、负数、NaN、Inf，`dt` 非正，
  `min_frames < 1`，`energy_policy` 越界，或观测中途改变链条拓扑，一律 `invalid`，
  且不会在后续观测中"复活"。

状态：`cons_status()` 返回 `ok` / `undetermined` / `violation` / `invalid`，
文字说明见 `cons_status_reason()`。

## 用法

```c
ConsCfg cfg = { .len_tol = 0.05f, .en_tol = 1e-2f, .min_frames = 50,
                .energy_policy = CONS_ENERGY_NONINCREASE,
                .dt = 0.016f, .gx = 0.0f, .gy = -9.8f };
ConsCheck chk;
if (!cons_init(&chk, &chain, &cfg)) { /* cfg 非法 */ }
for (each frame) { chain_step(&chain, 0.016f, 0.0f, -9.8f);
                   cons_observe(&chk, &chain); }
ConsStatus s = cons_status(&chk);
float peak; int bad = cons_worst_constraint(&chk, &peak);  /* -1 = 全局 */
```

明细缓冲上限 `MAX_FRAMES`（2048）只限制逐帧明细留存；计数、定位与最终判定不受影响，
超出时 `chk.capped = 1`。
