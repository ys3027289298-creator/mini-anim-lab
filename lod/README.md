# lod

按屏幕占比驱动的骨骼 LOD。`make test` 跑自检。

## 数据与流程

- `LodSkeleton`：父子层级 + 局部变换；父节点索引必须小于子节点，`lod_build` 先跑前向运动学得到全身姿态。
- `screen[i]`：每根骨的屏幕占比（必须严格大于 0，为零或负值拒绝，返回 `LOD_ERR_SCREEN`）。
- `LodConfig.importance[i]`：可配的重要性权重；排序分数 = 屏幕占比 × 重要性，降序，索引小者兜底。
- `LodTier`：按相机距离分档，每档可配 `mergeLevels`（允许塌缩的叶端层数）与 `maxError`（末端误差预算）。
  - 档位上界必须严格递增，末档上界可用 `INFINITY`。
  - `mergeLevels <= 0` 或超过根到叶的最大叶距一律拒绝（`LOD_ERR_LEVELS`），保证根与分支根不会被并掉。
- 合并：叶距（到最近后代叶子的边数）小于 `mergeLevels` 的非根骨被并入最近的保留祖先；
  分支场景下沿 `target` 指针向上走，不会跨分支串骨。
- 蒙皮权重按 `target` 重分配到保留骨，每行输入必须和为 1（±1e-5，否则 `LOD_ERR_WEIGHT`），
  压缩后再次校验和为 1，失败返回 `LOD_ERR_SUM`。
- `lod_compare`：只统计末端（叶子骨）的位置误差，按当前距离档位累计到 `LodReport`，
  超 `maxError` 置 `overBudget`；不启用 LOD 时不产生统计。
- `LodConfig.enabled = 0`：输出与输入骨骼同构（全部保留、target 指向自己），
  姿态逐字段等于全身 FK，权重原样拷贝——即不启用 LOD 时姿态与现在相同。
