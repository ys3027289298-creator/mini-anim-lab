# culling

层级 AABB 视锥剔除，`make test` 跑自检。

## 流程

1. `hc_validate` — 父索引/自环/成环检测（DFS 三色标记），成环直接报错拒绝剔除
2. `hc_rebuild` — 按 DFS 后序自底向上合并重建世界 AABB
3. `hc_cull` — 根向下逐层视锥裁剪：OUTSIDE 剪整棵子树，INSIDE 子树直接通过，INTERSECT 继续下探
4. `hc_reference` — 逐骨骼（各节点自身 AABB）粗测，作为真值可见集
5. `hc_compare` — 比对漏剔/误剔，并给出结构性上界
6. `HcCache` — 缓存上一帧包围盒与结论；盒未变而结论不同即检出 `HC_ERR_CACHE`

## 上界

- **误剔（wrong cull）上界 = 0**：父盒包含子盒时，自身盒与视锥相交的节点不可能被层级裁剪掉
- **漏剔（missed cull）上界 = 所有 INTERSECT 节点的子树规模之和**：保守过度覆盖只可能发生在这些区域内

## 视锥

2D 视锥用半平面集合表示：`frustum_rect` 为矩形可视区，`frustum_camera` 为带
near/far 和 FOV 侧棱的梯形相机视锥。
