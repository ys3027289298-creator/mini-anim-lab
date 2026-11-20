# skeleton

`make test` 跑自检。

## 统一旋转表示（unified）

默认走原来的角度累加路径（`unified == 0`），姿势与旧版逐位一致。
`skel_enable_unified(&s, 1)` 后切换为四元数表示：

- 世界旋转按四元数沿父子链复合（`world_q = world_q[parent] * local_q`），每帧归一化；
- 增量旋转用 `skel_rotate_local`（累积四元数，状态有界），绝对角度用 `skel_set_local_angle`；
- 四元数模长为零时 `skel_update` 返回 -1 并置 `quat_error`；
- 自检含 10 分钟（60fps × 36000 帧）漂移测试：与每帧用 double 重建矩阵的参考姿势比对，
  末端位置误差上界由 `unified_tip_bound` 给出，阈值 1e-3。
