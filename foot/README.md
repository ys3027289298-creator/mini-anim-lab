# foot

下坡脚步接触求解。`make test` 跑自检。

每帧管线（`foot_solve`）：

1. 由踝点、脚掌世界角 `foot_ang` 和脚掌长 `foot_len` 求期望落点（脚掌末端）。
2. 用分段线性高度场 `Ground` 查该 x 处地面高度，把落点 y 钉到地面；x 与朝向保持不变。
   查询越界或非有限数 → `FC_GROUND_ILLEGAL`。
3. 按同一朝向反解踝点。
4. 两骨 IK（与 `ik/` 相同的解析解与极向量约定）求膝、小腿；骨长非法或够不到 →
   `FC_OUT_OF_REACH`。
5. 小腿扭转 `twist = wrap_pi(foot_ang - shin_ang)`，即让脚掌保持 `foot_ang`
   所需的踝部局部转角，规范到 `[-pi, pi]`；`shin_ang` 是膝→踝的世界角。

手工脚位比对（`foot_compare`）：修正只重写落点 y，所以每帧的修正量等于该帧
手工落点相对地面的高度误差，`max_manual_ground` 即为接触点高度误差上界；
求解后的接触点 `ground_err` 应为 0。

拒绝码：`FC_OK=0`、`FC_BAD_INPUT=1`、`FC_GROUND_ILLEGAL=2`、`FC_OUT_OF_REACH=3`。
