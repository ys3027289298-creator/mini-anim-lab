# ik

`make test` 跑自检。

- `solve_two_bone`：无状态解析解，极向量定弯曲侧。
- `solve_two_bone_anim`：时序求解，跟踪上一帧弯曲侧，极点附近按滞回（`pole_eps` 死区）切换，切换时 smoothstep 插值过渡；`pole_eps <= 0` 返回 `err`。

