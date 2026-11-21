# wave

弹簧/约束链上的冲击传播。距离约束等效为刚度 k 的弹簧，长波解析波速
`c = h*sqrt(k/m)`，每帧按 CFL 条件 `c*dt_sub/h <= 1` 拆子步
（`wave_required_substeps`），velocity-Verlet 积分。

- `wave_speed` / `wave_arrival_time`：解析波速与前沿到达粒子 i 的时间
- 每个粒子记录波前首次到达的帧（`arrive_frame`）与子步插值时刻
  （`arrive_time`），并给出到达帧误差上界 `wave_frame_err_bound`
- 子步加倍（2→4→8→16）时实测波速偏差按 O(dt²) 收敛
- 刚度非正、子步数非正（含零）、质量/间距/帧时长非正均拒绝并返回错误码
- 子步不足违反 CFL 时能量发散；满足 CFL 时能量守恒、传播与帧率无关

`make test` 跑自检。
