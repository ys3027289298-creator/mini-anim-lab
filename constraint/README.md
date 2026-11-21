# constraint

`make test` 跑自检。

## 接口

- `chain_step(c, dt, gx, gy)`：单步 Verlet + 距离约束求解（固定步长用）。
- `chain_step_frame(c, frame_dt, gx, gy, &policy)`：渲染帧驱动，按
  `StepPolicy` 自适应子步，低帧率防炸：
  - `target_dt` 目标子步长，子步数 `ceil(frame_dt/target_dt)`，受
    `max_substeps` 上限保护；`frame_dt` 或 `target_dt` 非正时子步数为零，
    函数打印错误并返回 -1；
  - `max_disp` 限制每个子步单粒子最大位移；
  - `min_iterations` 求解迭代下限（`c->iterations` 低于它时按它算）。
- 阻尼按 1/60 s 定义、随 dt 缩放（`powf(damping, dt*60)`），不同帧率/
  子步数下行为一致。

自检里以 240 Hz 为参考、20 Hz 为低帧率：前 0.5 s 末端位置误差上界
0.02（实测约 0.010，naive 单步约 0.114），全程末端二阶差分抖动上界
0.05（实测约 0.030，naive 约 0.072）。摆动阶段链条是混沌的，超过
0.5 s 的逐帧姿态比对没有意义，故只比前 0.5 s。
