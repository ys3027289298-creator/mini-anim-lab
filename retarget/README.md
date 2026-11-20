# retarget

`make test` 跑自检。

- **提取**：`motion_extract` 把每帧拆成根轨迹与相对根的姿态轨迹（脚相对髋），可原样重组。
- **缩放**：`retarget_run` 用 `target_leg_length / source_leg_length` 缩放根与姿态，脚不会陷入地面。
- **相位**：按左脚触地区间测出源相位与周期，把触地重采样到 `target_contact_phase`；每绕过一个周期根骨推进一个步幅向量，支撑脚不打滑。
- **落点比对**：`retarget_compare` 把重定向后的落脚点与手工重定向结果逐次比对，`landing_error` 给出最大偏差；支撑脚漂移由 `slide` 度量，并对 `slide_limit` 设上界。
- **报错**（输出不写回）：骨架比例为 0 / 非有限 → `RET_ERR_SCALE`；相位差超过半周期（`phase_tolerance`，默认 0.5）→ `RET_ERR_PHASE`；滑步超限 → `RET_ERR_SLIDE`。
- **关闭重定向**：`enabled = 0` 时逐字节原样拷贝，播放与原动画一致。
