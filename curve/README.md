# curve

`make test` 跑自检。

`make analysis` 跑摆臂步频实验（`src/swing.c`）：按快走/慢走/跑三种步频
各录一段摆动角度序列（CSV 写到 `../artifacts/`），对每段的角速度序列
找突变点，并对照平滑后的曲率上界。

模块功能：

- 关键帧采样、三次埃尔米特切线、夹取（原有）。
- `curve_is_extremum`：相邻高度差变号判定极值关键帧。
- `curve_seg_curvature_bound` / `curve_curvature_bound`：曲率上界。
  Hermite 段的 v'' 对参数是线性的，端点取 max 即精确上界；内部关键帧
  若 tan_in != tan_out（切线断裂）则上界为 INFINITY。
- `curve_smooth_extrema(w, kappa_max, &skipped)`：极值处切线平滑。
  窗口 +-w 内割线平均出公共切线；相邻极值落入 2w 内说明窗口相对当地
  步频过宽，跳过并计入 skipped；平滑后任一相邻段会超出 kappa_max 的
  极值把切线钳回可行区间（k 关于切线是线性的，区间精确可解），不可行
  则压平为 0——这些被单独处理的极值计入返回值。
