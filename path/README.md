# path

`make test` 跑自检。

## 弧长重参数化

折线部分（`path_add` / `path_sample` / `path_closest` / `path_curvature`）保持不变，新增光滑曲线的一段：

| 函数 | 作用 |
|---|---|
| `arc_build_uniform` | 按参数等分采样，建累积弦长表 |
| `arc_build_adaptive` | 中点矢高超过 `tol` 就二分，节点密度随曲率增大 |
| `arc_lookup` | 距离 `s` → 段下标 + 段内比例（二分查找） |
| `arc_validate` | 表结构自检：长度严格递增、与点距一致 |
| `arc_position_bound` / `arc_speed_bound` | 位置误差与速度波动上界（仅用表数据） |
| `arc_length_rel_error` / `arc_max_position_error` / `arc_speed_fluctuation` | 与解析弧长比对 |
| `arc_fit_slope` | log-log 线性拟合收敛斜率 |

### 拒绝条件

- 段长为零（重复点 / 半径为零的曲线）：建表返回 `ARC_E_ZERO`，`arc_lookup` 返回 `ARC_E_ZERO`。
- 距离越界：`s < 0` 或 `s > total` 返回 `ARC_E_RANGE`；`NaN` 返回 `ARC_E_ARG`。
- 节点数超过上限：`ARC_E_CAP`。

### 自检结论

对圆（解析弧长 `r·u`）与抛物线 `y = x²/2`（解析弧长 `0.5(x√(1+x²) + asinh x)`）：
均匀采样密度 16…512 时，累积长度、位置、速度三项误差都按 **N⁻²** 下降，实测斜率偏离 −2 小于 0.005。
位置误差 ≤ 矢高 + 累积漂移上界，速度波动 ≤ `1 − cos(θ_max/2)` 上界，两条在每个采样密度下都成立。
曲率自适应在抛物线顶点附近节点密度更高（实测约 1.34 倍），其误差随节点数的斜率同样为 −2。
