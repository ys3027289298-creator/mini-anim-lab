# timewarp

按速度曲线做时间重映射：源时间均匀重采样 → 雅可比修正墙钟间隔 → 单调三次
（PCHIP）重拟合，消除慢放时的"卡带"感。`make test` 跑自检，输出可直接当报告看。

## 流水线

1. **速度曲线** `v(s)`：定义在源时间上，`v = 1` 为原速，`0.25` 为 1/4 速。
   关键帧之间用单调三次插值（Fritsch–Carlson），采样处处为正则插值处处为正，
   墙钟映射 `T(s) = ∫ ds/v(s)` 因而严格单调、可逆。
2. **重采样**：在 **源时间** 上等距取 `ds`（每个样本覆盖相同的动画行程）。
3. **雅可比修正**：相邻样本的墙钟间隔取 `dt = ds / v(s)`，而不是 `ds`。
   慢段落墙钟步长变大、帧数变多，快段落变小、帧数变少——这正是"慢放不卡带"
   与"按墙钟等距采样会丢帧"的区别。测试里 `wall_step_slow > wall_step_fast`。
4. **重拟合**：对不规则采样 `(t_i, s_i)` 再跑一次单调三次，得到可任意墙钟时刻
   求值的 `s(t)`，C1 连续、不过冲、不回退。

积分 `∫ds/v` 用 8 点 Gauss–Legendre 逐段计算并做前缀和，是解析参考；
`tw_source_time()` 用二分求逆，不依赖重采样网格。

## 自检覆盖

| 项 | 输出/断言 |
|---|---|
| 速度含负值 | `reject_negative_speed`，返回 `TW_ERR_SPEED` |
| 速度含零 | `reject_zero_speed` |
| 速度非有限 | `reject_nan_speed` / `reject_inf_speed` |
| 时间不可逆（相等或倒退） | `reject_nonmonotone_time` / `reject_time_reversal`，`TW_ERR_TIME` |
| 查询超出范围 | `inverse_rejects_out_of_range` / `..._past_end`，`TW_ERR_RANGE` |
| 恒速解析解 | `const_speed_wall_time_exact`、`const_speed_inverse_exact`（误差 < 1e-12） |
| 积分即 1/v | `integral_matches_inverse_speed` |
| 雅可比间隔恒定 | `jacobian_spacing_ds`（`dt·v == ds`，误差 < 1e-6） |
| **跳变点为零** | `jump_points_grid` 与 `jump_points_fit` 均为 0 |
| 拟合前后位置偏差上界 | `fit_drift <= fit_drift_bound = (8/9)·h_max·v_max` |
| 对参考积分的位置误差上界 | `position_err_cubic <= position_err_bound = h_max²/8 · v_max · max|v'|` |
| 三次不劣于线性 | `cubic_beats_linear` |
| 单调不回退 | `fit_is_monotone` |
| 一千关键帧重采样耗时 | `keycount,resample_us` 表（250/500/1000/2000/4000） |

`jump_points` 的定义：相邻两段割线斜率相对差超过 `0.05` 的内部样本数。
正确做雅可比修正并单调拟合后，割线斜率就是 `v(s)`（差 O(ds)），所以严格为 0。

## 复杂度与耗时

重采样走一遍关键帧 + 一遍样本，`O(n_keys + n_samples)`；样本数
`n_samples ≈ (s_end - s_start)/ds` 与关键帧数无关，所以耗时随关键帧数线性增长，
斜率很小。实测 `resample_us` 随 `keycount` 从 250→4000 近似线性，
一千关键帧在毫秒以下量级（具体数值见运行输出）。

## 接口

```c
TwStatus tw_speed_init(TwSpeed*, const TwKey* keys, int n);  /* 校验并预积分 */
double   tw_speed_at(const TwSpeed*, double s);              /* v(s) */
double   tw_wall_time(const TwSpeed*, double s);             /* T(s) 解析积分 */
TwStatus tw_source_time(const TwSpeed*, double t, double*); /* T^{-1}(t)，二分 */
TwStatus tw_resample(const TwSpeed*, double ds, TwRetime*); /* 重采样 + 雅可比修正 */
double   tw_retime_eval(const TwRetime*, double t);         /* 播放器查询 s(t) */
```

`const char* tw_strerror(TwStatus)` 给出负速度/时间倒退的中文可读原因（英文，
便于日志）。
