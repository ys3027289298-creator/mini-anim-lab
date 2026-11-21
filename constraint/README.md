# constraint

Verlet 粒子链：距离约束 + 迭代收敛，并新增关节锥形限位。

`make test` 跑自检。

## 锥形限位

2D 链中，内点 `i`（`1..n-2`）的锥轴取父段方向 `p[i]-p[i-1]`，子段
`p[i+1]-p[i]` 被限制在与该轴的夹角 `|θ| ≤ max_angle` 内。

| 接口 | 作用 |
|---|---|
| `chain_set_cone(c, joint, max_angle)` | 配置锥半角；`max_angle < 0` 或关节越界返回错误 |
| `chain_cone_limit(c, joint, &removed, &drift)` | 单关节投影 |
| `chain_project_cone_limits(c, &removed, &drift)` | 全部已配置关节投影 |
| `chain_apply_cone_limits(c, &removed, &drift)` | **先投影、再解链约束**（推荐入口） |
| `chain_cone_excess(c, joint)` | 距锥面的剩余越界角度，锥内为 0 |
| `chain_cone_limit_manual(...)` | 手工角度截断参考实现，供比对 |

投影方式：保持 `|dir|` 不变，仅把子段方向转回锥面（取最近边界），再把
关节之后的整条子链绕关节**刚性旋转**同一角度；`px/py` 同步旋转，Verlet
速度保持一致。由于旋转保长，之后的距离约束迭代不会再引入越界。

错误码：`CONE_ERR_NEG_ANGLE`（负锥角）、`CONE_ERR_ZERO_AXIS`（轴退化）、
`CONE_ERR_ZERO_DIR`（方向退化）、`CONE_ERR_BAD_JOINT`（关节下标越界）。

漂移上界：投影使 `|proj - old| ≤ removed * r_max`，其中 `r_max` 是子链
距关节的最大半径（弦长不超过弧长）。测试同时校验实现报告的 `drift` 不超
过该解析上界，且与手工截断结果逐点一致（< 1e-5）。
