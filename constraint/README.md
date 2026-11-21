# constraint

布料式位置约束（PBD）绳链：Verlet 积分 + 双边距离约束，`make test` 跑自检。

## 模型

每段 `p[i]..p[i+1]` 是一条距离约束，静止长度 `rest[i]`，刚度 `stiffness[i]`
（0 关闭该段，>1 夹到 1）。求解一次约束：

```
d      = |p[i+1] - p[i]|
diff   = (d - rest) / d * stiffness
p[i]   += (p[i+1]-p[i]) * diff * (w_i   / (w_i + w_{i+1}))
p[i+1] -= (p[i+1]-p[i]) * diff * (w_{i+1} / (w_i + w_{i+1}))
```

`w = inv_mass`。修正量按质量倒数分配：钉住的点 `w = 0` 不承担任何修正，其余
粒子按惯性比例分担，因此等质量自由段的质心严格守恒。坐标用 `double`，残差可以
敛到 1e-15 量级而不会因 float 舍入卡在 1e-7。

## 收敛控制

- `substeps`：一帧切成几个子步，每个子步重新积分 + 求解。峰值残差随子步数单调下降。
- `iterations`：每个子步内最多多少次对称 Gauss–Seidel 扫描（正序 + 逆序，消除方向偏置）。
- `tolerance`：扫描过程中残差降到该值以下就提前结束（0 表示不设）。
- `chain_solve` 记录 best-so-far 并恢复它，保证「迭代次数翻倍，最大残差单调不增」。
  单纯 Gauss–Seidel 在质量/刚度差异大时最大范数并不单调，自检里有
  `raw_gs_monotone=0 guarded_monotone=1` 的对照。

每次求解/每帧都给出残差，残差定义为**所有段里最大的绝对长度误差**（长度单位）：

- `chain_residual()` — 当前构型
- `last_residual` — 上次求解最优迭代的残差；`chain_step` 结束后是末子步的残差
- `peak_residual` — 上次 `chain_step` 内各子步残差的最大值
- `last_sweeps` — 实际执行的扫描次数（tolerance 提前退出时小于 iterations）

## 悬链线与误差上界

`catenary_a / catenary_node` 给出两端等高的解析悬链线，按等弧长取样（低频节点在
最低点两侧对称）。自检把 n=21、41 的绳子放松到静止后与解析解比对，允许误差

```
bound = h^2 * kappa / 4 + 2 * |g| * dt^2
```

- 离散项：等弦长节点在解析悬链线上，弦与曲线的最大偏离 ≤ h²κ/8；这里取 1/4 是
  保守放大（实测约为 h²κ·0.03–0.18，且网格越细越小）。
- 时间项：PBD 的不动点是 Verlet 步的投影，与静态平衡差 O(dt²)。实测系数约
  6.3·dt²（首端子步），2|g|dt² 是保守上界。

实测（span=4，len=5，dt=0.01，256 次扫描）：n=21 误差 1.9e-3 vs 上界 1.8e-2；
n=41 误差 1.1e-3 vs 上界 6.1e-3。同时检查误差是二阶的：
`err < h²κ/2 + 3dt²`。

## 输入校验

- `chain_set_mass`：拒绝 0、负数、NaN、Inf；钉住时改质量只更新 `mass`，`inv_mass` 保持 0。
- `chain_set_stiffness`：拒绝负数和 NaN/Inf；大于 1 夹到 1。
- `chain_init`：拒绝 n < 2、n > MAX_P、非正/非有限段长。
- `chain_set_mass / chain_set_stiffness / chain_pin`：越界索引返回 -1。
- `chain_step`：dt ≤ 0、NaN dt 或 NaN 重力不做积分，直接返回当前残差，位置不变。

## 自检

`make test` 输出逐项 `name=1`，最后 `all_ok=1`。覆盖：质量倒数分配与质心守恒、
输入校验、子步与迭代的收敛性、best-so-far 单调性、长跑不蠕变、悬链线误差上界。
