# constraint

`make test` 跑自检（`constraint` 与 `unified` 两个可执行）。

## unified：软硬约束统一求解

- 软约束 = 刚度项：`K` 由弹簧组装，`f` 为载荷。
- 硬约束 = 线性投影 `C d = b`，只在统一解上精确施加一次。
- 统一解 `d* = ds - K^-1 C^T S^-1 (C ds - b)`，其中 `S = C K^-1 C^T`、
  `ds = K^-1 f`。修正量由软解在硬流形上的残差驱动，硬约束不会把软块
  从其自身平衡位置拖走（`||K d* - f + C^T psi||` 在 `1e-14` 量级）。
- 软块按 `model_block` 的块划分独立 Cholesky 求解，并与稠密解逐分量比对；
  块划分切断刚度耦合、覆盖不全、含固定自由度都会拒绝。
- 对照解（分别求解再叠加）：`dn = ds + K^-1 C^T S^-1 b`。它只有在
  `C ds = 0` 时才满足硬约束，否则硬约束残差恰为 `C ds`；实测中为 `0.522`。
  两者的差 `delta = d* - dn` 满足
  `||delta||_K^2 = (C ds)^T S^-1 (C ds) <= ||C ds||^2 / lambda_min(S)`，
  同时给出更松的条件数上界 `||C ds|| sigma_max(C) sqrt(kappa(K)) / lambda_min(S)`。
- 报告 `cond(K)`、`cond(S)`；`K` 奇异或非正定、硬约束线性相关、刚度 `< 0`
  均拒绝并返回原因。
- `fixed` 自由度从系统中消去（位移保持 0）；所有范数均按自由自由度计。
