# constraint

`make test` 跑自检。

## 调试导出

`make test` 里的自检会顺带开一次导出，产物落在 `out/debug/`（`make clean` 会删掉）：

| 文件 | 内容 |
|---|---|
| `residuals.csv` | 每帧每条约束一行：`frame,constraint,particle_a,particle_b,length,rest,residual` |
| `worst.csv` | 每帧一行：残差最大的约束下标及其残差 |
| `heatmap.ppm` | 残差热图（P6），横轴帧、纵轴约束；绿→黄→红对应残差 0→`DBG_RAMP_MAX`，超出即饱和 |

残差定义为 `|length - rest|`，与 `chain_seg_error()` 共用 `chain_residual()`，所以单条数值和最差统计不会各算各的。

导出由 `DebugExport` 控制，`dbg_open()` 失败时求解照常进行；求解器从不读导出结构，打开/关闭导出粒子状态逐位相同（自检里用粒子内存的 FNV-1a 哈希比对）。

`dbg_verify(dir)` 把三个文件读回来对账：`|length - rest|` 必须等于记录的残差，热图每个格子必须是单色且能按色标解回同样的残差，`worst.csv` 必须和逐约束行一致；任何一项对不上都会打印 `verify: ...` 并返回非零。
