# mini-anim-lab

角色动画基础模块，C 实现。

| 目录 | 内容 |
|---|---|
| `skeleton/` | 骨骼层级、世界变换、骨段求取 |
| `ik/` | 双骨解析 IK、极向量、可达性判定 |
| `constraint/` | 距离约束链条、Verlet 积分、迭代收敛 |
| `path/` | 折线弧长表、按距离采样、曲率与最近点 |
| `curve/` | 关键帧采样、三次埃尔米特切线、夹取 |

## 统一自检

```sh
make test                 # 五个模块，一次跑完
make test-ik              # 只跑一个模块
make per-module           # 原来的 make -C <dir> test
make report               # 重看上次记录的证据
```

`make test` 对每个模块跑五种构建：原有自检、姿态数值、两档耗时、图像哈希、
复现命令。任一模块的自检项缺失/失败，或五样证据缺一样，退出码非零，并在
`[FAIL] <模块>` 行后列出原因。结果写入 `artifacts/selfcheck.json`（机器读）和
`artifacts/selfcheck.txt`（人读）。

### 五样证据

| 证据 | 形式 | 判定 |
|---|---|---|
| 通过项计数 | 每模块 `name=1` 行，逐项对照基线 | 原有检查项一个不少、全为 1 |
| 姿态数值 | `mal_pose k=… v=…` 加聚合哈希 | 数值与哈希可跨次比对 |
| 耗时两档 | 同源码 `-O0` / `-O2`，各跑 3 次取中位 | 两档迭代数相同、均为正测量 |
| 图像哈希 | 128×128 线光栅的 CRC-32 | 每模块一个稳定哈希 |
| 复现命令 | 带源码 sha256 的重建命令 | runner 当场执行，须复现自检 |

各模块源码在 `-DMAL_EVIDENCE_POSE / TIMING / IMAGE` 下额外打印证据，证据不计入
检查项，因此任何构建模式下检查项列表一致，可与基线严格对照。共享脚手架在
`common/harness.h`，runner 在 `tools/selfcheck.py`。

基线检查项：skeleton 7、ik 8、curve 10、constraint 6、path 13，共 44 项。
