# mini-anim-lab

角色动画基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `skeleton/` | 骨骼层级、世界变换、骨段求取 |
| `ik/` | 双骨解析 IK、极向量、可达性判定 |
| `constraint/` | 距离约束链条、Verlet 积分、迭代收敛 |
| `path/` | 折线弧长表、按距离采样、曲率与最近点 |
| `curve/` | 关键帧采样、三次埃尔米特切线、夹取 |

## 自检与可见产物

每个模块的 `make test` 除保留原有全部检查项外，还会把**可见产物**写入
`artifacts/run_<module>/`，并对落盘后的文件本身（而非内存状态）做断言：

| 模块 | 图像产物 (PGM) | 数值产物 (CSV) |
|---|---|---|
| `skeleton/` | `skeleton.pgm` 骨架栅格图 | `skeleton_pose.csv` 世界变换/骨段 |
| `ik/` | `ik_pose.pgm` 双骨姿态图 | `ik_pose.csv` 关节点/角度/可达性 |
| `constraint/` | `constraint.pgm` 链条图 | `constraint_chain.csv` 粒子坐标/约束误差 |
| `path/` | `path_length_map.pgm` 弧长场图 | `path_samples.csv` 等距采样 |
| `curve/` | `curve.pgm` 曲线图 | `curve_samples.csv` 采样值/夹取区 |

产物目录可用环境变量覆盖（默认 `artifacts/run_<module>/`）：
`ARTIFACT_DIR=/tmp/out make -C skeleton test`。

## 缺陷注入测试台

```
python3 -I tests/selftest.py
```

在隔离临时副本中干净构建，输出五块证据：基线通过计数、产物 sha256 哈希、
每模块 6 个注入缺陷的捕获率表、低于 80% 时的漏网清单、复现命令。
缺陷按构建失败 / 检查失败 / 产物哈希漂移 / 产物缺失四类判定捕获。

