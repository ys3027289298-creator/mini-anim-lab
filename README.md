# mini-anim-lab

角色动画基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `skeleton/` | 骨骼层级、世界变换、骨段求取 |
| `ik/` | 双骨解析 IK、极向量、可达性判定 |
| `constraint/` | 距离约束链条、Verlet 积分、迭代收敛 |
| `path/` | 折线弧长表、按距离采样、曲率与最近点 |
| `curve/` | 关键帧采样、三次埃尔米特切线、夹取 |

## 自检（旧导出格式，保持不变）

```
make -C skeleton test   # 每个模块同样适用：ik / curve / path / constraint
```

无参数运行时各模块向 stdout 打印 `name=value` 自检行，该输出是旧导出格式，
其逐字节基线保存在 `tests/golden/<module>.out`，不得变更。

## 统一姿态导出（新格式，可选开关）

`skeleton`、`ik`、`curve`、`path` 四个模块支持同一套新姿态导出格式，
共享实现见 `shared/pose_export.h`，通过 `--pose-export <outdir>` 开启：

```
./skeleton/skeleton --pose-export artifacts/poses/skeleton
```

每次导出在 `<outdir>` 下生成两份文件：

| 文件 | 内容 |
|---|---|
| `pose.json` | 姿态数值、图像哈希、两档耗时、复现命令、格式版本 |
| `pose.ppm` | 姿态的 P6 PPM 渲染图（320x240，确定性输出） |

`pose.json` 的四项必备数据：

1. `pose` — 姿态数值（各模块字段不同：骨骼链 / IK 关节角 / 曲线采样 / 路径采样）
2. `image.sha256` — `pose.ppm` 的 SHA-256，C 内置实现，验收脚本会外部复算核对
3. `timing` — 两档耗时：`compute_pose_us`（1000 次姿态计算，微秒）
   与 `render_image_us`（单次图像渲染，微秒）
4. `reproduce` — 实际调用的完整命令行，可原样重跑复现同一图像

不传 `--pose-export` 时走旧自检路径，stdout 与基线逐字节一致。

## 统一验收

```
python3 -I tests/acceptance.py
```

对四个模块逐一验证：旧输出逐字节一致、四项数据齐全、哈希外部核对通过、
两次导出图像与姿态值完全一致（确定性）、复现命令可重跑且图像哈希相同。
四份数据缺一即判未完成，脚本以非零码退出。

`tests/render_contact.py` 用标准库把四张 PPM 拼成
`artifacts/contact_sheet.png` 供目视检查：

```
python3 -I tests/render_contact.py artifacts/poses/*/pose.ppm
```
