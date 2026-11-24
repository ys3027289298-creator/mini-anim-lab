# mini-anim-lab

角色动画基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `skeleton/` | 骨骼层级、世界变换、骨段求取 |
| `ik/` | 双骨解析 IK、极向量、可达性判定 |
| `constraint/` | 距离约束链条、Verlet 积分、迭代收敛 |
| `path/` | 折线弧长表、按距离采样、曲率与最近点 |
| `curve/` | 关键帧采样、三次埃尔米特切线、夹取 |

## 统一姿态导出格式

`skeleton`、`ik`、`curve`、`path` 四个模块共用一套导出格式（实现见
`shared/pose_format.h`），取代此前各写一套的写法。文件是行式文本：

```
mal-pose 1
kind <skeleton|ik|curve|path>
count <记录数>
<记录名> <字段>=<值> <字段>=<值> ...
```

`magic/version/kind/count` 与记录语法四处一致，只有记录名和字段集按模块不同：

| 模块 | 记录 | 字段 |
|---|---|---|
| `skeleton` | `bone`（每根骨一条） | `index parent lx ly lang wx wy ang` — 层级（父索引 + 局部变换）与世界变换 |
| `ik` | `ik`（单条） | `root_x root_y target_x target_y pole_x pole_y l1 l2 ok joint_x joint_y a1 a2` — 末端目标、关节位置与两个关节角 |
| `curve` | `key`（每个关键帧一条） | `t v tin tout` — 时间、值与入/出切线 |
| `path` | `meta` + `pt` | `total`；`x y cum` — 顶点与累计弧长（弧长表） |

浮点写成 `%.9g`（binary32 的 `FLT_DECIMAL_DIG`），读回按 float 解释后与原值
逐位相等，因此往返比对用的是精确相等而非容差。

### 命令

每个模块以参数选择动作，无参数时行为不变（见下节）：

```
./skeleton/skeleton --pose-export  out.pose   # 导出新格式
./skeleton/skeleton --pose-import  out.pose   # 读回新格式
./skeleton/skeleton --pose-roundtrip out.tmp  # 导出后读回，与导出前内部状态比对
./skeleton/skeleton --legacy-export out.csv   # 导出旧格式
./skeleton/skeleton --legacy-check  old.csv   # 读旧文件，与当前内部状态比对
```

旧格式（迁移前各模块自己的写法，均为无头的 CSV）读取路径保留，读入结果与
改动前一致：

| 模块 | 旧格式 |
|---|---|
| `skeleton` | `index,parent,lx,ly,lang,wx,wy,ang` 每行一根骨 |
| `ik` | 单行 `root_x,...,l2,ok,joint_x,joint_y,a1,a2` |
| `curve` | `t,v,tan_in,tan_out` 每个关键帧一行 |
| `path` | 首行 `total <total>`，其后每行 `x,y,cum` |

`tests/format_check.sh` 对四个模块分别做：构建、新格式往返、新格式样例比对、
旧样例读回、旧默认输出逐字节比对。样例在 `samples/`（新）与
`samples/legacy/`（旧），旧默认输出的基线在 `tests/samples/baseline/`。

## 自检

各模块 `make -C <module> test` 仍向 stdout 打印原有的 `name=value` 行，
内容与本次改动前逐字节一致。
