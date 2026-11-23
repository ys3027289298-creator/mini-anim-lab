# mini-anim-lab

角色动画基础模块，C 实现，每目录独立构建可自测。

| 目录 | 内容 |
|---|---|
| `common/` | 统一约束求解入口 `mal_solve()`、错误码与四个模块的公共头 |
| `skeleton/` | 骨骼层级、世界变换、骨段求取 |
| `ik/` | 双骨解析 IK、极向量、可达性判定 |
| `constraint/` | 距离约束链条、Verlet 积分、迭代收敛 |
| `path/` | 折线弧长表、按距离采样、曲率与最近点 |
| `curve/` | 关键帧采样、三次埃尔米特切线、夹取 |
| `pose/` | 四模块链接为一个二进制，跨模块姿态一致性与耗时 |
| `acceptance/` | 验收：基线逐字节对照、非法参数、姿态数值、两档耗时 |

## 统一求解入口

所有求解（骨骼更新/加骨、IK、曲线求值/加帧、约束链初始化/步进）都经过
`common/include/solver.h` 声明的 `mal_solve()`：

```c
mal_request r = { MAL_MODULE_IK, MAL_OP_IK_SOLVE, &req, &resp };
char why[128];
mal_status st = mal_solve(&r, why, sizeof why);  /* 失败时 why 以模块名开头 */
```

- 入口负责跨模块校验：空请求、模块 id 越界、空载荷。
- 模块负责领域校验：非法参数返回 `MAL_ERR_PARAM`，且错误信息以所属模块
  开头（如 `ik: bone lengths must be positive`），便于定位。
- 每个二进制只链接需要的模块，通过 `mal_dispatch_table` 注册；`pose/`
  链接全部四个模块并提供合并表。

### 关闭统一入口（回归原始实现）

模块编译开关 `MAL_UNIFIED`，默认 `1`。置 `0` 时自检直接调用原有 static
函数，数学代码两条路径共用，完全不经过 `mal_solve`：

```
make -C ik clean test MAL_UNIFIED=0
```

`acceptance/run_all.sh` 会自动把该路径输出与改动前基线逐字节 `diff`。

## 一键验收

```
acceptance/run_all.sh
```

包含五部分，缺一不可：

1. 四个模块各自独立自检，统一入口/原始路径两种模式；
2. 统一入口关闭时，四处 stdout 与 git HEAD 基线**逐字节一致**；
3. 23 个非法参数用例全部报错并指出所属模块（`acceptance/invalid`）；
4. 跨模块姿态一致性：同一目标与极向量下 IK、骨骼、约束链肘部/末端重合；
5. 每处改动前/改后姿态数值对照（`acceptance/compare_poses.sh`，原始代码
   取自 git HEAD）与两档耗时（单轮 + 批量）。

单模块独立验证：`acceptance/verify_module.sh <skeleton|ik|curve|constraint>`。
