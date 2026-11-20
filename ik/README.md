# ik

双骨解析 IK，外加任意节数链条的 FABRIK 迭代求解。`make test` 跑自检。

## API

```c
#define IK_MAX_BONES 16

typedef struct {
    int iterative;        /* 0 = 旧路径：两骨走解析解，多节链不动 */
    int max_iterations;   /* 迭代开启时必须 > 0 */
    float tolerance;      /* 迭代开启时必须 > 0 */
} IkOptions;

typedef struct {
    int bone_count;                  /* 1..IK_MAX_BONES */
    V2 root;
    float length[IK_MAX_BONES];      /* 每节骨长，必须 > 0 */
    V2 pole;                         /* 仅解析路径用：弯曲方向 */
    V2 pose[IK_MAX_BONES + 1];       /* 迭代种子姿态；pose[0] 会被 root 覆盖 */
} IkChain;

typedef struct {
    IkStatus status;
    int iterations_used;
    float residual;                  /* 最后一次迭代的最大关节位移 */
    int used_analytic;
    V2 joints[IK_MAX_BONES + 1];     /* 解出的 bone_count+1 个点，首点为根 */
} IkSolveResult;

IkSolveResult ik_solve(const IkChain* chain, V2 target, IkOptions opt);
IkOptions ik_default_options(void);  /* iterative=1, max_iterations=64, tolerance=1e-5 */
```

## 迭代算法

FABRIK：前向把末端钉到目标并沿链回推，后向把根钉住并沿链外推。收敛判据是**一次完整迭代内所有关节的最大位移**（`residual`），`residual <= tolerance` 即收敛。超过 `max_iterations` 仍不满足则返回 `IK_ERR_NOT_CONVERGED`，不会把未收敛的姿态当成功返回。

不可达目标（距离 > 骨长总和，或距离 < `2*max_len - total`）会让链条朝目标伸直/折叠，最终稳定在极限姿态，仍可正常收敛并返回 `IK_OK`；此时末端与目标存在物理上无法消除的残差。

## 状态码

| 值 | 名称 | 含义 |
|---|---|---|
| 0 | `IK_OK` | 成功 |
| 1 | `IK_ERR_BAD_CHAIN` | 节数越界或输入含 NaN/Inf |
| 2 | `IK_ERR_BAD_LENGTH` | 骨长 <= 0 或非有限 |
| 3 | `IK_ERR_BAD_TOLERANCE` | `tolerance <= 0` 或非有限 |
| 4 | `IK_ERR_BAD_ITERATIONS` | `max_iterations <= 0` |
| 5 | `IK_ERR_ITERATION_DISABLED` | 关闭迭代且链长 > 2，姿态原样返回 |
| 6 | `IK_ERR_NOT_CONVERGED` | 达到迭代上限仍在容差外 |
| 7 | `IK_ERR_NO_SOLUTION` | 解析路径无解（不可达、过近、零长） |

## 关闭迭代

`iterative = 0` 时：

- 两根骨：走原来的 `solve_two_bone`，输出的关节、末端姿态与改动前逐位相同（自检 `off2_pose_identical`）。
- 三节及以上：不做任何求解，`joints` 原样拷贝自 `pose`，返回 `IK_ERR_ITERATION_DISABLED`。

## 与解析双骨解对比

对 `l1=1.0, l2=0.8`，在 16 个方向 × 5 个可达比例（0.30/0.50/0.70/0.90/0.98）共 80 组用例上，比较 FABRIK 收敛解与解析解的关节位置与末端位置，取较大者：

```
two_bone_cases=80
two_bone_converged=80
two_bone_err_max=5.11776261e-06      (tolerance = 1e-6)
two_bone_err_bound=9.99999975e-06    (10 × tolerance)
```

**位置误差上界**：对严格可达的目标，迭代解与解析解的误差 ≤ `10 × tolerance`（实测 4.98–5.61 倍）。该倍数作为 `TWO_BONE_ERR_RATIO` 断言在自检里。

## 迭代次数与残差

`./ik --curve` 输出四节链（骨长全 1，目标 `(1.0, 2.4)`）每次迭代后的残差，CSV 格式，可直接画曲线：

```
iteration,residual
1,2.01151371
2,0.0419625826
3,9.27145593e-05
4,1.1920929e-07
5,1.1920929e-07
6,0
```

可配置：`./ik --curve --curve-tol 1e-3 --curve-iters 200`。容差参数非法会以退出码 2 拒绝。

容差与所需迭代数、最坏误差的关系（同一双骨网格）：

```
tol,max_iterations,converged,cases,max_err,err_over_tol
1e-03,17,80,80,0.00497543672,4.98
1e-04,30,80,80,0.000545660674,5.46
1e-05,44,80,80,5.23831659e-05,5.24
1e-06,58,80,80,5.11776261e-06,5.12
1e-07,74,73,80,7.11676705e-07,7.12
```

`1e-7` 一档触及 float32 精度下限：近伸直构型下残差会停在约 7e-7 不再下降，7 组用例如实返回 `IK_ERR_NOT_CONVERGED`。自检对低于 `1e-6` 的容差只要求“收敛的仍在误差界内、失败的必须显式报未收敛”，不会静默放过。

## 已知限制

- FABRIK 在近乎伸直的奇异区收敛慢。随机扫描 2–8 节链时，`距离/总长` 接近 1 的目标需要更多迭代：把 `max_iterations` 从 256 提到 2000 后，绝大多数此类用例能在 ~1300 次内达到 1e-5 容差；极个别（如 8 节链、`距离/总长 ≈ 0.9998`）2000 次仍未达到，按设计如实返回 `IK_ERR_NOT_CONVERGED`。
- 可达但接近“极限伸展”的目标，末端误差与 `residual` 不是同一个量：`residual` 是最后一次迭代的关节位移，会先于末端误差降到容差以下。双骨网格上两者最大约 5.6 倍容差，长链上可达 `1e-4` 量级。需要严格末端精度时应把 `tolerance` 设得比目标误差小一个量级。
- 解析路径的 `pole`（弯曲分支）在迭代路径中没有对应概念：FABRIK 的弯曲分支由种子姿态 `pose` 决定。
