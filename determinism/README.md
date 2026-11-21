# determinism

绳链求解的确定性校验。同一份输入跑两次，形状必须逐位一致。

`make test` 跑自检。

## 为什么形状会不一样

绳子是 PBD（位置约束）链：Verlet 积分 + Gauss–Seidel 逐段投影。只要下面任意一条
不稳定，两次运行的最后一个 bit 就会不同，随后 200 帧的浮点误差被放大成肉眼可见
的形状差异：

1. **迭代顺序**：Gauss–Seidel 的修正写回共享状态，从 i=0 扫到 i=n-2 和从
   i=n-2 扫到 i=0 是两种不同的离散格式，结果不同；并行扫描时写回顺序取决于线程
   调度，同一程序两次也不一样。
2. **求和/归约顺序**：`(a+b)+c != a+(b+c)`。并行归约的 tree 归约与串行从左到右
   归约不等价，阈值切在哪里还会随线程数变化。
3. **并行归约本身**：跨线程的 `+=` 是非结合浮点累加，线程完成顺序一变结果就变。

本模块把这三条都钉死，并把「违反」变成可检测的错误，而不是随机形状。

## 求解器约定（src/determinism.c）

- 坐标用 `double`；所有三角函数走 libm 的同一条代码路径，不依赖 fast-math。
- **固定迭代顺序**：`sweep_forward`（i = 0..n-2）与 `sweep_backward`
  （i = n-2..0）是两种显式模式；`iterations` 次扫描按 `fwd,bwd,fwd,...` 交替，
  与 `iterations` 奇偶无关地可复现。
- **固定求和顺序**：每段约束只改两个端点，天然没有归约；诊断量
  （最大长度误差、质心）一律按索引从左到右累加，写成 `sum += x[i]` 的循环，
  不做 tree 归约、不排序、不动点迭代的退出判据也用 `>`（严格）避免并列。
- **禁用并行归约**：`parallelism` 参数只允许 1。>1 时求解器直接返回
  错误码 `DET_ERR_PARALLEL`，不做「尽力而为」的并行，所以不存在调度随机性。
- **禁止线程调度引入的不确定**：整个模块单线程，自检里用 `DET_THREADS=1`
  强制断言（见下）。

## API

```c
DetConfig cfg = { .particles = 8, .seg = 1.0, .iterations = 8,
                  .substeps = 1, .parallelism = 1, .dt = 0.016 };
DetRope a, b;
det_rope_init(&a, &cfg);          // 非法配置返回 -1
det_rope_init(&b, &cfg);
det_rope_reset(&a);               // 回到 t=0 初始构型（逐位可复现）
for (int f = 0; f < 240; f++) {
    det_rope_step(&a, 0.0, -9.8);
    det_rope_step(&b, 0.0, -9.8);
}
DetDiff d = det_rope_compare(&a, &b, cfg.particles);
d.bit_differences;  // 必须为 0
```

`det_rope_snapshot` 把每帧的 `x,y` 以原始 IEEE-754 bit 形式拷出来，配合
`det_snapshot_write/read` 可以跨进程逐帧比对（自检里就是两个独立进程比对）。

## 顺序打乱必须报错

`det_rope_init` 时记录的 `order_id` 是 `(fwd/bwd 顺序, iterations, substeps)`
的规范化指纹；`det_rope_step` 每次先校验 `order_id` 与当前扫描计划一致，
不一致返回 `DET_ERR_ORDER` 且**不推进状态**。自检用
`det_rope_force_order` 模拟错误的扫描顺序，断言 `step` 返回 `DET_ERR_ORDER`，
并且此时两次运行的位置依然逐位为零差异（因为被拒绝的那次没有写入任何状态）。

注意一个刻意保留的语义：**注入的顺序变化本身**（`DET_ORDER_FWD` vs
`DET_ORDER_ALT`）不是错误——它是两个合法的确定性配置。自检
`order_change_quantified` 量化它造成的差异，并断言该差异在重复运行时逐位相同；
只有「运行中途把计划改掉、又让它继续跑」才报错。

## 并行度差异的量化

并行度 >1 被拒绝，但 `parallelism` 作为**诊断维度**仍然可以比较：两种配置
各自串行运行，`det_rope_quantify` 报告

```
max_abs   = max_i |a[i] - b[i]|        （坐标差的最大绝对值）
l2        = sqrt(sum (a-b)^2)          （固定左到右求和）
changed   = 坐标 bit 不同的个数
```

`parallelism` 从 1 改到 4 时，因为 4 被拒绝、求解器不运行，差异必须是 **0**，
而不是随线程数抖动的随机值——自检断言 `max_abs == 0.0 && changed == 0`。如果将来
真的实现了并行路径，这两个量就是量化指标，且必须能解释来源（比如 tree 归约与串行
归约的差）；不允许「同一配置跑两次得到不同差异」。

## 自检项

`make test` 输出逐项 `name=1`，最后 `all_ok=1`。覆盖：

- `same_run_bitwise`：同配置双实例 240 帧，逐帧逐位 0 差异；
- `snapshot_roundtrip`：快照写盘再读回，bit 完全一致；
- `cross_process_bitwise`：`--trace` 子进程各跑一次，两个独立进程的输出
  与本进程逐帧逐位一致；
- `order_violation_errors`：打乱扫描顺序 `step` 返回 `DET_ERR_ORDER` 且状态不变；
- `parallelism_zero_diff`：parallelism 1 vs 4，差异量化为 0；
- `iterations_order_permutation`：不同 `iterations` 不是错误，但同配置下结果
  仍然是确定的（与逐位比较一致）；
- `illegal_config`：n<2、非有限段长、dt<=0、parallelism 0 或 >1 均返回 -1 / 错误码。
