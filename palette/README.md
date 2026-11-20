# palette

蒙皮矩阵调色板批处理：脏标记增量重算、双缓冲 + seqlock 发布、上传前完整性校验。

- 脏标记：只重算变化骨，未变项从 front 缓冲带到 back，统计实际重算项数
- 双缓冲：back 构建期间渲染线程读 front；每个缓冲带序号（偶=稳定，奇=写入中），
  校验通过后解封并原子交换
- 完整性校验：有限性/每槽代戳检查，再与全量重算逐元素 `memcmp`（位级零误差）；
  漏标返回 `PAL_ERR_DIRTY_MISS`，坏数据返回 `PAL_ERR_BAD_DATA`，拒绝发布并强制下一帧全量重建
- 撕裂检出：渲染线程两阶段读序号并复核校验和，半更新快照计入 `tears` 后重试，
  超过重试上限返回 `PAL_ERR_RETRY`
- 容量不足：加骨/输出缓冲过小均返回 `PAL_ERR_CAPACITY`

`make test` 跑确定性自检（含零误差比对、漏标与容量报错、半更新读取的检出与计数）。
`make thread-test` 跑 pthreads 生产/消费压力测试。
`make tsan-test` 在 ThreadSanitizer 下跑并发测试。

## 错误码与统计

| 错误码 | 触发条件 |
|---|---|
| `PAL_ERR_CAPACITY` | 加骨超过容量，或读取输出缓冲小于骨数 |
| `PAL_ERR_BAD_INDEX` | 脏标记/访问的骨索引越界 |
| `PAL_ERR_DIRTY_MISS` | 校验门发现某槽与全量重算不一致（脏标记漏标） |
| `PAL_ERR_BAD_DATA` | 非有限值（NaN/Inf）或槽从未写入 |
| `PAL_ERR_TORN` | 读取快照期间序号变化或校验和不符（半更新） |
| `PAL_ERR_RETRY` | 连续撕裂超过 `PAL_READ_RETRIES` 预算 |

校验失败时**不交换缓冲**：front 保持上一帧的安全数据，并置 `force_full`，
下一帧强制全量重建后恢复。所有关键计数（重算项、上传、撕裂、重试、校验失败）
均为原子操作，渲染线程与动画线程可并发读取。

## 关于 TSan 抑制

调色板是 seqlock：写者构建 back 缓冲时，渲染读者可能正在快照另一块缓冲，
序号 + 校验和会检出任何重叠并重试，撕裂快照从不返回。ThreadSanitizer 无法
建模 seqlock 协议，会把这种**故意重叠**的载荷访问报成 data race。代码已用
`__tsan_acquire/__tsan_release` 标注一致快照路径，剩余的撕裂态访问通过
`tsan-suppressions.txt` 精确豁免（仅限三个 seqlock 函数）；所有共享标量
（front 索引、序号、统计）都是原子访问。WSL2/新内核上 TSan 还需
`setarch --addr-no-randomize`，已内建在 `make tsan-test`。

