# acceptance

统一入口的验收套件。入口脚本：

- `run_all.sh`：五部分全量验收，结果存 `results/`。
- `verify_module.sh <module>`：单模块独立验证（两种调用路径 + 基线 diff）。
- `compare_poses.sh`：改动前/改后姿态数值逐项对照。

目录说明：

| 路径 | 作用 |
|---|---|
| `baseline_out/` | 改动前（git HEAD）四个模块的 stdout/stderr 基线，入库保存 |
| `reference/orig/` | 从 git HEAD 导出的未改动模块源码 |
| `reference/ref_*.c` | 直接 include 原始源码、调用其 static 函数的参考驱动，输出改前姿态 |
| `src/invalid.c` | 23 个非法参数用例，校验报错状态码与模块前缀 |
| `bench/bench_*.c` | 每模块姿态快照 + 统一/直连对照 + 两档耗时（单轮、批量） |
| `legacy_out/`、`unified_out/`、`results/` | 运行生成物，已 gitignore |

## 验收项与需求对应

1. **每处独立可跑**：四个模块目录 `make test` 依旧独立；`verify_module.sh`
   分别验证。
2. **关闭统一入口逐字节不变**：`MAL_UNIFIED=0` 构建产物的 stdout 与
   `baseline_out/*.stdout.txt` 用 `diff` 比对，必须 IDENTICAL。
3. **非法参数报错并指出模块**：`invalid` 每个用例打印
   `<module>.<case> status=<n> msg="<module>: ..." ok=1`。
4. **改前/改后姿态数值对照**：`compare_poses.sh` 从参考驱动取改前值、从
   bench 取改后值，逐值匹配；constraint 额外用 `memcmp` 断言粒子块
   逐字节一致。
5. **两档耗时**：bench 同时输出单次调用耗时与大批量（1 万~10 万次）
   总耗时/单次均值。
