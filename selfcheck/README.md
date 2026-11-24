# 统一自检入口 / unified self-check

五个模块（`skeleton` `ik` `curve` `constraint` `path`）的自检合成一个入口：

```sh
./selfcheck/run.sh
```

- 用各模块**自己的** Makefile 构建、**原样**运行各自的测试二进制，只解析 stdout
  里的 `name=1` / `name=0` 行。不改模块源码、Makefile，也不改二进制行为。
- 每条结果都带模块前缀，如 `[skeleton] world_chain=1`，可看出原来属于哪一处。
- 期望的自检项从各模块源码里的 `check("...")` 调用读出；某条没跑出来（崩溃、
  提前返回）会报 `MISSING`，不会被悄悄漏掉。
- 末尾打印每处 pass/fail/missing 和总计；任何失败、缺失、构建失败都会以非零退出。
- `make -C selfcheck` 是 `./selfcheck/run.sh` 的等价入口。

产物留档在 `artifacts/selfcheck/`。
