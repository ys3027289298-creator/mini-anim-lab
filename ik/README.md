# ik

`make test` 跑自检（默认经统一入口 `mal_solve()`）。
`make test MAL_UNIFIED=0` 走原始 static 函数路径，输出与改动前逐字节一致。
独立验证：`../acceptance/verify_module.sh ik`。
