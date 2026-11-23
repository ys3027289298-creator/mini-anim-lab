# common

统一约束求解层。

- `include/solver.h`：`mal_solve()` 入口、模块/操作枚举、状态码
  （`MAL_OK` / `MAL_ERR_PARAM` / `MAL_ERR_UNSUPPORTED` / `MAL_ERR_STATE`）、
  请求结构与分发表 `mal_dispatch_entry`。
- `include/mal_err.h`：模块前缀错误信息工具 `mal_fail()` 与有限数检查。
- `include/{skeleton,ik,curve,constraint}.h`：各模块公共类型与包装函数。
- `src/solver.c`：入口实现，只做跨模块校验与分发表查询，不含领域数学。

独立模块二进制各自定义单条目的 `mal_dispatch_table`；`pose/` 与
`acceptance/` 下的整合程序定义四条目合并表。
