# constraint

`make test` 跑自检。

求解器拆分为 `src/chain.c` / `src/chain.h`，自测入口在 `src/constraint.c`，
拆分前后自检输出逐字节一致。`chain.h` 暴露默认 NULL 的只读观察钩子
（`ChainObserver`），供 `../audit` 链接同一份求解器代码做守恒审计；
不挂钩子时求解器行为与原单文件版本完全相同。
