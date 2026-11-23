# lighting

光照贴图(间接光+烘焙直接光)与实时直接光的混合,修复直接光重复计亮导致的漏光。

方案:烘焙时额外输出直接光分量图,运行时 `final = albedo * (max(0, lightmap - direct_lm) + rt_direct)`。

`make test` 跑自检,输出四项验收数据:

1. **漏光像素表**:8x8 块计数表 + `out/leak_before.csv` / `out/leak_after.csv` 逐像素清单
2. **混合前后误差**:相对参考图的 max/mean 误差,混合后上界 0.05
3. **图像哈希**:FNV-1a 64,黄金哈希锁定在源码中(`GOLDEN_NAIVE` / `GOLDEN_BLEND`)
4. **复现命令**:`make -C lighting test`;关闭混合 `./lighting --no-blend --out out`,输出与现状逐字节一致
