# skeleton

`make test` 跑自检（旧导出格式，stdout 基线见 `tests/golden/skeleton.out`）。

新姿态导出格式（详见根目录 README）：

```
./skeleton/skeleton --pose-export <outdir>
```

生成 `<outdir>/pose.json`（姿态数值、图像哈希、两档耗时、复现命令）
与 `<outdir>/pose.ppm`。不传参数时旧自检输出逐字节不变。
