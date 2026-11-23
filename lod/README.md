# lod

骨骼 LOD 与视锥剔除联动模块。`make test` 跑自检（20 项断言：焊接零误差、跨界段
保留、锚点规则、整段离屏剔除、迟滞、无洞、旧管线缺陷复现、逐像素确定性等）。

## 两条管线

- `linked=1`（联动）：
  - 合并 LOD 末端焊接（weld）到当前姿态的 LOD0 末端效应器；
  - 合并胶囊半径按「弦高 + 中间关节骨半径」膨胀，直线段像素级覆盖弯曲轮廓；
  - 每帧按当前 LOD 重建剔除数据，线段对扩张视锥 AABB 求交（余量 0.12）；
  - 锚点规则：连接端在扩张视锥内则保留（肢体不脱节），否则整段剔除；
  - 切换距离 9.0 / 13.0，迟滞 0.8，边界来回不抖动。
- `linked=0`（旧管线，行为不变）：合并段沿首子骨方向保留编写长度、不焊接，
  剔除只做两端点判定。

## 离屏模拟

```
./lod sim <输出目录> <0|1> [帧数]
```

确定性、无随机源。输出 P6 PPM 画面、联系表、切换三联图与三张 CSV：
`frame_metrics.csv`、`error_table.csv`、`switch_log.csv`。

## 测试与复现

- `tests/debug_frame.c`：单帧几何/像素诊断（`-DLOD_NO_MAIN` 方式包含本模块）；
- `tests/ppm2png.py`：纯标准库 PPM→PNG（`python3 -I` 可运行）；
- `tests/reproduce.sh`：端到端验收复现，产物与四样交付物说明见
  `../artifacts/lod_acceptance/ACCEPTANCE.md`。
