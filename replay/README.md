# replay

动画回放一致性校验。`make test` 跑自检。

回放管线：关键帧曲线采样（同 `curve/`）→ 每骨骼局部姿态 → 骨骼层级世界变换（同 `skeleton/`），逐帧记录世界姿态。

## 用法

```
replay check  <anim> <ref.rec>   # 校验输入；同进程回放两次逐骨骼逐帧比对；
                                 # 再与手算参考姿态比对，给出误差上界
replay record <anim> <out.rec>   # 校验输入并回放，姿态记录写入 out.rec
replay diff   <a.rec> <b.rec>    # 比对两份记录，差异定位到骨骼与帧号
```

- 输入含 NaN/Inf、关键帧时间为负、帧起始时间为负：报错并以退出码 2 终止。
- 差异报告格式：`diff frame=<帧号> bone=<骨骼号> chan=<x|y|ang> ...`，
  并汇总 `max_abs_err` 及其位置；与参考比对时误差上界为 `REF_TOL=1e-4`。
- 样例文件（`samples/`、`tests/`）只以只读方式打开；唯一会写的文件是
  `record` 的输出路径。`make test` 运行前后对样例做 sha256 校验，确认未被改动。

## 文件格式

动画样例（`.anim`）：

```
bone <parent> <rest_x> <rest_y> <rest_ang>
key  <bone> <x|y|ang> <t> <v> <tan_in> <tan_out>
frames <count> <dt> <t0>
```

姿态记录 / 参考（`.rec` / `.ref`）：

```
rec <frames> <bones>
pose <frame> <bone> <x> <y> <ang>
```

`samples/chain3.ref` 为手算参考：根骨骼角度 `ang(t)=t`（埃尔米特切线均为 1，
退化为线性），子骨骼刚性跟随，世界姿态闭式解按双精度计算。
