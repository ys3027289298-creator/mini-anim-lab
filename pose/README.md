# pose

把 skeleton / ik / curve / constraint 四个模块链接进同一个二进制，验证
统一入口下的跨模块姿态一致性。

`make test` 运行的场景（全程只调用 `mal_solve()`）：

1. 曲线（curve）按时间生成移动目标点；
2. 解析 IK 求出双骨肘部世界坐标与两个关节角；
3. 用 IK 关节角搭建骨架（skeleton），世界变换得到的肘部必须与 IK 重合；
4. 同样两根骨长的 Verlet 距离约束链（constraint），锚点钉在根、末端钉在
   目标、中粒子受同一极向量方向的弯曲牵引，收敛后的肘部必须与 IK 重合。

输出 `d_ik_skel` / `d_ik_chain` / `d_end` 三个偏差均为 0，并给出两档耗时
（单轮 `time_one_round_us`、2000 轮总量与单轮均值）。

模块源码以 `-DMAL_NO_MAIN -DMAL_UNIFIED=1` 编译：不带各自 main，只贡献
包装函数与 dispatch；合并分发表由 `src/pose.c` 提供。
