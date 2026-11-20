# retarget

镜像重定向：把一套姿势映射到左右对称的另一侧。

- `fit_symmetry_axis`：对标记为矢状面（root/脊柱/颈/头）的骨骼做正交最小二乘直线拟合，取散布矩阵的次特征向量为镜面法向；残差为这些点到直线的 RMS 距离，超过容忍值时返回 `MIRROR_WARN` 供上层报告。
- `plane_reflect`：`p' = p - 2((p-o)·n)n`，平移只翻法向分量，镜两次回到原处。
- `conj_angle`：旋转用 `M R M` 共轭，2D 下等价于局部角度取反（符号修正）。
- `mirror_pose`：同一骨架的几何镜像，严格对合（`mirror_twice_identity`）。
- `mirror_retarget`：把镜像后的世界旋转写到 partner 骨骼上，目标骨架保留自身骨长；与手工镜像姿势（`manual_mirror`）逐骨比对一致；末端误差以两侧骨长差之和为上界。

`make test` 跑自检。
