# ik

`make test` 跑自检。

## API

- `solve_two_bone(root, target, l1, l2, pole)` — 原始解析解。不可达 / 过近 / 零长骨 / 零距离时
  返回 `ok = 0`，行为与修改前一致。
- `reach_level(root, target, l1, l2)` — 解析可达性判定，返回 `IK_REACHABLE` /
  `IK_TOO_CLOSE`（内环以内）/ `IK_UNREACHABLE`（外环以外），阈值与 `solve_two_bone` 共用。
- `ik_chain_init(&chain, root, l1, l2, pole)` — 初始化策略层。
- `ik_chain_solve(&chain, request, dt)` — 每帧求解。请求在环外时先钳制到可达环上，再按
  `speed * dt` 沿环平滑走过去；请求可达时直接交给原始解，姿势与修改前逐位一致。
- `ik_chain_status(&chain)` — 查询上一帧结果（按值返回，链为空也安全）。
- `ik_chain_target(&chain, &out)` — 查询上一帧实际使用的目标。
- `ik_chain_set_smooth_speed(&chain, speed)` — 过渡速度（单位/秒），仅接受正数。

## IkStatus

| 字段 | 含义 |
|---|---|
| `ok` | 本帧是否产出姿势；`err` 非零时为 0 |
| `status` | 请求相对可达环的位置（`IkReach`） |
| `err` | `IK_ERR_DEGENERATE_LIMB` 零长骨 / `IK_ERR_DEGENERATE_TARGET` 目标与根重合 |
| `clamped` | 请求被钳制过 |
| `smoothing` | 钳制后的目标仍在过渡中 |
| `t` | 当前过渡进度 0..1（1 = 已稳定） |
| `dist` | `‖request - root‖` |
| `request` / `applied` | 原始请求 / 本帧实际求解用的目标 |

`err` 优先于 `ok` 与 `status`：零长骨时不会给出任何姿势，也不会做钳制。

## 行为约定

- 可达请求：`applied == request`，姿势与 `solve_two_bone` 逐位相同。
- 不可达请求（环外或环内洞）：钳制点保持请求方向，落在环上；过渡每帧位移不超过
  `speed * dt`，路径在极坐标下插值，半径始终位于 `[|l1-l2|, l1+l2]`，不会穿过内洞。
- 目标稳定在不可达位置时，过渡结束后 `smoothing = 0` 且目标不再移动（无抖动）。
- `dt <= 0`（含 NaN）按标称 1/60 帧处理，`dt` 上限 0.1 秒。
