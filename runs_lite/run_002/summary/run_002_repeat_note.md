# run_002 — Re-run Note

本目录是修复 `ddc_model_a.py` 中 `simulate_single_cell` 干预分支的记录 bug 后重新运行的结果。

## 变更说明

- **时间线**：
  - 旧版运行时间：2026-06-03（原始 run_002，已移至 `run_002_deprecated`）
  - 代码修复时间：2026-06-11
  - 新版重新运行时间：2026-06-11（即本目录）

- **Bug**: 旧版 `ddc_model_a.py` 中 `simulate_single_cell` 的干预分支缺少 `X_traj[t + 1] = X` 赋值语句。该分支先调用 `apply_intervention` 修改当前状态 X，然后计算 `A @ X + b` 并 clamp 得到下一步状态 `X_next`，但忘记将 `X_next` 写入轨迹数组对应位置。后果：

  - `X_traj[501]`（干预后第一步，t=500 的下一状态）保留预分配 tensor 的初始零值，未被覆盖
  - 后续主循环从干预步的 `continue` 之后继续，重复执行 `A @ X + b` → clamp → `X_traj[t+1] = X`，因此 `X_traj[502]` 实际记录的是本应在 `[501]` 位置的值，`[503]` 记录的是 `[502]` 的值，以此类推
  - 即：`old[t] == correct[t-1]` 对 t ≥ 502 全部成立 —— 纯属**时序记录偏移 1 位**，模拟的动力学过程本身（`A @ X + b` → clamp 的计算结果）完全相同
  - 视觉影响：绘图时在 x=501 处出现一个 1 步宽的 V 形尖刺（所有基因瞬时跌至 y=0 后反弹），为明显的人工痕迹

  **旧版（buggy）代码**：
  ```python
  if intervention_time is not None and t == intervention_time:
      X = apply_intervention(X, intervention_config)
      X_traj[t] = X
      clip_count[t] = 0
      # ... intervention_history ...
      raw = A @ X + b
      X_next = torch.clamp(raw, min=0.0)
      X = X_next
      # ❌ 漏写: X_traj[t + 1] = X
      # ❌ 漏写: clip_count[t + 1] = ...
      continue
  ```

  **新版（修复后）代码**：
  ```python
  if intervention_time is not None and t == intervention_time:
      X = apply_intervention(X, intervention_config)
      X_traj[t] = X
      clip_count[t] = 0
      # ... intervention_history ...
      raw = A @ X + b
      X_next = torch.clamp(raw, min=0.0)
      X = X_next
      X_traj[t + 1] = X              # ✅ 补上
      clip_count[t + 1] = (raw != X).sum().item()  # ✅ 补上
      continue
  ```
  差异仅在于缺失的两行赋值（记录层），干预状态 `X` 和下一步的计算 `A @ X + b → clamp` 完全一致。

- **修复**: 在干预分支中补上 `X_traj[t + 1] = X` 一行。修复后 `X_traj` 所有时间步各归其位，与干预前轨迹 t=0..500 完全无缝衔接。
- **旧版备份**: `../run_002_deprecated/`（原始 run_002 已删除，仅保留本地备份）。

## 与 run_002_deprecated 的一致性

| 数据 | 状态 |
|---|---|
| trajectories（.json） | 完全一致 |
| world_metadata（.json） | 完全一致 |
| all_world_results.json | 完全一致 |
| topo_strength_analysis.json | 完全一致 |
| perturbation .pt 文件 | 差 1 步索引偏移（`old[t] == new[t-1]` for t≥502） |
| perturbation_summary.json | recovery_time 等字段差 ±1，最终分类一致 |

## 新增分析

- out-degree distribution vs stability（柱状图 + 拓扑图）
- 所有 in/out-degree 函数和输出文件统一命名
