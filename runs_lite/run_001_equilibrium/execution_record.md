# 执行记录 — Run 001

Level 0 Minimal TF / Equilibrium Characterization

---

## 代码版本

- **分支**: `learnability-level0-implementation`
- **Commit**: `58dcb78` — commit 6: add equilibrium characterization utilities
- **模块**: `src/ddc_lite.py`
- **日期**: 2026-05-18

---

## 系统配置

| 参数 | 值 |
|-----------|-------|
| 基因总数 (G) | 20 |
| TF 基因数 (N_TF) | 4（基因 0-3） |
| n（Hill cooperativity） | 2.0（固定） |
| gamma（翻译速率） | 1.0（固定） |

---

## 默认参数范围（来自 ddc_lite.sample_world）

| 参数 | 范围 |
|-----------|-------|
| rho（转录速率） | 0.5 ~ 2.0 |
| K（Hill 半饱和常数） | 1.0 ~ 5.0 |
| delta_x（mRNA 降解率） | 0.1 ~ 0.5 |
| delta_p（蛋白质降解率） | 0.05 ~ 0.3 |

---

## 模拟配置

| 参数 | 值 |
|-----------|-------|
| N_WORLDS | 20 |
| M_CELLS_PER_WORLD | 3 |
| BASE_SEED | 1000 |
| T_SIM | 200 |
| World seeds | 1000 ~ 1019 |
| Cell seeds | 1020 ~ 1079 |

---

## 轨迹配置

| 参数 | 值 |
|-----------|-------|
| 轨迹长度 (T_SIM) | 200 步 |
| 初始状态采样 | X0 ~ uniform(0, 1)；P0 = gamma × X0 |
| 调控模型 | Signed single-input Hill response |
| 更新顺序 | compute R → update X → update P |

---

## 收敛判定

| 参数 | 值 |
|-----------|-------|
| EPS（收敛阈值） | 0.01 |
| WINDOW（需连续满足的步数） | 10 |
| COLLAPSE_THR | 1e-3 |
| DIVERGENCE_THR | 1e3 |

**判定规则**: 若连续 `WINDOW` 步同时满足 `||X(t+w+1)-X(t+w)|| < EPS` 且 `||P(t+w+1)-P(t+w)|| < EPS`，则该 trajectory 判定为 converged。

---

## Per-world 聚合

每个 world 运行 3 个 cell。该 world 的判定标准：
- **converged**: n_converged >= 2（3 个 cell 中至少 2 个收敛）
- **bounded**: 全部 3 个 cell 均 bounded
- **collapsed**: n_collapsed >= 2

---

## 参数扫描

### sweep: default
- rho: 0.5 ~ 2.0
- delta_p: 0.05 ~ 0.3
- 其余参数保持默认范围

### sweep: rho0.2-0.8
- rho: 0.2 ~ 0.8
- delta_p: 0.05 ~ 0.3

### sweep: dp0.1-0.4
- rho: 0.5 ~ 2.0
- delta_p: 0.10 ~ 0.4

### sweep: rho0.2-0.8_dp0.1-0.4
- rho: 0.2 ~ 0.8
- delta_p: 0.10 ~ 0.4

---

## 结果汇总

| Sweep | Converged | Convergence Rate | Bounded Rate | Collapse Rate |
|-------|:---------:|:----------------:|:------------:|:-------------:|
| default | 6 / 20 | 30.0% | 100.0% | 0.0% |
| rho0.2-0.8 | 11 / 20 | 55.0% | 100.0% | 0.0% |
| dp0.1-0.4 | 7 / 20 | 35.0% | 100.0% | 0.0% |
| rho0.2-0.8_dp0.1-0.4 | 14 / 20 | 70.0% | 100.0% | 0.0% |

---

## 图拓扑

所有 sweep 共享同一张调控图（seed=1000）。

### TF 子图
```
TF0 ← TF2 (activation)
TF1 ← TF0 (activation)
TF2 ← TF1 (repression)
TF3 ← TF2 (repression)
```

三节点环路 `TF2 → TF0 → TF1 ⊣ TF2`（两正一负，奇数条负边构成负反馈环）是系统的振荡核心。

---

## 主要发现

1. 全部 20 个 world 在所有 sweep 中均为 **bounded**（无发散）。
2. 所有 sweep 中均无 **collapse**（趋近零表达）。
3. 仅降低 rho（0.2~0.8）即可将收敛率从 30% 提升至 55%。
4. 仅抬升 delta_p（0.1~0.4）在 T=200 时效果有限（35%），因为高降解需更长时间才能收敛。
5. **同时降低 rho + 抬升 delta_p** 在 T=200 时达到 70% 收敛率。
6. 剩余不收敛的 world（W1000, W1001, W1007, W1010, W1013, W1014）源于 TF 环路中关键位置 TF 基因的 delta_p 抽到过低值，配合非微小 rho 导致 feedback gain 仍然偏高。
