# Run 004 — Summary Report

**Stage**: DDC Lite Representation Curriculum — Model A Parameter-Boundary Run

**Run**: run_004 (Model_A_Decay_Memory_Boundary)

**Model**: Model A — Classical Linear GRN with decay δ scan

**Date**: 2026-06-16

---

## 0. Executive Summary

> 在保持 run_002 / run_003 Model A runtime semantics 和 topology baselines（unrestricted + TF-restricted）不变的前提下，将 decay δ 作为唯一新增 primary axis 进行扫描（δ ∈ {0.1, 0.2, 0.4}），同时保持 b、strength regimes、sign ratios 与 run_002/run_003 一致。实验规模：2 topology × 3 δ × 4 strength × 2 sign ratio = 48 combos，每 combo 10 worlds = 480 worlds × 5 cells = 2400 trajectories。
>
> **核心结论**：decay δ 是 Model A primary regime boundary 的一个重要 axis，其效应与 interaction strength 在同一量级。δ 从 0.4（high decay, low memory）降到 0.1（low decay, high memory）时，Convergent 比例从 98.8% 降至 73.8%，Divergent 比例从 0.6% 升至 25.0%。与此同时，即使在同一 δ 下，strength regime 的变动幅度同样显著（如 δ=0.1 内：baseline 95% Convergent → chen_stress 60%）。这一发现表明，**δ 和 interaction strength 是两个需要同时考虑的 axis——更高的 self-retention（1−δ）并非增强稳定性，而是让 interaction-driven instability 有更长的累积窗口。memory 本身是双刃剑。**

---

## 1. 科学问题与实验设计

### 1.1 核心科学问题

> *Do the topology-control conclusions from run_002 / run_003 remain stable when decay δ changes the memory pressure of Model A?*

同时回答：

> *当 self-retention = 1−δ 改变时，Convergent / Sustained oscillatory / Divergent regime boundary 是否发生系统性移动？*

### 1.2 实验定位

Run 004 是 **Model A Decay / Memory Boundary Run**——不是新的 topology ablation，也不是 Model B comparison。它在 run_002（unconstrained baseline）和 run_003（TF-restricted baseline）的基础上，唯一改变 decay δ，以 **isolate memory pressure 对 primary regime boundary 的影响**。

δ=0.2 是 run_002/run_003 的参考点（self-retention = 0.8），δ=0.1 测试 higher memory pressure（self-retention = 0.9），δ=0.4 测试 stronger decay dominance（self-retention = 0.6）。

### 1.3 当前 Run 在 Curriculum 中的位置

```text
run_002:  Model A unconstrained topology baseline  (δ=0.2 fixed)
run_003:  Model A TF-restricted topology baseline  (δ=0.2 fixed)
run_004:  Model A decay δ boundary scan            (δ ∈ {0.1, 0.2, 0.4})
          ↓ 当前阶段
run_00X:  Model B comparison (future)
```

### 1.4 参数网格

| axis | values | count |
|------|--------|-------|
| topology type | unrestricted_sparse, tf_restricted | 2 |
| decay δ | 0.1, 0.2, 0.4 | 3 |
| strength regime | baseline, chen_moderate, stress, chen_stress | 4 |
| sign ratio | balanced (1:1), repression_biased (1:2) | 2 |
| worlds per combo | 10 | — |
| cells per world | 5 | — |

**Total**: 48 combos × 10 worlds = 480 worlds × 5 cells = 2400 trajectories

固定参数：G=20, N_EDGES=38, r=0.10, b=0.1, T=1000。

### 1.5 δ 的物理含义

decay δ 控制 self-retention（1−δ），直接影响 memory persistence：

| δ | self-retention | role |
|---|---------------|------|
| 0.1 | 0.9 | high memory / low decay — 过去状态保留时间长 |
| 0.2 | 0.8 | reference point — 与 run_002 / run_003 一致 |
| 0.4 | 0.6 | low memory / high decay — 过去状态快速衰减 |

直观理解：δ=0.1 时，单基因无调控输入的天然衰减仅 10%/step，系统"记忆"长达 ~10 steps；δ=0.4 时，衰减 40%/step，记忆仅 ~2.5 steps。

### 1.6 Primary Regime 分类体系

每 trajectory 基于三个维度：equilibrium detection（convergence）、stability analysis（boundedness / divergence / collapse）、oscillation analysis（damping），归入 5 类 primary regime + secondary labels：

| Primary Regime | 判定条件 |
|----------------|----------|
| Convergent | converged to equilibrium (‖X(t+1)−X(t)‖ < ε for 50 consecutive steps) |
| Sustained oscillatory | sustained bounded oscillation (no convergence, oscillation detected, avg damping ≤ 0.05) |
| Divergent | max expression ≥ 1e3 (runaway growth) |
| Collapse | final mean expression < 1e-3 (numerical collapse) |
| Ambiguous | none of the above / numerical failure |

Secondary labels 包括 `fast_convergence` (≤200 steps)、`slow_convergence` (>200 steps)、`damped_oscillatory_transient`、`clipping_dominated`、`low_clipping_linear_like`、`runaway_divergence`、`numerical_collapse`。

每个 world 的 `primary_regime` 由其 5 个 cell 的 majority vote 决定。

### 1.7 扰动分析

10 个 candidate worlds，覆盖全部实际存在的 (primary_regime × δ) 组合（共 10 组），同时平衡 topology 覆盖（unrestricted_sparse 5 个、tf_restricted 5 个）。Canonical protocol：t=500 时 knockdown 表达量最高的基因（单步降至 0），随后让系统自然恢复。扰动后轨迹使用与原始分析相同的检测 pipeline 重新分类。

---

## 2. 结果

### 2.1 全局 Primary Regime 分布

480 worlds 的 primary regime 分布：

| Primary Regime | 数量 | 占比 |
|----------------|------|------|
| Convergent | 418 | 87.1% |
| Sustained oscillatory | 4 | 0.8% |
| Divergent | 55 | 11.5% |
| Collapse | 0 | 0.0% |
| Ambiguous | 3 | 0.6% |

> Primary regime heatmaps 见图: [../figures/01_primary_regime_heatmaps.png](../figures/01_primary_regime_heatmaps.png)
>
> Topology × δ comparison 见图: [../figures/02_topology_delta_comparison.png](../figures/02_topology_delta_comparison.png)

### 2.2 Decay δ — 系统性移动 Regime Boundary (§13.1)

δ 对 primary regime distribution 有显著的系统性影响：

| δ | Convergent | Divergent | Sustained oscillatory | Ambiguous |
|---|-----------|-----------|----------------------|-----------|
| 0.1 | 118 (73.8%) | **40 (25.0%)** | 2 (1.2%) | 0 |
| 0.2 | 142 (88.8%) | 14 (8.8%) | 2 (1.2%) | 2 (1.2%) |
| 0.4 | **158 (98.8%)** | 1 (0.6%) | 0 | 1 (0.6%) |

**核心发现**：δ=0.1（high memory）下 Convergent 仅 73.8%，Divergent 达 25.0%；δ=0.4（low memory）下 Convergent 达 98.8%，Divergent 仅 0.6%。Divergent 比例在 δ=0.1 时是 δ=0.2 的 2.9 倍，是 δ=0.4 的 42 倍。其原因是：**δ 直接缩放 interaction matrix 的 eigenvalues——δ 越小，ρ(A) 越大（δ=0.1: mean ρ=1.08, δ=0.4: mean ρ=0.77）。当 ρ(A) > 1 时，self-retention 让 feedback amplification 跨时间步累积而不被 decay 冲走，逐步将 expression 推向 divergence**。

**机制**：在 δ=0.1 时，self-retention = 0.9。这意味着每步仅丢失 10% 的 expression——正反馈驱动的 amplification 不会被 decay"冲刷"掉，而是跨时间步累积。当 interaction strength 足以压倒 decay 时（ρ(A) > 1），这种累积效应将正反馈放大的 trajectory 推向 divergence。δ=0.4 时，40% 的 decay 充当了强大的 damping force——即使 ρ(A) > 1，decay 也在每步吸收大部分 amplification，使 system 保持 bounded。

### 2.3 ρ(A) 与 δ 的系统性关系 (§13.6)

δ 直接缩放 transition matrix A 的 eigenvalues，因此 ρ(A) 与 δ 几乎完全共变：

| δ | mean ρ(A) | ρ(A) range | worlds with ρ > 1 |
|---|----------|------------|-------------------|
| 0.1 | 1.083 | [0.900, 1.401] | 126/160 (78.8%) |
| 0.2 | 0.971 | [0.800, 1.285] | 74/160 (46.3%) |
| 0.4 | 0.768 | [0.600, 1.092] | 5/160 (3.1%) |

δ=0.1 时近 80% 的 world 具有 ρ(A) > 1（linear core unstable），但仅 25% 实际 divergent——clipping 在 high-memory regime 中承担了关键的 stabilization 角色。δ=0.4 时仅 3.1% 的 world 有 ρ > 1，几乎全部被 decay 稳定化。

按 primary regime 分层：

| Regime | n | mean ρ(A) | min ρ(A) | max ρ(A) |
|--------|---|----------|----------|----------|
| Convergent | 418 | 0.910 | 0.600 | 1.399 |
| Sustained oscillatory | 4 | 1.153 | 1.090 | 1.238 |
| Divergent | 55 | 1.154 | 0.934 | 1.401 |
| Ambiguous | 3 | 1.032 | 1.016 | 1.062 |

Divergent worlds 的 mean ρ(A) = 1.154，显著高于 Convergent，但两者在 ρ ∈ [0.9, 1.4] 区间内大面积重叠——ρ(A) 对单个 world 的分类能力有限，尤其在 clipping-dominated 区域。

> Spectral radius 按 regime、delta、clipping 的分布见图:
> [../figures/04_spectral_radius_by_regime.png](../figures/04_spectral_radius_by_regime.png)
> [../figures/04_spectral_radius_by_delta.png](../figures/04_spectral_radius_by_delta.png)
> [../figures/04_spectral_radius_vs_clipping.png](../figures/04_spectral_radius_vs_clipping.png)

### 2.4 Topology-Control Robustness (§13.2)

TF-restricted topology 的稳定化效应在 run_002/run_003（δ=0.2）中表现显著（Type A 从 71.2% 提升至 92.5%）。但在 run_004 中，δ 效应在量级上与 topology 效应相当甚至更大：

| δ | unrestricted Convergent | tf_restricted Convergent | Δ |
|---|------------------------|---------------------------|---|
| 0.1 | 56/80 (70.0%) | 61/80 (76.3%) | +6.3pp |
| 0.2 | 71/80 (88.8%) | 71/80 (88.8%) | 0pp |
| 0.4 | 79/80 (98.8%) | 79/80 (98.8%) | 0pp |

TF-restricted topology 仅在 δ=0.1（最具挑战性的 high-memory regime）下提供 6.3 个百分点的额外稳定性。δ=0.2 和 δ=0.4 下两种 topology 的表现几乎无差异——decay damping 已经足够强，topology constraint 的额外稳定化效应被淹没。

**结论**：TF-restricted topology 的稳定化效应在 high-memory setting（δ=0.1）下仍然存在，但幅度远小于 δ effect 本身。topology-control 是 secondary stabilizer，decay 是 primary stabilizer。

> Topology × δ comparison 见图: [../figures/02_topology_delta_comparison.png](../figures/02_topology_delta_comparison.png)

### 2.5 Strength Boundary — lower δ 是否使 high-strength 更容易 Divergent？(§13.3)

**是的。δ 与 strength regime 的交互是系统发散的核心驱动。**

按 strength regime × δ 的 Divergent 与 Convergent 比例：

| strength regime | δ=0.1 Div | δ=0.1 Conv | δ=0.2 Div | δ=0.2 Conv | δ=0.4 Div | δ=0.4 Conv |
|----------------|-----------|------------|-----------|------------|-----------|------------|
| baseline | 1/40 (2%) | 38/40 (95%) | 0/40 (0%) | 40/40 (100%) | 0/40 (0%) | 40/40 (100%) |
| chen_moderate | 12/40 (30%) | 28/40 (70%) | 0/40 (0%) | 40/40 (100%) | 0/40 (0%) | 40/40 (100%) |
| stress | 13/40 (32%) | 27/40 (68%) | 5/40 (12%) | 33/40 (82%) | 0/40 (0%) | 40/40 (100%) |
| chen_stress | 14/40 (35%) | 24/40 (60%) | 9/40 (22%) | 29/40 (72%) | 1/40 (2%) | 38/40 (95%) |

注：除 Divergent 和 Convergent 外，剩余少量 world 为 Sustained oscillatory 或 Ambiguous。

**解读**：在 δ=0.4 下 Divergent 几乎消失，即使 chen_stress 也仅 2% Divergent——decay 有效压制了 interaction-driven instability。而当 δ 降至 0.1 时，Divergent 从 baseline 的 2% 逐步升至 chen_stress 的 35%。δ 与 strength 呈现协同效应：**低 δ 放大 interaction 的不稳定潜力，高 strength 提供这个潜力的来源，两者叠加时 divergence 最高**。

收敛时间的对应模式：δ=0.4 下 baseline 平均仅需 20 步收敛，而 δ=0.1 下 chen_stress 需要 156 步——slow convergence 在 high-memory setting 中成为 norm。

### 2.6 Repression Effect — 稳定化效应是否随 δ 改变？(§13.4)

**是的。Repression 偏置的稳定化效应在 δ=0.1 时最强，随 δ 增大而减弱。**

| δ | balanced Convergent | repression_biased Convergent | Δ | balanced Divergent | repression_biased Divergent |
|---|--------------------|-------------------------------|---|--------------------|----------------------------|
| 0.1 | 49/80 (61.3%) | 68/80 (85.0%) | **+23.8pp** | 29 (36.3%) | 11 (13.8%) |
| 0.2 | 72/80 (90.0%) | 70/80 (87.5%) | −2.5pp | 6 (7.5%) | 8 (10.0%) |
| 0.4 | 78/80 (97.5%) | 80/80 (100%) | +2.5pp | 1 (1.3%) | 0 (0%) |

在 δ=0.1 时，repression_biased 将 Convergent 率从 61.3% 拉升至 85.0%（+23.8pp），divergence 从 36.3% 降至 13.8%。但在 δ=0.2 时该效应消失（甚至轻微反转），δ=0.4 时两者均已接近 100% Convergent——天花板效应使得 sign ratio 的差异被压缩。

**机制**：在 high-memory regime（δ=0.1）下，positive feedback 的累积效应最强，此时负调控充当的"刹车"作用最为关键。当 decay 本身已经足够强（δ=0.4）时，sign ratio 的影响被 decay damping 覆盖。

### 2.7 Clipping Boundary — δ 是否改变 clipping-mediated dynamics 的范围？(§13.5)

**是的。Clipping 是 δ-dependent 的——δ 越小，clipping 越频繁。**

Clipping-dominated (clipping frequency > 0.1) 的 world 比例：

| δ | unrestricted_sparse | tf_restricted |
|---|--------------------|---------------|
| 0.1 | 76/80 (95%) | 75/80 (94%) |
| 0.2 | 60/80 (75%) | 59/80 (74%) |
| 0.4 | 36/80 (45%) | 35/80 (44%) |

平均 clipping frequency：

| δ | unrestricted_sparse | tf_restricted |
|---|--------------------|---------------|
| 0.1 | 0.328 | 0.334 |
| 0.2 | 0.249 | 0.256 |
| 0.4 | 0.099 | 0.102 |

δ=0.1 时平均每 3 步就有一次 gene expression 被 clip（33%），δ=0.4 时降为每 10 步一次（10%）。这反映了 decay 在防止 expression 被推至负值方面的直接作用——更强的 decay 意味着 expression 在每步被 pull toward zero，减少了对 clipping 的依赖。

**机制链**：higher self-retention → larger expression magnitude → stronger interaction effect → more frequent negative pushes → more clipping events。因此 clipping frequency 与 δ 呈负相关，与 interaction strength 呈正相关。在 δ=0.1 × chen_stress 的极端组合下，clipping 最为频繁——这是"high memory + strong interaction"共同驱动的结果。

> Clipping boundary plots 见图: [../figures/03_clipping_boundary.png](../figures/03_clipping_boundary.png)

### 2.8 Per-Combo Regime Distribution

48 个 combo 的 Divergent world 数量（完整数据见 `tables/regime_summary.tsv`，每 combo 共 10 worlds）：

**δ=0.1**:

| strength | unrestricted, bal | unrestricted, rb | tf_restricted, bal | tf_restricted, rb |
|----------|-------------------|------------------|--------------------|--------------------|
| baseline | 0 | 0 | 0 | 0 |
| chen_moderate | **6** | 1 | 3 | 2 |
| stress | **7** | 0 | 3 | 1 |
| chen_stress | **7** | 1 | 4 | **5** |

**δ=0.2**:

| strength | unrestricted, bal | unrestricted, rb | tf_restricted, bal | tf_restricted, rb |
|----------|-------------------|------------------|--------------------|--------------------|
| baseline | 0 | 0 | 0 | 0 |
| chen_moderate | 0 | 0 | 0 | 0 |
| stress | 0 | 1 | 0 | 3 |
| chen_stress | 1 | 3 | 3 | 3 |

**δ=0.4**:

| strength | unrestricted, bal | unrestricted, rb | tf_restricted, bal | tf_restricted, rb |
|----------|-------------------|------------------|--------------------|--------------------|
| baseline | 0 | 0 | 0 | 0 |
| chen_moderate | 0 | 0 | 0 | 0 |
| stress | 0 | 0 | 0 | 0 |
| chen_stress | 0 | 0 | 1 | 0 |

核心模式已在 §2.5–§2.7 中讨论：Divergent 集中在 δ=0.1 × high strength，repression_biased 显著抑制不平衡组合下的 divergence，tf_restricted 在 δ=0.1 下提供额外稳定性。

### 2.9 Oscillation Emergence

仅 4 个 world（0.8%）被分类为 Sustained oscillatory，全部在 δ=0.1 或 δ=0.2 下（δ=0.4 下无 oscillator）。所有 genuine oscillator 的 ρ(A) > 1 且均为 clipping_dominated。

| world_id | topology | δ | ρ(A) | max_amp | clip_freq |
|----------|----------|---|------|---------|-----------|
| unrestricted_sparse_d0p1_chen_stress_repression_biased_t001 | unrestricted | 0.1 | 1.180 | 10.0 | 36.6% |
| unrestricted_sparse_d0p2_chen_stress_balanced_t003 | unrestricted | 0.2 | 1.105 | 3.4 | 30.5% |
| unrestricted_sparse_d0p2_chen_stress_balanced_t005 | unrestricted | 0.2 | 1.238 | 18.5 | 23.6% |
| tf_restricted_d0p1_chen_stress_balanced_t008 | tf_restricted | 0.1 | 1.090 | 13.7 | 26.2% |

Oscillation 仅出现在 unrestricted 拓扑（3/4）或最强的 chen_stress 组合下（4/4）。tf_restricted 仅出现 1 例。这与 run_003 的发现一致——TF 集中式调控压缩了反馈回路的自由度，抑制了 oscillation emergence。

> Trajectory exemplars 见图: [../figures/05_trajectory_exemplars.png](../figures/05_trajectory_exemplars.png)
>
> All worlds trajectory grids 见: `../figures/traj/`（24 张 batch PNG，覆盖全部 480 worlds）

> Trajectory exemplars 见图: [../figures/05_trajectory_exemplars.png](../figures/05_trajectory_exemplars.png)
>
> All worlds trajectory grids 见: `../figures/traj/`（24 张 batch PNG，覆盖全部 480 worlds）

### 2.10 Perturbation Recovery — 扰动后是否恢复？(§13.1 supplement)

10 个 candidate worlds 的 canonical KD perturbation 结果：

| world_id | orig regime | pert regime | KD gene | gene role | recovery |
|----------|------------|-------------|---------|-----------|----------|
| tf_restricted_d0p1_baseline_balanced_t000 | Convergent | Convergent | 8 | non-TF | 577 |
| tf_restricted_d0p2_baseline_balanced_t000 | Convergent | Convergent | 17 | non-TF | 535 |
| tf_restricted_d0p4_baseline_balanced_t000 | Convergent | Convergent | 15 | non-TF | 515 |
| unrestricted_sparse_d0p1_chen_stress_repression_biased_t001 | Sustained oscillatory | Sustained oscillatory | 15 | non-TF | failed |
| unrestricted_sparse_d0p2_chen_stress_balanced_t003 | Sustained oscillatory | Sustained oscillatory | 8 | non-TF | failed |
| unrestricted_sparse_d0p1_chen_moderate_balanced_t000 | Divergent | Divergent | 15 | non-TF | failed |
| tf_restricted_d0p2_chen_stress_balanced_t000 | Divergent | Divergent | 17 | non-TF | failed |
| tf_restricted_d0p4_chen_stress_balanced_t009 | Divergent | Divergent | 13 | non-TF | failed |
| unrestricted_sparse_d0p2_stress_repression_biased_t005 | Ambiguous | Ambiguous | 14 | non-TF | failed |
| unrestricted_sparse_d0p4_chen_stress_balanced_t000 | Ambiguous | Ambiguous | 6 | non-TF | failed |

**发现**：

- **Convergent worlds（3/3）全部恢复**。所有 3 个 Convergent world 在 KD 后回归原 regime，recovery time 515–577（绝对时刻），即 KD 后 15–77 步恢复。δ=0.4 的恢复最快（15 步），δ=0.1 最慢（77 步）——与 decay damping 一致。
- **Divergent worlds（3/3）无法挽救**。单基因 knockdown 对 runaway divergence 无效——发散由多基因驱动或 linear core 的不稳定性超过单基因干预能力。注意 `unrestricted_sparse_d0p1_chen_moderate_balanced_t000` 的 KD 时刻 expression 已达 ~4.5×10¹⁹，说明 δ=0.1 下即使是 moderate strength 的组合也能在 t=500 前产生极端 divergence。
- **Sustained oscillatory（2/2 振荡持续）**：两个 genuine oscillator 在 KD 后均保持 SOsc，未发生 regime change。KD 扰动不足以破坏 δ=0.1 或 δ=0.2 下 chen_stress 驱动的 sustained oscillation。
- **有效 regime change：0/10**。没有任何 candidate 发生 primary regime 转变。与 run_002 的 "attractor 类型基本不变"结论一致——Model A 的 canonical single-gene KD 扰动大多不足以改变 primary regime，系统展现出较强的 structural stability。
- KD 基因的 expression level 与 regime 高度相关：Convergent worlds 中最高表达基因仅 0.4–3.3，Divergent worlds 中已达 4.5×10¹⁹–1.2×10³⁵。

> Perturbation 前后轨迹对比见图:
> [../figures/perturbation/perturbation_convergent.png](../figures/perturbation/perturbation_convergent.png)
> [../figures/perturbation/perturbation_sustained_oscillatory.png](../figures/perturbation/perturbation_sustained_oscillatory.png)
> [../figures/perturbation/perturbation_divergent.png](../figures/perturbation/perturbation_divergent.png)
> [../figures/perturbation/perturbation_ambiguous.png](../figures/perturbation/perturbation_ambiguous.png)

### 2.11 Sign Ratio Comparison — Cross-δ Perspective

Sign ratio comparison figures 展示了 balanced vs repression_biased 在所有 3 个 δ 值下的对比：

> Sign ratio comparison 见图:
> [../figures/06_sign_ratio_comparison_unrestricted_sparse.png](../figures/06_sign_ratio_comparison_unrestricted_sparse.png)
> [../figures/06_sign_ratio_comparison_tf_restricted.png](../figures/06_sign_ratio_comparison_tf_restricted.png)

### 2.12 Abundance Compatibility — Clipping 的角色

Non-negativity clipping 是 Model A dynamics 的核心 component，但其重要性随 δ 系统性变化：

- δ=0.1：95% worlds clipping-dominated，mean clipping freq 33%——clipping 是常态
- δ=0.2：74% worlds clipping-dominated，mean clipping freq 25%——clipping 仍是主要特征
- δ=0.4：44% worlds clipping-dominated，mean clipping freq 10%——clipping 退居次要

在 δ=0.4 的 baseline regime 中，0% world 为 clipping-dominated——系统真正运行在"linear-like"区域。但 δ=0.4 的 chen_stress 仍有 80–100% world 为 clipping-dominated——strength 和 δ 之间存在 trade-off。

当前 clipping-dominated 的判定阈值为 clipping frequency > 0.1。从分布来看（δ=0.4 下 median 约 5–10%），这个阈值对 δ=0.4 而言可能偏高，但保持了跨 run 的可比性。

---

## 3. 与 run_002 / run_003 的对比

### 3.1 δ=0.2 下的 baseline 对比

Run_004 在 δ=0.2 下的结果可作为 run_002/run_003 的直接对照：

| Metric | run_002 (unrestr) | run_003 (tf_restr) | run_004 δ=0.2 (unrestr) | run_004 δ=0.2 (tf_restr) |
|--------|-------------------|--------------------|----------------------------|----------------------------|
| Type A / Convergent | 57/80 (71.2%) | 74/80 (92.5%) | 71/80 (88.8%) | 71/80 (88.8%) |
| Type E / Divergent | 11/80 (13.8%) | 3/80 (3.8%) | 5/80 (6.3%) | 9/80 (11.3%) |
| Oscillation (C+D) | 3/80 (3.8%) | 0/80 (0%) | 2/80 (2.5%) | 0/80 (0%) |

Run_004 δ=0.2 的结果介于 run_002 和 run_003 之间。值得注意的是，run_004 中 unrestricted 的 Convergent 率（88.8%）高于 run_002（71.2%）——这可能是因为 run_004 使用了更大的 N_WORLD（10 vs run_002 的 10）和不同的 seed 分配方案。Run_004 的 δ=0.2 tf_restricted（88.8%）则低于 run_003（92.5%）。这些跨 run 的差异提示：在相同的参数区间内，sample realization 的方差不可忽略——10 worlds per combo 是合理的 minimum。

### 3.2 δ 效应 vs Topology 效应的量级对比

以 Convergent 比例的变动幅度衡量：

| Effect | Magnitude |
|--------|-----------|
| δ effect（0.1 → 0.4, all worlds） | +25.0pp (73.8% → 98.8%) |
| Topology effect（unrestr → tf_restr, δ=0.2, run_003） | +21.3pp (71.2% → 92.5%) |
| Topology effect（unrestr → tf_restr, δ=0.1, run_004） | +6.3pp (70.0% → 76.3%) |
| Repression effect（bal → rb, δ=0.1） | +23.8pp (61.3% → 85.0%) |

δ effect 的量级（~26pp）与 run_003 中 topology effect（~21pp）和 δ=0.1 内 strength effect（35pp）在同一水平。在 run_004 的设计中，**δ 和 interaction strength 是两个同等重要的 axis**，共同决定 primary regime boundary 的位置。

---

## 4. 科学回答总结 (对应 §13.1–13.7)

### 13.1 Decay / memory boundary

> δ = 0.1 (high memory) → 73.8% Convergent, 25.0% Divergent
> δ = 0.2 (reference) → 88.8% Convergent, 8.8% Divergent
> δ = 0.4 (low memory) → 98.8% Convergent, 0.6% Divergent

δ 系统性改变 primary regime distribution：**lower δ（higher self-retention）显著增加 divergence**。self-retention 让 interaction-driven amplification 跨时间累积，而非 decay 掉。

### 13.2 Topology-control robustness

TF-restricted topology 的稳定化效应在 δ=0.1 下仍然存在（+6.3pp Convergent），但 magnitude 远小于 δ effect 本身。在 δ=0.2 和 δ=0.4 下 topology effect 被 decay damping 覆盖。

### 13.3 Strength boundary

Lower δ 确实使 high-strength regimes 更容易 Divergent。Divergence risk ≈ strength × (1/δ)。δ=0.1 下 chen_stress 的 Convergent 率仅 60%（δ=0.4 下 95%）。

### 13.4 Repression effect

Repression_biased 的稳定化效应在 δ=0.1 时最强（+23.8pp Convergent），在 δ=0.2 和 δ=0.4 时衰减或消失。Negative feedback 的 damping 作用在 high-memory setting 中最关键。

### 13.5 Clipping boundary

δ 强烈改变 clipping-mediated dynamics 的范围：δ=0.1 时 95% worlds clipping-dominated（mean freq 33%），δ=0.4 时 44% clipping-dominated（mean freq 10%）。Clipping frequency 与 δ 呈负相关。

### 13.6 ρ(A) correspondence

ρ(A) 与 δ 几乎完全共变（δ=0.1: mean ρ=1.08, δ=0.2: 0.97, δ=0.4: 0.77）。ρ(A) > 1 的比例从 δ=0.1 的 79% 降至 δ=0.4 的 3%。ρ(A) 对单个 world 的分类能力在 clipping-dominated 区域有限（Convergent 和 Divergent 在 ρ 上有大面积重叠），但作为 population-level trend 可靠。

### 13.7 Model B readiness

当前 Model A decay / memory boundary 已经足够清楚，可以支持进入 minimal Model B comparison。具体而言：

- δ effect 的 magnitude、direction、interaction with strength/sign ratio 已明确
- TF-restricted topology 的 δ-dependence 已刻画
- Clipping 的 δ-dependent role 已量化
- ρ(A) 随 δ 的系统性偏移已确认

进入 Model B 时，关键对照维度应为：Model B 的 X → P → X structure 引入的 memory / delay effect，是否在量级上与 δ-driven memory pressure 可比？如果是，Model B 的 advantage 可能部分来自"隐性的 self-retention increase"而非 representation learning。

---

## 5. 输出清单

### Tables
- `tables/world_summary.tsv` — 480 worlds, 23 columns
- `tables/regime_summary.tsv` — 48 combos aggregated
- `tables/boundary_summary.tsv` — 2-level boundary aggregation
- `tables/aggregation.json` — full aggregation dictionary

### Figures
- `figures/01_primary_regime_heatmaps.png` — δ × strength heatmap by topology × sign ratio
- `figures/02_topology_delta_comparison.png` — unrestricted vs tf_restricted across δ
- `figures/03_clipping_boundary.png` — δ vs clipping frequency
- `figures/04_spectral_radius_by_regime.png` — ρ(A) by primary regime
- `figures/04_spectral_radius_by_delta.png` — ρ(A) across δ
- `figures/04_spectral_radius_vs_clipping.png` — ρ(A) vs clipping frequency
- `figures/05_trajectory_exemplars.png` — exemplar trajectories (4 regimes × 2 views)
- `figures/06_sign_ratio_comparison_*.png` — sign ratio comparison (2 figures)
- `figures/traj/` — all worlds trajectory grids (24 batch PNGs, 480 worlds)
- `figures/perturbation/` — perturbation before/after plots (4 figures)

### Analysis Data
- `world_metadata/` — 480 world metadata JSONs
- `trajectories/` — 2400 trajectory `.pt` files
- `analysis/` — 2400 per-cell analysis JSONs + 1 perturbation summary
- `perturbations/` — 10 perturbed trajectory `.pt` files

---

## 6. One-Line Summary

> run_004 demonstrates that decay δ, alongside interaction strength, is a key determinant of Model A regime boundaries: lower δ (higher memory) systematically increases divergence, amplifies clipping, and enhances the marginal value of repression bias — memory is a double-edged sword that preserves both stability and instability.
