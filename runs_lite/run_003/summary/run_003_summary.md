# Run 003 — Summary Report

**Stage**: DDC Lite Representation Curriculum — Model A Topology-Control Run

**Run**: run_003 (TF_Restricted_Topology_Control)

**Model**: Model A — Classical Linear GRN with TF-restricted topology

**Date**: 2026-06-08

---

## 0. Executive Summary

> 在保持 run_002 所有 parameter regime（4 strength × 2 sign ratio = 8 combos，每 combo 10 worlds，共 80 worlds）不变的前提下，将 topology 生成规则从"任意 gene 可作为 regulator"切换为"仅 5 个 TF gene 可作为 regulator"（TF-restricted topology），同时保持 edge count 与 run_002 一致（N_EDGES=38，r=0.1）。
>
> **核心结论**：TF-restricted topology 使系统显著更稳定——Type A（stable equilibrium）从 run_002 的 71.2% 提升至 92.5%，divergence（Type E）从 13.8% 降至 3.8%，oscillation（Type C/D）完全消失。system 处于高度 clipping-dominated 状态（71.3% worlds），non-negativity clipping 在 TF 集中式调控结构上起到非线性 stabilizer 的作用。repression 偏置依然将 divergence 转化为 slow convergence，而非完全消除不稳定性。

---

## 1. 科学问题与实验设计

### 1.1 核心科学问题

> *在尽量保持 run_002 runtime / parameter semantics 不变时，TF-restricted topology 是否改变 stability、oscillation、convergence 或 perturbation recovery？*

同时必须回答：

> *在 edge count / effective edge density 与 run_002 匹配时，regulator eligibility restriction 是否改变 clipping-mediated dynamics？*

### 1.2 实验定位

Run 003 是 run_002 的 **topology-control 对照实验**。run_002 建立了 unconstrained Chen-style sparse GRN 的 dynamical baseline；run_003 在完全相同的 parameter grid 上，唯一改变 topology 生成规则——所有 regulatory edge 必须从固定的 5 个 TF gene 出发。这使得我们可以**隔离 topology constraint 的效应**，回答以下问题：将调控权限制在少数 TF gene 上，是增强还是削弱系统稳定性？

### 1.3 参数网格（与 run_002 完全一致）

沿两个维度扫描：**interaction strength**（4 档）和 **sign ratio**（2 档），共 8 个 parameter combos，每 combo 10 个独立 topology，共 80 worlds × 5 cells = 400 trajectories。

固定参数：δ=0.2, b=0.1, G=20, N_EDGES=38, r=0.10。

| strength regime | a ∈ [min, max] | 设计意图 |
|---|---|---|
| `baseline` | [0.02, 0.10] | 保守稳定基线，strength 远小于 decay（δ=0.2），必然 decay-dominated |
| `chen_moderate` | [0.10, 0.20] | 对应 Chen (2019) 的模拟强度区间 |
| `stress` | [0.10, 0.30] | 压力测试，strength 上界已超过 decay |
| `chen_stress` | [0.20, 0.40] | 极端条件，strength 可达 decay 两倍 |

| sign ratio | act:rep | 含义 |
|---|---|---|
| `balanced` | 1:1 | activation 与 repression 各半 |
| `repression_biased` | 1:2 | repression 是 activation 的两倍 |

### 1.4 TF-restricted Topology 的核心变化

| 属性 | Run 002 (unconstrained) | Run 003 (TF-restricted) |
|---|---|---|
| regulator 来源 | 任意 20 个 gene | 仅 5 个 TF gene |
| candidate edge pool | C(20,2) = 380 | 5 × 20 − 5 = 95 |
| edge 总数 | 38 | 38（与 run_002 匹配） |
| TF core edge 数 | 不适用 | ~9（TF-TF 之间） |
| TF → non-TF edge 数 | 不适用 | ~29 |
| TF pool internal density | 不适用 | ~0.4 |
| 平均 in-degree | ~1.9 | ~1.9 |

关键约束：**non-TF gene 不能作为任何 gene 的 regulator**。这意味着网络的因果结构被压缩为"TF hub → 全体 gene"的双层架构，feedback loop 只能通过 TF-TF 互调控形成。

> TF-restricted 拓扑的代表性结构见图: [../figures/topology_exemplars.png](../figures/topology_exemplars.png)

### 1.5 与 run_002 的 world 配对规则

`world_summary.tsv` 中包含两个字段用于关联 run_002：

- **`run_002_comparison_label`**：对所有 world 取值为 `edge-density-matched (N_EDGES=38)`，表明两个 run 在 edge 总数层面匹配（均采样 38 条 edge），对照的是"相同连接密度下，regulator eligibility 是否改变 dynamics"。
- **`paired_run002_world_id`**：按 `(combo_key, topo_idx)` 确定 run_002 中对应的 world。即第 i 个 parameter combo 下的第 j 个 topology 在 run_002 中的 world_id。需要注意，这里的配对是 **grid 位置配对**而非**拓扑配对**：run_002 和 run_003 采用不同的 topology 生成规则（unconstrained vs TF-restricted），即使 `topo_idx` 相同，两者生成的邻接矩阵也完全不同。因此 §2.9 中的 pair-wise 对比本质上比较的是"相同 parameter grid 位置上，两种 topology 生成策略下的 attractor 差异"，而非同一拓扑在两个 run 下的前后变化。

---

## 2. 结果

### 2.1 全局 Attractor 分布

80 worlds 的 primary attractor 分布：

| Type | 名称 | 数量 | 占比 | Run 002 对比 |
|---|---|---|---|---|
| A | Stable equilibrium | 74 | 92.5% | 71.2% ↑ |
| B | Slow convergence | 3 | 3.8% | 11.2% ↓ |
| C | Damped oscillation | 0 | 0% | 2.5% ↓ |
| D | Sustained oscillation | 0 | 0% | 1.2% ↓ |
| E | Runaway divergence | 3 | 3.8% | 13.8% ↓ |
| F | Numerical collapse | 0 | 0% | 0% — |

### 2.2 按 Combo 的 Attractor 分布

| strength | sign ratio | Type A | Type B | Type E |
|---|---|---|---|---|
| baseline | balanced | 100% | — | — |
| baseline | repression_biased | 100% | — | — |
| chen_moderate | balanced | 90% | 10% | — |
| chen_moderate | repression_biased | 100% | — | — |
| stress | balanced | 90% | — | 10% |
| stress | repression_biased | 90% | 10% | — |
| chen_stress | balanced | 80% | — | 20% |
| chen_stress | repression_biased | 90% | 10% | — |

> 按 strength regime 的 attractor 分布见图: [../figures/strength_vs_regime_type.png](../figures/strength_vs_regime_type.png)
>
> 各 attractor type 的 exemplar trajectory 见图: [../figures/canonical_trajectories.png](../figures/canonical_trajectories.png)
>
> convergence time 分布见图: [../figures/convergence_time.png](../figures/convergence_time.png)

### 2.3 Stability Geometry — 哪些 parameter geometry 更容易稳定？（§8.1）

**interaction strength 是稳定性的主控轴。**

- **baseline**（a ∈ [0.02, 0.1]）：20/20 worlds 全为 Type A。strength 远小于 decay，系统被 decay-dominated，无论 sign ratio 如何都必然收敛。
- **chen_moderate**（a ∈ [0.1, 0.2]）：95% Type A。仅 1 个 balanced world 出现 Type B（slow convergence）。
- **stress**（a ∈ [0.1, 0.3]）：稳定性开始退降。balanced 下出现 10% Type E divergence。
- **chen_stress**（a ∈ [0.2, 0.4]）：balanced 下 Type E 升至 20%。

与 run_002 的关键差异：run_002 在 stress 和 chen_stress 下已有大量 Type E（分别 20% 和 35%），run_003 的 divergence 率仅为 run_002 的 1/2 至 1/4。TF 集中式调控结构显著提高了同等 strength 下的稳定性。

这一发现与 run_002 §3.9 的拓扑分析形成呼应。run_002 发现，in-degree 分布越不均匀（hub-likeness 越强）的 topology，Type A 比例越高、Type E 比例越低：在 hub-likeness 最高的 Q4 区间，Type A 达到 18/20 且 Type E 完全消失，而 hub-likeness 较低或中等的区间（Q1–Q3）反而出现更多 divergence 和 oscillation。当时给出的解释是——hub gene 作为 dominant sink/source 迫使系统偏向 single-equilibrium 结构。run_003 的 TF-restricted topology 正是将这种"hub-centric"结构推向了设计上的极限：不再是自然界涌现 hub，而是硬性规定所有 38 条 edge 必须从仅有的 5 个 TF gene 出发，hub-likeness 被结构性最大化。结果与 run_002 Q4 的趋势完全延续——Type A 升至 92.5%，Type E 降至 3.8%，oscillation 完全消失。这进一步支持了 run_002 的推断：**调控源越集中，dynamics 越偏向稳定均衡**。

**机制**：随着 interaction strength 增大，transition matrix A 的 spectral radius ρ(A) 接近并超过 1，linear core 出现 unbounded growth 的可能。但 TF-restricted topology 将 feedback 集中在 TF-TF 子图（仅 5 个 gene 参与互调控），限制了正反馈放大的扩散范围。同时，non-negativity clipping 在 TF 集中式调控结构上起到 stabilizer 作用——一旦某个 TF gene 被 clip 到 0，其下游所有 non-TF target 的驱动信号立即中断。

> in-degree 分布与 stability 的关系见图: [../figures/degree_distribution_vs_stability.png](../figures/degree_distribution_vs_stability.png)
>
> ρ(A) 按 attractor type 的分布见图: [../figures/spectral_radius.png](../figures/spectral_radius.png)

### 2.4 Oscillation Emergence — 是否出现 oscillation？（§8.2）

**80 worlds 中未观察到任何 oscillation（Type C/D = 0%）。**

run_002 中出现了 3.8% 的 oscillation（Type C 衰减振荡 2 例 + Type D 持续振荡 1 例），但在 run_003 中完全消失。

**机制**：oscillation 需要延迟负反馈环（delayed negative feedback loop）。在 run_002 的 unconstrained topology 中，任意 gene 之间可以形成长程反馈回路，这为振荡提供了结构基础。而在 run_003 中，所有 regulatory edge 必须从 TF gene 出发，网络简化为"TF → TF ∪ non-TF"的双层结构。负反馈只能通过 TF-TF 互调控形成，而 TF-TF 子图只有 5 个 node——如此小的子图无法支撑足以产生持续振荡的多步延迟反馈。

### 2.5 Repression Effect — repression-biased 系统是否更稳定？（§8.3）

**是的。repression 偏置将 divergence 转化为 slow convergence，而非简单还原为 Type A。**

| regime | balanced | repression_biased | 效应 |
|---|---|---|---|
| stress | 10% Type E | 0% Type E → 10% Type B | divergence 消失，转为 slow convergence |
| chen_stress | 20% Type E | 0% Type E → 10% Type B | divergence 消失，转为 slow convergence |

在出现 divergence 的两个 regime 中，将 sign ratio 从 1:1 切换到 1:2 后，divergence 完全消除。但系统并非直接回归 Type A——一部分 world 进入 Type B（slow convergence）。

**机制**：增加 repression 比例意味着 interaction matrix A 中有更多负值项。负调控在 dynamics 中起到"刹车"作用：正反馈驱动的 amplification 被负反馈部分抵消，使 spectral radius 降低。然而，当负调控过强时，gene expression 频繁被推至负值区间、触发 non-negativity clipping——clipping 虽然防止了数值发散，但也破坏了 convergence 的 smoothness，使收敛变慢。

这与 run_002 的发现定性一致，但 run_003 中 repression 的"稳定化但不完全恢复"效应更加干净：run_002 的 repression_biased 仍有 7.5% 的 Type E，而 run_003 为 0%。

> repression ratio 与 stability 的关系见图: [../figures/repression_ratio_vs_stability.png](../figures/repression_ratio_vs_stability.png)

### 2.6 Spectral Radius Correspondence — ρ(A) 与 runtime dynamics 是否一致？（§8.4）

**结论取决于 world 是否处于 clipping-dominated regime。**

先看整体分布：

| Attractor Type | Mean ρ(A) | ρ 范围 | n |
|---|---|---|---|
| Type A | 0.938 | [0.800, 1.264] | 74 |
| Type B | 0.991 | [0.981, 1.005] | 3 |
| Type E | 1.105 | [1.042, 1.164] | 3 |

进一步按是否 clipping-dominated 拆分后，模式更加清晰：

**非 clipping-dominated worlds（23 个，全部为 Type A）**：

23/23 的 ρ(A) < 1.0（最高仅 0.955）。ρ(A) 在这里是**完美的 stability predictor**：在 linear-like regime 下，只要 linear core 稳定（ρ < 1），runtime dynamics 就必然收敛。这是因为 dynamics 主要遵循 pure linear update，non-negativity clipping 几乎不参与。

**Clipping-dominated worlds（57 个）**：

| ρ(A) | Type A | Type B | Type E | 合计 |
|---|---|---|---|---|
| ρ < 1 | 34 | 2 | 0 | 36 |
| ρ > 1 | 14 | 1 | 3 | 18 |

- ρ < 1 时（36 个）：**仍然全部收敛**（34 Type A + 2 Type B），ρ(A) 继续保持预测力——即使 clipping 频繁发生，linear core 稳定时 system 不会发散。
- ρ > 1 时（18 个）：**predictability 破裂**。多数 world（14 Type A）仍由 clipping 稳定下来，但 3 个 Type E 的 world 表明 linear core 的不稳定性最终压倒了 clipping 的 stabilizer 能力。是否发散取决于 clipping 强度与 linear core 不稳定性之间的竞争。

**总结**：ρ(A) 的预测力与 clipping 角色直接相关：
- 非 clipping-dominated：ρ < 1 ⇔ 稳定，ρ > 1 未出现（完美 predictor）
- Clipping-dominated 且 ρ < 1：线性层面稳定，system 必然收敛（仍可靠）
- Clipping-dominated 且 ρ > 1：不确定——clipping 可能兜底（Type A），也可能兜不住（Type E）。ρ > 1 是 divergence 的必要条件而非充分条件。

> spectral radius 按 attractor type 的分布见图: [../figures/spectral_radius.png](../figures/spectral_radius.png)

### 2.7 Abundance Compatibility — non-negativity constraint 是否显著影响 dynamics？（§8.6）

**是的。系统是高度 clipping-dominated 的。**

- **57/80 worlds（71.3%）**被分类为 clipping-dominated。
- Clipping frequency 与 interaction strength 单调正相关：

| strength | balanced | repression_biased |
|---|---|---|
| baseline | 3.0% | 2.5% |
| chen_moderate | 19.4% | 21.4% |
| stress | 27.1% | 34.9% |
| chen_stress | 33.0% | 45.8% |

- 仅 baseline regime 保持 linear-like（clipping freq < 5%）。chen_moderate 及以上，clipping 成为 dynamics 的常态。
- 最高 clipping 出现在 chen_stress repression_biased（45.8%）：接近一半的 time step × gene 组合触发了 non-negativity projection。
- repression_biased 始终比 balanced 产生更高的 clipping frequency——强负调控频繁将 expression 推至负值。

**TF 集中式调控结构放大了 clipping 的效应**：因为所有 non-TF gene 的 regulator 都是 TF，一旦某个 TF 被 clip，其下游所有 target 立即失去驱动信号。这种"集中式 clipping"在 run_002 的分布式调控中不存在，是 TF-restricted topology 独有的稳定化机制。

**关于 clipping-dominated 判定标准的讨论**：当前以"任何一个 cell 的 clipping 次数超过 T×G×10%"作为 world 的 clipping-dominated 判据。10% 的阈值是否过于宽松？从 clipping frequency 的分布来看——中位数 25%、均值 23.4%、仅 20 个 world 低于 5%——阈值提升至 15% 或 20% 依然有超过半数 world 被分类为 clipping-dominated（分别 49/80 和 44/80）。退一步说，即使将 10–20% clipping freq 的 world 称为"clipping-active"而非"clipping-dominated"，也不改变核心结论：non-negativity projection 是绝大多数非 baseline world 中不可忽略的动力学成分。值得指出的是，baseline 以外的 regime 中 clipping freq 中位数为 ~27%，意味着每 4 个 time step 就有一次 gene expression 被 clip——这已经不是 edge case，而是 dynamics 的常规特征。

> clipping frequency 分布见图: [../figures/clipping_frequency.png](../figures/clipping_frequency.png)
>
> 各 regime 的 clipping-dominated 比例见图: [../figures/clipping_dominated_by_regime.png](../figures/clipping_dominated_by_regime.png)

### 2.8 Perturbation Recovery — 哪些 regimes 更容易恢复？（§8.5）

本 run 在 Type A、Type B、Type E 三类 attractor 中各选取 1 个 exemplar world，执行了 canonical perturbation 实验（t=500 时 knockdown 最高表达基因）：

| world | 原始类型 | knockdown gene | gene 类型 | 扰动后类型 | 恢复时间 |
|---|---|---|---|---|---|
| baseline_balanced_t000 | Type A | gene 14 | non-TF | **Type A（不变）** | 37 steps |
| chen_moderate_balanced_t003 | Type B | gene 17 | non-TF | **Type A（恢复为稳定）** | 51 steps |
| stress_balanced_t003 | Type E | gene 17 | non-TF | **Type E（不变，继续发散）** | 未恢复 |

**Type A world**：knockdown 后迅速恢复，attractor type 不变。system 在 baseline strength 下具有足够的 stability margin 来吸收单基因扰动。

**Type B world**：扰动后轨迹从 knockdown 后的状态起步，经历了与原始轨迹不同的瞬态过程，最终在 51 步内收敛，因此被分类为 Type A。这与 run_002 中 Type B→A 的现象一致：Type A/B 的区分仅在于收敛速度（≤200 vs >200 步），扰动后段以 knockdown 后状态为起点独立分类，其瞬态过程与原始轨迹的瞬态过程是两个不同的事情——前者恰好收敛更快，并不代表原始 Type B 的慢收敛被"治愈"，只说明 knockdown 后的新起点恰好较快地到达了 equilibrium。

**Type E world**：发散状态无法通过单基因 knockdown 挽救——linear core 的不稳定性（ρ=1.11）超越了 knockdown 的干预能力。

**局限**：knockdown target 的选取规则是 t=500 时 expression 最高的基因，因此 3 次 perturbation 恰好都命中了 non-TF gene，缺少对 TF gene knockdown 的观测。TF-restricted 拓扑下，knockdown TF gene 应产生远更强的 cascade 效应（因为所有 non-TF 的 regulator 都依赖该 TF），这是下一轮需要优先补充的实验。

> perturbation 前后轨迹对比见图:
>
> Type A: [../figures/perturbation/perturbation_type_a.png](../figures/perturbation/perturbation_type_a.png)
>
> Type B: [../figures/perturbation/perturbation_type_b.png](../figures/perturbation/perturbation_type_b.png)
>
> Type E: [../figures/perturbation/perturbation_type_e.png](../figures/perturbation/perturbation_type_e.png)

### 2.9 Topology Extension Effect — 与 run_002 的全面对比（§8.7）

**TF-restricted topology 系统性地增加了稳定性，减少了动力学多样性。**

| 指标 | Run 002 | Run 003 | 变化 |
|---|---|---|---|
| Type A (stable) | 71.2% | **92.5%** | ↑ 21.3pp |
| Type B (slow convergence) | 11.2% | 3.8% | ↓ 7.5pp |
| Type C/D (oscillation) | 3.8% | 0% | ↓ 3.8pp |
| Type E (divergence) | 13.8% | 3.8% | ↓ 10.0pp |
| Convergence rate | — | 95.0% | — |
| Bounded rate | — | 96.3% | — |

**按 strength regime 的对比**：

| strength | run_002 Type A | run_003 Type A | run_002 Type E | run_003 Type E |
|---|---|---|---|---|
| baseline | 100% | 100% | 0% | 0% |
| chen_moderate | 80% | 95% | 0% | 0% |
| stress | 50% | 90% | 20% | 5% |
| chen_stress | 55% | 85% | 35% | 10% |

在每个 regime 上，run_003 的 Type A 比例均高于 run_002，Type E 比例均低于 run_002。差距在 stress 和 chen_stress 上最为显著——这正是 run_002 中大量出现 divergence 和 oscillation 的 regime。

更精细的 pair-wise 对照来自按 `(combo_key, topo_idx)` 匹配的"同位置"比较。需要注意：run_002 和 run_003 采用不同的 topology 生成规则（unconstrained vs TF-restricted），因此即使 `world_id` 或 `(combo_key, topo_idx)` 相同，两组实际也对应着**不同的拓扑结构**。这里的比较不是在"同一拓扑"上的前后变化，而是在相同 parameter grid 位置上，两种 topology 生成策略所产生的 attractor 差异：

> run_002 → run_003 attractor 转移矩阵见图: [../figures/regime_shift_map.png](../figures/regime_shift_map.png)

| run_002 → run_003 | Type A | Type B | Type E | 合计 |
|---|---|---|---|---|
| Type A (stable) | **52** | 3 | 2 | 57 |
| Type B (slow) | **9** | 0 | 0 | 9 |
| Type C (damped osc) | **2** | 0 | 0 | 2 |
| Type D (sustained osc) | **1** | 0 | 0 | 1 |
| Type E (divergence) | **10** | 0 | 1 | 11 |

这个矩阵揭示了两种 topology 生成策略下的 attractor 分布差异：

- run_002 中 **全部 Type B/C/D 的 grid 位置（12 个），在 run_003 的 TF-restricted 生成策略下均变为 Type A**。
- run_002 中 Type E 的 grid 位置，**10/11（91%）在 run_003 中变为 Type A**。唯一仍为 Type E 的位置对应最强发散的参数条件（ρ=1.11）。
- 反向而言，2 个 run_002 Type A 的 grid 位置在 run_003 中变为 Type E——提示 TF 集中式调控结构并非在所有参数区间都有利，少数位置在失去非 TF 调控源的分散性后反而更不稳定。

总体 66.2%（53/80）的 grid 位置保持了相同的 attractor type，但变化方向**严重不对称**：27 个位置在 TF-restricted 下趋向更稳定（非 A → A），仅 5 个趋向更不稳定（A → 非 A）。这是对 TF-restricted topology 系统性稳定化效应最直接的 pair-wise evidence。

> run_002 vs run_003 多维度对比:
>
> attractor 分布对比见图: [../figures/run002_vs_run003_attractor.png](../figures/run002_vs_run003_attractor.png)
>
> ρ(A) 对比见图: [../figures/run002_vs_run003_spectral_radius.png](../figures/run002_vs_run003_spectral_radius.png)
>
> clipping 对比见图: [../figures/run002_vs_run003_clipping.png](../figures/run002_vs_run003_clipping.png)
>
> convergence 对比见图: [../figures/run002_vs_run003_convergence.png](../figures/run002_vs_run003_convergence.png)
>
> 各 topology × strength regime 的 attractor 热力图见图: [../figures/topo_strength_interaction.png](../figures/topo_strength_interaction.png)

**机制总结**：将 regulator pool 从 20 个 gene 压缩到 5 个 TF gene 产生了三重效应：
1. **减少 feedback loop 复杂性**：可能的 cycle 数量大幅减少，oscillatory motif 无法形成
2. **增大 effective diagonal dominance**：每条 edge 的 target 数量不变，但 source 高度集中，使 interaction matrix 更接近对角占优
3. **集中式 clipping stabilization**：TF 集中式调控结构使 clipping 的稳定化效应集中在少数关键 node 上

---

## 3. 讨论与局限

### 3.1 Perturbation 实验的覆盖范围

本 run 已执行 perturbation 实验（见 §2.8），但仅覆盖了 3 类 attractor type 各 1 个 exemplar world，全部 knockdown target 为 non-TF gene。最关键的缺失是**TF gene knockdown**——TF-restricted 拓扑的核心特征，下一轮应在更多 regime 下测试 TF vs non-TF knockdown 的差异化恢复行为。

### 3.2 Oscillation 的缺失

本 run 未观察到 oscillation，但这不意味着 TF-restricted topology 一定不能产生 oscillation。当前 parameter range 下所有 world 的 decay（δ=0.2）和 basal transcription（b=0.1）均为定值，oscillation 的出现可能与 gene 自身动力学参数（δ, b）和调控强度 a 之间的相对关系有关——例如当某些 gene 的 decay 较小、basal transcription 较低，同时又受到较强的 signed regulation 时，是否会出现更丰富的动态行为，是值得继续探究的方向。

### 3.3 核心发现总结

1. **TF-restricted topology 系统性地提升稳定性**：Type A 从 run_002 的 71.2% 升至 92.5%，divergence（Type E）从 13.8% 降至 3.8%，oscillation 完全消失。将调控源限制在少数 TF 上的结构效应与 run_002 §3.9 的自然 hub-likeness 趋势一脉相承——调控越集中，dynamics 越偏向稳定均衡。

2. **Interaction strength 是 regime 切换的主控轴，repression 将 divergence 转化为 slow convergence**：baseline（a ≤ 0.1）100% Type A，chen_stress（a ≤ 0.4）出现 20% Type E。repression_biased 在所有高强度 regime 中将 Type E 降为零，但代价是部分 world 进入 Type B。

3. **系统高度 clipping-dominated，clipping 在 TF 集中式调控结构上起非线性 stabilizer 作用**：71.3% worlds 处于 clipping-dominated，clipping frequency 随 interaction strength 单调上升。TF 集中式结构放大了 clipping 的效应——单个 TF 被 clip 即可中断其所有下游 target 的驱动信号。

4. **ρ(A) 的预测力取决于 clipping regime**：非 clipping-dominated world 中 ρ(A) 是完美 predictor（ρ < 1 ⇔ 稳定）；clipping-dominated 且 ρ > 1 时 predictivity 破裂——clipping 可能兜底（14 例 Type A）也可能兜不住（3 例 Type E）。

---

## 4. 输出清单

### 4.1 表格

| 文件 | 内容 |
|---|---|
| `tables/world_summary.tsv` | 80 worlds × 31 columns，覆盖 §7.1 全部字段 |
| `tables/regime_summary.tsv` | 8 combos × 15 columns，strength × sign → attractor 分布 |

### 4.2 图像

| 文件 | 内容 |
|---|---|
| `figures/canonical_trajectories.png` | 各 attractor type 的 exemplar trajectory |
| `figures/spectral_radius.png` | ρ(A) 按 attractor type 的分布 |
| `figures/strength_vs_regime_type.png` | interaction strength vs attractor regime |
| `figures/repression_ratio_vs_stability.png` | repression ratio vs stability |
| `figures/clipping_frequency.png` | clipping frequency 分布 |
| `figures/convergence_time.png` | convergence time 分布 |
| `figures/degree_distribution_vs_stability.png` | in-degree vs attractor type |
| `figures/topology_exemplars.png` | 代表性 TF-restricted 拓扑结构 |
| `figures/clipping_dominated_by_regime.png` | clipping-dominated 比例按 regime |
| `figures/run002_vs_run003_attractor.png` | run_002 vs run_003 attractor 对比 |
| `figures/run002_vs_run003_spectral_radius.png` | run_002 vs run_003 ρ(A) 对比 |
| `figures/run002_vs_run003_clipping.png` | run_002 vs run_003 clipping 对比 |
| `figures/run002_vs_run003_convergence.png` | run_002 vs run_003 convergence 对比 |
| `figures/regime_shift_map.png` | run_002 attractor → run_003 attractor 转移矩阵（per-world 混淆矩阵） |
| `figures/topo_strength_interaction.png` | 各 topology × 各 strength regime 的 attractor 热力图 |
| `figures/perturbation/perturbation_type_a.png` | Type A world 的 perturbation 前后轨迹 |
| `figures/perturbation/perturbation_type_b.png` | Type B world 的 perturbation 前后轨迹 |
| `figures/perturbation/perturbation_type_e.png` | Type E world 的 perturbation 前后轨迹 |
