# Run 002 — Summary Report

**Stage**: DDC Lite Representation Curriculum — Model A Foundation
**Run**: run_002 (Chen_Style_Sparse_GRN_Foundation)
**Runtime Version**: v0.2
**Date**: 2026-05

---

## 0. Executive Summary

> 在 Chen-style sparse topology（G=20, r=0.1）下，对 Classical Linear GRN 进行了系统性的 interaction strength × sign ratio 力学扫描：4 个强度区间（从保守基线到 Chen-inspired 压力水平） × 2 种正负调控比例 = 8 个 parameter combos，共 80 个 world。
>
> **核心结论**：interaction strength 是 attractor regime 切换的主控轴；repression-biased 系统在同等 strength 下显著更稳定；即使是最简单的 linear GRN，在强 interaction + 平衡正负调控条件下，也能涌现出衰减振荡、持续振荡、慢收敛和发散等多种动力学行为。

---

## 1. 科学问题与实验设计

### 1.1 核心科学问题

> *What kinds of dynamical regimes emerge in Classical Linear GRNs under controlled parameter geometry?*

### 1.2 实验定位

Run 002 是 **Model A Foundation Run**——不是 toy baseline，而是整个 DDC representation curriculum 的参考坐标系。目标是先建立一个 clean and interpretable 的 Classical Linear GRN dynamical universe，后续再研究 representation changes（Model B/C）如何偏离该 universe。

### 1.3 核心假设（来自 00_Problem_Definition §11）

| # | 假设 | 实验证据 | 结论 |
|---|------|----------|------|
| H1 | δ 固定为 0.2 的 decay-dominated regime 下，低 interaction 更稳定 | a_max ≤ 0.1 时，无论 sign ratio 如何，40/40 worlds 全为 Type A | **成立** |
| H2 | repression 比例提高会增强系统稳定性 | 全局 Type A 率：repression-biased (87.5%) > balanced (55%)；发散率：rb (7.5%) < balanced (20%) | **成立** |
| H3 | 纯 linear GRN 也可出现慢收敛和振荡行为 | 数据中出现了 Type B（慢收敛）9 例、Type C（衰减振荡）2 例、Type D（持续振荡）1 例 | **成立** |
| H4 | attractor type 主要由 parameter geometry 决定，而非仅由 topology 决定 | 60% (12/20) 的 topology 在扫过 4 个 strength regime 时切换了 attractor type | **成立** |

---

## 2. 方法

### 2.1 Model A 定义

单层 abundance-compatible piecewise-linear GRN：

$$
x_i(t+1) = \max\bigl(0,\; (1-\delta_i)x_i(t) + b_i + \sum_{j} a_{ij}x_j(t)\bigr)
$$

其中 `δ = 0.2`（self-retention = 0.8），`b = 0.1`。`a_ij = sign_ij × strength_ij`。non-negativity clipping（ReLU）使 Model A 在概念上成为 piecewise-linear system。

### 2.2 World 生成

- **Gene count**: G = 20
- **Topology**: Chen-style sparse directed graph, r = 0.1（38 条边 / 380 candidate，均匀抽样）
- **Edge signs**: 按 sign_ratio 分配 activation (+1) / repression (-1)
- **Edge strengths**: `a_ij ~ Uniform(a_min, a_max)`，sign 与 magnitude 分开采样
- **Basal / decay**: 全局固定 δ=0.2, b=0.1
- **Initial states**: 每个 world 5 个随机 cell（N_init=5），共 400 trajectories
- **Runtime length**: T = 1000 steps

### 2.3 Parameter Sweep

参数扫描沿两个维度展开：**interaction strength**（每条边的调控强度 |a_ij| 的采样区间，共 4 档）和 **sign ratio**（activation 与 repression 边的比例，共 2 档），交叉产生 8 个 parameter combos。每个 combo 下独立生成 10 个拓扑，保持 G=20、r=0.1、δ=0.2、b=0.1 不变。

#### 2.3.1 为什么要扫 strength？

核心问题：*interaction strength 从"远小于 decay"到"可与 decay 抗衡"时，系统动力学如何变化？*

δ=0.2（self-retention=0.8），对单基因而言，无调控时的天然衰减率为 0.2。当所有 incoming interaction 的总和远小于 0.2 时，系统被 decay-dominated——必然稳定。当 interaction 超过 decay 后，系统进入"调控与衰减竞争"的区域，可能出现 divergence、oscillation 和 slow convergence。

因此强度区间设计为一个 **从保守到极端的梯度**：

| regime | a 范围 | 设计意图 |
|--------|--------|----------|
| `baseline` | [0.02, 0.10] | 保守稳定基线。所有 edge 的强度总期望低于单基因的 decay（0.2），系统必然 decay-dominated。用于验证"在最温和条件下是否全为 Type A"。 |
| `stress` | [0.10, 0.30] | 压力测试。strength 上界（0.3）已超过 decay。用于测试：interaction 进入可与 decay 抗衡的范围后，系统何时失去稳定。 |
| `chen_moderate` | [0.10, 0.20] | 与 Chen (2019) 的模拟强度区间（a∈[0.1275, 0.1675]）对应。用于回答：Chen 的连续时间结论是否与我们的离散时间结果定性可比？ |
| `chen_stress` | [0.20, 0.40] | 超越 decay 两倍的极端条件。用于探测：在最强的 Chen-inspired 参数下，系统会完全崩溃还是仍有 survivor？ |

四个区间不是孤立的——它们形成了一条从 `baseline → chen_moderate → stress → chen_stress` 的强度递增链。每个链环节比前一环的强度期望高出约 50%-100%，使得我们可以观察到 **从全稳定到混合再到以发散为主的逐步退化**。

#### 2.3.2 为什么要扫 sign ratio？

核心问题：*增加 repression 比例是否能系统性地提升稳定性？*

两个条件来自 Chen (2019) 的实际数据：

- **balanced (1:1)**：activation 与 repression 各 50%。与 Chen 报道的 yeast (0.504:0.496) 和 human (0.53:0.47) 的实测比例一致。
- **repression_biased (1:2)**：repression 是 activation 的两倍。用于测试：当系统偏向负调控时，是否更强的"动态阻尼"效应？

通过固定 strength 区间、对比两个 sign ratio 下的 attractor 分布，可以直接量化 repression 偏置对稳定性的贡献。

#### 2.3.3 为什么每个 combo 扫 10 个 topology？

每个 combo 不是 1 个 world，而是 **10 个独立随机的 topology realization**。这 10 个 world 共享完全相同的参数（a_min, a_max, sign_ratio, δ, b），仅在**拓扑连接结构**和**每条边的具体 strength 值**上不同。

这意味着 combo 内的 variation 纯粹来自 topology 采样噪声，combo 间的差异才是 parameter 效应。10 个样本足够计算 attractor 分布比例、mean ρ(A)、convergence rate 等统计指标，同时保持计算量可控。

#### 2.3.4 实验设计的整体逻辑

```
                  interaction strength →
                  baseline  chen_moderate  stress  chen_stress
sign ratio        [.02,.1]  [.1,.2]       [.1,.3] [.2,.4]
  ↓
balanced (1:1)    ..........  ..........  ..........  ..........
rep_biased (1:2)  ..........  ..........  ..........  ..........
                  各 10 worlds  per cell
```

- **水平方向**：interaction strength 递增 → 观察 attractor 分布的退化模式
- **垂直方向**：同一个 strength 区间下对比 balanced vs repression_biased → 量化 repression 效应
- **对角线方向**：baseline_balanced → chen_stress_balanced 的极端对比 → 全参数空间的稳定性跨度

这个设计覆盖了从"几乎必然稳定"到"大部分发散"的全谱段，使得 parameter-regime mapping 既有覆盖面又具有 interpretable 的梯度结构。

### 2.4 Attractor 分类体系

每 trajectory 基于三个维度：equilibrium detection（convergence）、stability analysis（divergence/collapse）、oscillation analysis（damping），最终归入 7 类：

| Type | 名称 | 判定条件 |
|------|------|----------|
| A | Stable equilibrium | converged, conv_time ≤ 200 |
| B | Slow convergence | converged, conv_time > 200 |
| C | Damped oscillation | oscillatory, avg_damping > 0.05 |
| D | Sustained oscillation | oscillatory, avg_damping ≤ 0.05 |
| E | Runaway divergence | max expression ≥ 1e3 |
| F | Numerical collapse | final mean < 1e-3 |
| G | Others (unclassified) | 防御性兜底 |

每个 world 的 `primary_attractor` 由其 5 个 cell 的 majority vote 决定。

### 2.5 扰动分析

选取 5 种 attractor type 各 1-2 个 world（共 9 个），在 t=500 对其原始轨迹中表达量最高的基因执行**瞬时 knockdown**（该基因表达量设为 0，一个 timestep），然后让系统自然恢复。扰动后轨迹分析使用相同的 equilibrium / stability / oscillation 检测 pipeline。

| # | world_id | Type | regime | Knocked-down gene |
|---|----------|------|--------|-------------------|
| 1 | stress_balanced_t005 | D | stress | 15 |
| 2 | chen_stress_balanced_t008 | C | chen_stress | 4 |
| 3 | chen_stress_repression_biased_t001 | C | chen_stress | 15 |
| 4 | stress_repression_biased_t006 | B | stress | 10 |
| 5 | chen_moderate_balanced_t007 | B | chen_moderate | 6 |
| 6 | stress_repression_biased_t007 | E | stress | 10 |
| 7 | stress_balanced_t000 | E | stress | 15 |
| 8 | chen_stress_balanced_t005 | A | chen_stress | 15 |
| 9 | chen_stress_repression_biased_t007 | A | chen_stress | 3 |

---

## 3. 结果

### 3.1 Stability Geometry — 哪些 parameter geometry 更容易稳定？

**全局 attractor 分布**（80 worlds）：

| Type | 数量 | 占比 |
|------|------|------|
| A — Stable equilibrium | 57 | 71.3% |
| B — Slow convergence | 9 | 11.3% |
| C — Damped oscillation | 2 | 2.5% |
| D — Sustained oscillation | 1 | 1.3% |
| E — Runaway divergence | 11 | 13.8% |

**按 combo 的 attractor 分布**：

| combo_key | A | B | C | D | E | conv% | clip% | mean ρ(A) |
|-----------|---|---|---|---|---|-------|-------|-----------|
| baseline_balanced | **10** | 0 | 0 | 0 | 0 | 100% | 10% | 0.871 |
| baseline_repression_biased | **10** | 0 | 0 | 0 | 0 | 100% | 0% | 0.877 |
| chen_moderate_repr_biased | **10** | 0 | 0 | 0 | 0 | 100% | 100% | 0.969 |
| chen_moderate_balanced | **6** | 4 | 0 | 0 | 0 | 100% | 90% | 0.972 |
| stress_repression_biased | **7** | 1 | 0 | 0 | 2 | 80% | 100% | 1.035 |
| stress_balanced | 3 | **4** | 0 | 1 | 2 | 70% | 100% | 1.025 |
| chen_stress_repr_biased | **8** | 0 | 1 | 0 | 1 | 90% | 100% | 1.138 |
| chen_stress_balanced | 3 | 0 | 1 | 0 | **6** | 40% | 100% | 1.115 |

**结论**：interaction strength 是主导 regime 切换的轴。a_max ≤ 0.1（baseline）时无论 sign ratio 如何，100% Type A。当 a_max ≥ 0.3（chen_stress balanced）时，仅 30% 保持 Type A，60% 发散。

### 3.2 Oscillation Emergence — linear GRN 在什么条件下出现振荡？

- **Type C（衰减振荡）**：仅 2 例，均出现在 `chen_stress` regimes（a_max=0.4），ρ(A) = 1.030 和 1.083
- **Type D（持续振荡）**：仅 1 例（`stress_balanced_t005`），ρ(A)=1.083。该 world 在 4 个 strength regime 中仅在 stress 层级出现持续振荡，baseline 和 chen_moderate 下均为 Type A

**机制解释**：δ=0.2（self-retention 0.8）提供了时间持续性。当 interaction strength 足够大以对抗 decay damping 且 ρ(A) 接近 1 时，系统处于临界区域，signed 调控的正负反馈回路可维持振荡行为。由于是 linear GRN，振荡的存在完全来自 signed 调控的回路效应 + abundance clipping，而非非线性 Hill 函数。

### 3.3 Repression Effect — repression-biased 是否更稳定？

| 指标 | balanced (40 worlds) | repression_biased (40 worlds) |
|------|---------------------|-------------------------------|
| Type A 占比 | 22 (55.0%) | **35 (87.5%)** |
| Type E 占比 | 8 (20.0%) | 3 (7.5%) |
| Type B/C/D 占比 | 10 (25.0%) | 2 (5.0%) |

在高 stress regime 下对比尤为显著：

| regime | balanced | repression_biased |
|--------|----------|-------------------|
| chen_stress | A=3, E=6, C=1 | **A=8**, E=1, C=1 |
| stress | A=3, B=4, D=1, E=2 | **A=7**, B=1, E=2 |

**结论**：repression-biased 系统显著更稳定。在等价 interaction strength 级别下，rb 的 Type A 率比 balanced 高出 30+ 个百分点（87.5% vs 55%），发散率不到后者的一半。机制上，repressive 调控充当了"动态阻尼器"——负反馈回路将表达限制在 bounded region，防止 positive feedback runaway。

### 3.4 Spectral Radius Correspondence — ρ(A) 的预测能力取决于 clipping regime

| Attractor | n | mean ρ(A) | min ρ(A) | max ρ(A) | clip |
|-----------|----|-----------|----------|----------|------|
| Type A（non-clipping） | 20 | **0.877** | 0.800 | 0.962 | 0/20 |
| Type A（clipping） | 37 | 1.024 | 0.800 | 1.201 | 37/37 |
| Type B | 9 | 1.033 | 0.958 | 1.113 | 9/9 |
| Type C | 2 | 1.056 | 1.030 | 1.083 | 2/2 |
| Type D | 1 | 1.083 | — | — | 1/1 |
| Type E | 11 | 1.106 | 0.995 | 1.216 | 11/11 |

**Non-clipping 世界（baseline combo，20/20）**：全部 ρ < 1，max = 0.962，且全部为 Type A。ρ(A) < 1 在此区间是稳定的完美判据——无一反例。

**Clipping-dominated 世界（其余 60 个 world）**：ρ(A) 的预测力系统性失效。37 个 clipping Type A 中 20 个（54%）ρ > 1（最高 1.201），按线性判据应发散却稳定；反过来，1 个 Type E（`stress_repression_biased_t007`，ρ = 0.995）ρ < 1 却发散。ρ(A) 的范围在 Type A clip（0.800–1.201）与 Type E（0.995–1.216）之间大面积重叠，对单个 world 几乎无分类能力。

**原因**：clipping（ReLU 钳制至 0）引入了一个非线性边界条件。ρ(A) 描述的是 unbounded linear core 的扩张/收缩倾向，而 clipping 可以在 ρ > 1 时充当 nonlinear damping 将系统稳定化，也可以在 ρ 接近 1 时通过 clipping-and-rebound 诱发振荡甚至发散。因此，**ρ(A) 的有效性严格取决于系统是否运行在非 clipping 区域**：一旦进入 clipping-dominated regime，ρ 不能被直接解读为稳定判据。

### 3.5 Perturbation Recovery — 哪些 regimes 更容易恢复？

**一句话总结**：对表达量最高的基因在单一时步进行"表达量设为 0"的瞬时 knockdown（t=500）后，系统的 attractor 行为类别基本不变——是什么类型，扰动后还是什么类型。

扰动后的 attractor 判定方式：以扰动时刻（t=500）作为分析的"初始时刻"（即只看扰动后的轨迹段 X[500:]），使用与原始分析完全相同的 equilibrium / stability / oscillation 检测 pipeline 重新分类。因此，`perturbed_attractor` 反映的是扰动后系统的独立动力学倾向，而非与原轨迹的差异。

| world_id | 原类型 | 扰动后 | recovery_time | 注释 |
|----------|--------|--------|---------------|------|
| chen_stress_repr_biased_t007 | A | A | **545** | 稳定恢复 |
| chen_stress_balanced_t005 | A | A | **546** | 稳定恢复 |
| chen_moderate_balanced_t007 | B | **A** | **537** | 稳定恢复 |
| stress_repression_biased_t006 | B | **A** | **558** | 稳定恢复 |
| chen_stress_balanced_t008 | C | C | — | 保持衰减振荡 |
| chen_stress_repr_biased_t001 | C | C | — | 保持衰减振荡 |
| stress_balanced_t005 | D | D | — | 保持持续振荡 |
| stress_repression_biased_t007 | E | E | — | 继续发散 |
| stress_balanced_t000 | E | E | — | 继续发散 |

**发现**：
- A/B 类型世界（原已收敛者）扰动后全部恢复，recovery_time 在 537-558（绝对时间），即扰动后约 40-60 步重新达到均衡。其中 2 个原分类为 Type B 的 world 扰动后被重新判定为 Type A——这是因为扰动后轨迹从"被扰动打乱的状态"起步，经历了与扰动前不同的瞬态过程，最终同样恢复到均衡状态，且收敛速度较快（≤200 步），因此被分类为 Type A。
- Type C/D/E world 保持了原有动力学行为类别——扰动未改变 attractor regime。

**5 张扰动前后轨迹图**（`figures/perturbation/perturbation_type_[a-e].png`）提供了视觉证据：
- **Type A**：在 `chen_stress_balanced_t005` 中，扰动前后轨迹完全一致（系统对基因 15 免疫）；在 `chen_stress_repr_biased_t007` 中，扰动触发瞬态尖峰但 50 步内完全恢复
- **Type B**：2 个原 B world 扰动后均被重新判为 A，但图中可见扰动瞬间的瞬态扰动比 Type A 更明显（高幅尖峰）
- **Type C**：扰动后出现一轮新的衰减振荡，但最终回到原稳态（attractor basin 稳定）
- **Type D**：扰动前后轨迹**完全相同**（基因 15 显然对该持续振荡非必要）——可作为后续 "essential gene for oscillation" 研究的负样本
- **Type E**：repression_biased 扰动（`stress_repression_biased_t007`）将发散量级从 1.4×10⁸ 降低到 3.8×10⁷（约 1 阶），但 balanced 扰动（`stress_balanced_t000`）几乎无效（1.3×10²⁵ vs 1.2×10²⁵）——提示**某些 high-density 发散系统的 knockdown 难以缓解**，可能因为发散由多基因共同驱动

### 3.6 Abundance Compatibility — non-negativity constraint 显著影响吗？

| regime | clipping-dominated worlds |
|--------|---------------------------|
| baseline | 5% (1/20) |
| chen_moderate | 95% (19/20) |
| chen_stress | 100% (20/20) |
| stress | 100% (20/20) |

"clipping-dominated" 是文档（02_Output_Spec.md §6.1）要求的输出项，但未规定精确定义。本 run 使用的实现判据：clip 事件总数 > total_steps × G × 0.1（G=20, T=1000 即 clip 次数 > 2000，等价于每基因每步平均 clip 超过 0.1 次）。

[图: figures/clipping_dominated_by_regime.png](file:///home/zhanghl/projects/ddc_github/runs_lite/run_002/figures/clipping_dominated_by_regime.png) 进一步按 sign_ratio 拆分，可看到 8 组细分：
- baseline：balanced 10% (1/10) vs repression_biased 0% (0/10) —— **balanced 反而略高**
- chen_moderate：balanced 90% (9/10) vs repression_biased 100% (10/10) —— 两者都已接近饱和
- stress / chen_stress：两者均 100%

**结论**：当 interaction strength 达到 chen_moderate 及以上（a_max ≥ 0.2）时，几乎所有 world 都进入 **clipping-dominated regime**。linear eigenvalue theory（如 ρ(A)）在此区间仅提供方向性参考，不能替代 nonlinear clipped dynamics 的完整解释。系统处于 "piecewise-linear" 状态：大部分时间在 linear 区域运行，但频繁触发 zero-boundary（ReLU clipping），这实质上等价于一种 nonlinear constraint。

需要注意：上图的 8 组数据呈现"0% / 10% / 90% / 100%"的极端分布，提示**对于 clipping-dominated 的判定条件**过于严苛——当前的 0.1 阈值把低 strength 端几乎全归为非 clipping-dominated、把高 strength 端几乎全归为 clipping-dominated，导致中间的差异被压扁。未来可以考虑把阈值降低（例如 0.05 或基于 percentile 自适应），以便更细致地区分不同 regime 间的 clipping 行为差异。

### 3.7 Topology-Strength 交互

[图: figures/topo_strength_interaction.png](file:///home/zhanghl/projects/ddc_github/runs_lite/run_002/figures/topo_strength_interaction.png) 展示了 20 个独立 topology 在 4 个 strength regime 下的 attractor type 分布。横轴从左到右为 baseline → stress → chen_moderate → chen_stress，纵轴上半 10 行是 10 个独立采样的 balanced topology（topo0 - topo9），下半 10 行是 10 个独立采样的 repression_biased topology。颜色：绿=A、蓝=B、橙=C、红=D、紫=E。

- **baseline 列（a ∈ [0.02, 0.1]）**：两组全部 20 个 topology 均为 Type A
- **chen_moderate 列（a ∈ [0.1, 0.2]）**：仅 balanced 组 4 个 topology（topo0, topo2, topo6, topo7）出现 Type B，rb 组仍全部 Type A
- **stress 列（a ∈ [0.1, 0.3]）**：balanced 组 7/10 切换为非 A（5 个 B、1 个 D、2 个 E，其中 topo5 的 Type D 是全图唯一的持续振荡）；rb 组仅 2/10 出现非 A（topo6=B, topo7=E）
- **chen_stress 列（a ∈ [0.2, 0.4]）**：balanced 组 6 个 topology 落入 Type E（topo0, topo2, topo4, topo6, topo7, topo9）+ 1 个 Type C（topo8）；rb 组仅 1 例 Type E（topo8）+ 1 例 Type C（topo1），其余 8 个保持 Type A

两组全程 Type A 的 topology 数：balanced 2 个（topo1, topo3）vs repression_biased 6 个（topo0, topo2, topo3, topo4, topo5, topo9）。这一对比直接量化了 §3.3 的核心结论：repression 偏置在中等与最强强度 regime 下都提供了显著的稳定性保护。stability 是个多维度问题：topology-specific 的稳定性与 sign_ratio 偏置共同决定。

### 3.8 Clipping Frequency — clipping 在不同 attractor / sign 下的实际发生率

[图: figures/clipping_frequency.png](file:///home/zhanghl/projects/ddc_github/runs_lite/run_002/figures/clipping_frequency.png) 展示了不同 attractor type 与 sign_ratio 下的 clipping frequency（0-0.5 范围）分布。结合 §3.6 的"是否 clipping-dominated"二元判据，这里进一步量化了**实际发生率**：

- **Type A**：balanced 的 clipping 频率中位数 ~0.125（很多 world clip 比例很低），repression_biased 中位数 ~0.3（普遍更高）。这是一个反直觉的现象——尽管 rb 总体更稳定（更多 Type A），但其进入 Type A 路径上需要更频繁地触发 ReLU clipping。
- **Type B**：balanced 中位数 ~0.25，repression_biased 几乎全部压缩在 0.3 附近（低方差）。
- **Type C/D**：clipping 频率均集中在 0.22-0.35 区间，方差小。
- **Type E**：repression_biased 中位数高达 0.45，最大值 0.5（全图最高）；balanced 中位数 0.3，但有一个低 outlier（~0.15）。

**解读**：clipping 越频繁反而越可能与发散相关——Type E 的 clipping 频率中位数最高，提示过度 clipping 反映系统在 "linear growth ↔ zero-boundary" 间的快速震荡，本身是不稳定信号。

### 3.9 Topology 拓扑特征 — hub-likeness 与 stability 的关系

[图: figures/degree_distribution_vs_stability.png](file:///home/zhanghl/projects/ddc_github/runs_lite/run_002/figures/degree_distribution_vs_stability.png) 展示了 80 个 world 按 in-degree 标准差（衡量 hub-likeness）分四分位后的 attractor 分布：

| Hub-likeness | Type A | Type B | Type C | Type D | Type E |
|---|---|---|---|---|---|
| Q1 (低) | 13 | 3 | 0 | 0 | 4 |
| Q2 | 14 | 0 | 0 | 0 | 2 |
| Q3 (中高) | 12 | 6 | 0 | 1 | 5 |
| Q4 (高) | **18** | 0 | **2** | 0 | 0 |

**反直觉发现**：hub-likeness 最高的 Q4 反而有最多的 Type A（18 个）且 **Type E 完全消失**。这与"hub 越多越不稳定"的直觉相反——可能的原因是：sparse G=20 拓扑中 hub 基因提供了 dominant sink/source，迫使系统偏向 single-equilibrium 结构，从而压制了发散行为。Q3（中等 hub-likeness）反而是 diversity 最高的区间（同时出现 A/B/D/E）。这提示**stability 与拓扑特征的关系不是单调的**——可能存在一个"最优 hub-likeness"区间，过低（Q1）和中等（Q3）都更容易出现非 A。

### 3.10 收敛时间 CDF — sign ratio 对 convergence speed 的量化对比

[图: figures/convergence_time.png](file:///home/zhanghl/projects/ddc_github/runs_lite/run_002/figures/convergence_time.png) 是 Type A + Type B 子集的收敛时间 CDF（n=31 balanced, n=37 repression_biased）：

- **中位收敛时间**：repression_biased ~65 步，balanced ~85 步
- **100 步时累计收敛率**：repression_biased 90% vs balanced 60%
- **收敛时间尾部分布**：repression_biased 在 200 步附近基本完成（1 个 ~925 步的极端 outlier 除外）；balanced 的尾一直拖到 ~780 步

这量化了 §3.4 表格中 Type A 收敛时间（mean 69）vs Type B（mean 367）的差异——repression_biased 不仅更少落入 B/C/D/E，而且**即便落入 Type A 也收敛更快**。sign ratio 对 convergence speed 也有系统性影响。

---

## 4. 与 Chen (2019) 的对照

| 维度 | Chen (2019) | run_002 | 一致性 |
|------|-------------|---------|--------|
| Network size | N=500 (sim), N=356-746 (GRN) | G=20 | 缩小（可解释性优先）|
| Edge density | r=0.1, r=0.076, r=0.031 | r=0.1 | ✅ 一致 |
| Regulation sign ratio | ~1:1（yeast 0.504:0.496）| balanced 1:1 | ✅ 一致 |
| Interaction strength | a∈[0.1275,0.1675] (sim) | chen_moderate [0.1,0.2] | ✅ 可比 |
| Time formalism | Continuous-time (ODE) | Discrete-time | ⚠ 不同 |
| Stability criterion | max Re(λ(M)) < 0 | ρ(A) < 1 | 对应但不对等 |
| Runtime | Linear ODE | Piecewise-linear (ReLU clip) | ⚠ 不同 |

**关键差异**：Chen 使用 continuous-time local linearization，其 Jacobian stability 判据为 `max Re(λ(M)) < 0`。run_002 使用 discrete-time abundance-compatible runtime，对应判据为 `ρ(A) < 1`。此外，Chen 模型中不存在 non-negativity clipping——其 abundance 由 ODE 自然维持在 positive region。run_002 的 ReLU clipping 在高 interaction 时成为主导因素，使系统本质上成为 piecewise-linear system。

---

## 5. 讨论

### 5.1 主要发现

1. **Interaction strength 是 regime 切换的主控轴**。在 δ 固定时，仅需 a_max 从 0.1 变化到 0.4，系统就从 100% stable equilibrium 切换到 60% divergent（balanced）或 10% divergent + 10% oscillatory（rb）。

2. **Repression is a universal stabilizer**。在同等 interaction strength 级别下，repression-biased 系统的 Type A 率全面高于 balanced 系统。负反馈充当了 dynamical damping 机制。

3. **Linear GRN 的动力学出人意料地丰富**。尽管没有 nonlinear Hill 函数、没有 protein layer、没有 delay，系统仍然能产生慢收敛、衰减振荡和持续振荡——这些行为完全来自 signed regulation + abundance clipping 的组合。

4. **ρ(A) 的诊断价值依赖 clipping 状态**。在 non-clipping 世界中 ρ < 1 ⇔ 稳定，线性理论完美成立（20/20 无反例）；在 clipping-dominated 世界中线性谱判据系统性失效——Type A clip 中 54% ρ > 1，Type E 中也有 ρ < 1 的个案。换言之，ρ(A) 仅在系统运行于 linear regime 时才是可靠的稳定指标；在 piecewise-linear regime（clipping-dominated）中，ρ 与 attractor type 几乎脱钩。

### 5.2 方法论反思

- **δ 固定为 0.2**：当前无法区分 "decay 效应" 和 "interaction 效应"。`decay vs interaction strength` 图在文档（02_Output_Spec.md）中被要求绘制，但因本 run 的 δ 为定值（仅 0.2），该图暂未生成——需要未来 δ-scan run 来补齐。
- **5 cells/world 的 majority vote**：偶尔掩盖 cell-to-cell 差异。在 `stress_repression_biased_t006` 中，cell 0 为 Type A，但 5 个 cell 的 majority 为 Type B。
- **扰动分析规模有限**（9 cases）：未来可扩展到每种 combo 各选代表。
- **Type F/G 未出现**：不代表 future runs（更大 N、更 extreme regimes）不会出现。

### 5.3 后续方向

- δ-scan run：探索 decay strength 与 interaction strength 的交互效应
- Model A → Model B 比较：引入 explicit X/P layering + delayed protein propagation
- 扩大扰动分析：覆盖更多 regime × sign 组合，引入不同类型的 perturbation（activation pulse, multi-gene knockdown）
- Topology correlation：哪种拓扑特征（diameter, motif count, degree variance）最能预测 attractor behavior

---

## 6. 结论

Run 002 成功建立了 **Chen-style sparse Classical Linear GRN 的 parameter-regime dynamical map**：

- **80 worlds, 400 trajectories**，覆盖 4 strength regimes × 2 sign ratios
- 发现了 5 种 distinct attractor regimes（A-E），各类型的行为清晰可区分
- 验证了 4 条核心假设：decay-dominated 稳定、repression 增强稳定性、linear GRN 可振荡、parameter geometry 主导 regime
- 扰动分析表明 convergent systems（A/B）对瞬时 knockdown 具有弹性恢复能力
- 满足文档需求的输出全部到位（详见 Appendix A），唯一的已知缺口是 `decay vs interaction strength` 图（因 δ 固定为定值，见 §5.2）

**One-line summary**: *Run 002 establishes the interpretable discrete-time Chen-style sparse Classical Linear GRN dynamical foundation — the reference coordinate system for the DDC representation curriculum.*

---

## Appendix A: 输出文件清单

| 目录 | 文件 | 数量 | 说明 |
|------|------|------|------|
| `tables/` | `world_summary.tsv` | 80 rows, 22 cols | per-world summary |
| `tables/` | `regime_summary.tsv` | 8 combos, 13 cols | per-combo aggregation |
| `analysis/` | `all_world_results.json` | 80 entries | 完整 world 数据（含 eigenvalues） |
| `analysis/` | `aggregated_analysis.json` | 8 combos | combo-level aggregation |
| `analysis/` | `topo_strength_analysis.json` | 20 topologies | topology × strength 交互 |
| `analysis/` | `perturbation_summary.json` | 9 entries | 扰动前后对比 |
| `world_metadata/` | `*.json` | 80 files | per-world topology + parameters |
| `trajectories/` | `*.pt` | 400 files | full X(t) per cell |
| `perturbations/` | `*_perturb.pt` | 9 files | 扰动后轨迹 |
| `figures/` | `canonical_trajectories.png` | 1 file | 5 类 attractor 示例 |
| `figures/` | `spectral_radius.png` | 1 file | ρ(A) per regime |
| `figures/` | `convergence_time.png` | 1 file | 收敛时间分布 |
| `figures/` | `clipping_frequency.png` | 1 file | 裁剪频率分布 |
| `figures/` | `topo_strength_interaction.png` | 1 file | topology × strength 热力图 |
| `figures/` | `repression_ratio_vs_stability.png` | 1 file | sign ratio 效应 |
| `figures/` | `degree_distribution_vs_stability.png` | 1 file | degree 方差效应 |
| `figures/` | `clipping_dominated_by_regime.png` | 1 file | clipping-dominated 占比按 regime × sign 拆分 |
| `figures/traj/` | `all_trajectories_*.png` | 6 files | 全部 world 的 cell 0 轨迹 |
| `figures/perturbation/` | `perturbation_type_[a-e].png` | 5 files | 扰动前后对比 |
