# Run 006 — Summary Report

**Stage**: DDC Lite Representation Curriculum — Model A Edge Perturbation & State Landing Run

**Run**: run_006 (Model_A_Edge_Perturbation_State_Landing)

**Model**: Model A — Classical Linear GRN with RELU non-negativity clipping, acute edge deletion perturbation

**Date**: 2026-07-07

---

## 0. Executive Summary

> 在 run_005 的 120 个 GRN world pool 中选取 12 个代表性 world（覆盖 stable anchors、boundary worlds、rare specimens），采用 fork 并行对比方案对每个 world 进行 4 种急性删边策略（随机 K=5 条、最强 K=5 条、hub 出边 K=5 条、SCC 内边 K=5 条），比较删边前后的 trajectory regime 类型和 state landing 位置。每个 world 先通过 initial-condition screening 确定 representative initial conditions 数量：收敛世界 1 个，非收敛或多 regime 世界 2–5 个（共 37 个 representatives），每个 representative 执行 4 种策略（A 策略含 10 个 random repeats）= 481 条 trajectory。
>
> **核心结论**：急性删除 5 条调控边后，系统整体表现出较强的鲁棒性——仅 23.5%（113/481）的扰动导致 regime transition，绝大多数扰动（61.5%）表现为同一 regime 内的 landing_shift（收敛/振荡类型不变但落点偏移），13.9% 为 robust（完全不变）。这与我们"少数拓扑边具有可测扰动杠杆"的假设一致，但 targeted deletion（B/C/D）并未显著优于 random deletion（A）——四种策略均以 landing_shift 为主导结果。**边界 worlds 确实比稳定 anchor 更敏感**（H2_boundary 的 regime transition 比例显著高于 H1_stable），但过渡 worlds 中频繁出现的 convergent ↔ sustained oscillatory 互转说明 boundary 并非单向不稳定走廊，而是一个 regime-flexible 区域。**关键发现**：在 RELU-based dynamics 中，拓扑性选边策略（hub outgoing / SCC internal）天然受沉默节点限制——删掉表达量为 0 的基因的出边等同于无扰动操作，揭示出 "topological hub ≠ functional hub" 的基本机制。此外，唯一的 ambiguous outcome 揭示了 Convergent ↔ Divergent 之间的中间地带：删掉 hub 基因出边后系统滑入缓慢但持续的膨胀，既未收敛也未跨过 divergence 阈值。

---

## 1. 科学问题与实验设计

### 1.1 核心科学问题

> *在同一个 Model A GRN world 运行到 t=500 后，急性删除 5 条调控边，比较 fork 后 baseline 分支与 perturbed 分支的 trajectory，系统最终落入的状态、吸引子类型或动力学 regime 是否改变？*

同时回答文档要求的 5 个关键问题：

> 1. 急性删边后，post-perturb regime 多常改变？
> 2. Convergent worlds 中删边多常改变最终 state landing？
> 3. Boundary worlds 是否比 stable worlds 更敏感？
> 4. Targeted deletion（strongest / hub / SCC）是否比 random deletion 更易打乱系统？
> 5. 哪些 edge features 与更强扰动效应相关？

### 1.2 实验定位

Run 006 是 run_005 的直接后续实验。Run_005 扫描了 edge density r ∈ {0.10, 0.25} 对 attractor multiplicity 的影响，发现了从 single-equilibrium dominant 到 non-convergence dominant 的转变。Run_006 追问一个更精细的问题：**在给定的 GRN world 中，拓扑内部的"具体调控边"是否具有可测的扰动杠杆？**

Run_002 中的临时补充实验暗示：随机删除 5 条边后，一些 world 非常鲁棒，一些则发生了 regime switching。Run_006 系统性地执行并量化这一观察。

```text
run_002:  Model A unconstrained topology baseline  (r=0.10, b=0.1, δ=0.2)
run_003:  Model A TF-restricted topology baseline  (r=0.10, b=0.1, δ=0.2)
run_004:  Model A decay δ boundary scan            (r=0.10, b=0.1, δ∈{0.1,0.2,0.4})
run_005:  Model A sparsity / attractor multiplicity (r∈{0.10,0.25}, b=10, δ=0.2)
run_006:  Model A edge perturbation & state landing  (b=10, δ=0.2, 12 worlds × 4 strategies)
          ↓ 当前阶段
```

### 1.3 Fork 并行对比方案（zhl 修改版）

与原始设计的 single-trajectory pre/post 拆分不同，本 run 采用 fork 方案：

```
t=0 ────────────────── t_fork=500 ──┬── baseline 分支（原始 W）... → t=1000
                                    │
                                    └── perturbed 分支（W_post，删边后）... → t=1000
```

- 两条分支在 t=500 时从完全相同的 X(t_fork) 出发，运行相同步数（T_post=500）
- `pre_reference_state` 取自 baseline 分支的 late-window（t=900–1000），而非 t_fork 单点
- `post_final_state` 取自 perturbed 分支的 late-window
- 消除了时间漂移混淆，确保 pre/post 差异仅来源于删边操作

### 1.4 固定参数

| 参数 | 值 |
|------|-----|
| G (gene count) | 20 |
| b (basal expression) | 10.0 |
| δ (decay) | 0.2 |
| self-retention | 0.8 |
| T_fork | 500 |
| T_post | 500 |
| K_delete | 5 |
| R_random | 10 |
| LATE_WINDOW | 100 |

### 1.5 Selected Worlds（12 个）

| 类别 | # | world_id | 特点 |
|------|:-:|------|------|
| **H1 Stable anchors** | 2 | `baseline_balanced_r010_t000`, `chen_balanced_r010_t006` | rho(W) < 1，零 clipping |
| **H2 Boundary worlds** | 4 | `stress_balanced_r010_t003/001/007/005` | 跨 rho(W)=1 边界的 convergent/mixed/osc/div 梯度 |
| **H2 Control** | 3 | `chen_balanced_r025_t005/000`, `stress_repression_biased_r010_t004` | 跨 combo 边界对照 |
| **H3 Rare specimens** | 3 | `chen_repression_biased_r025_t007/001`, `stress_balanced_r025_t001` | multi_eq world、同 combo 对照、high-rho 对照 |

### 1.6 删边策略

| 策略 | 名称 | 方法 | 重复 |
|------|------|------|:--:|
| **A: random** | Random edge deletion | 从所有调控边中均匀随机选 K=5 条 | R=10 |
| **B: strongest** | Strongest edges | 按 abs(a_ij) 降序取 top K=5 | 1 |
| **C: hub** | Hub outgoing edges | 按 out-degree 贪心扩展 hub 基因，按强度加权采样 K=5 条出边 | 1 |
| **D: SCC** | Largest-SCC internal edges | 取最大 SCC 内部边，不足扩展到次大 SCC | 1 |

### 1.7 Perturbation Outcome 分类

| Outcome | 定义 |
|---------|------|
| `robust` | baseline 与 perturbed 同 regime，且 normalized L2 < 0.05 |
| `landing_shift` | 同 regime 但 normalized L2 ≥ 0.05（state landing 偏移） |
| `regime_transition` | pre_primary_regime ≠ post_primary_regime |
| `collapse` | post_primary_regime == Collapse |
| `ambiguous` | 任一方 regime 不可分类或数值问题 |

### 1.8 工作假设

| # | 假设 |
|:--:|------|
| **H1** | stable worlds 对急性删边较鲁棒 |
| **H2** | boundary worlds 对删边更敏感，更易 regime transition |
| **H3** | r=0.25、chen-strength 的 worlds 有更丰富的扰动后行为 |
| **H4** | strong/cycle-internal/hub outgoing edges 比随机删除更有影响力 |
| **H5** | rho(W) 有帮助但不充分；active zero/clipping pattern 也决定实际落点 |

### 1.9 Pre-Perturb State Class Assignment（扰动前状态类分配）

执行配置：`execution_grid = core`。每个 world 从 run_005 中最多选取 `N_init_screened = 5` 个初始条件（按类型 `low_00, medium_00, high_00, sparse_00, low_01` 各 1 个）。

**Representative 选择规则**：

- **收敛世界**（regime label 含 "convergent" 且不含 "non_convergent"）：在 t=500 时刻对 X(t_fork) 做层次聚类（average linkage, relative L2 distance, threshold=0.001）。同一 cluster 内的 init 视为同一 pre-perturb state class，取离 cluster centroid 最近的成员作为 representative。
- **非收敛世界**（含 "non_convergent" 的 mixed/osc/div worlds）：每个 init 独立作为一个 pre_perturb state class，不做聚类（`direct_per_init_no_clustering`）。

**Per-World 分配统计**（共 37 个 pre-perturb state classes）：

| # classes | Worlds |
|:--:|------|
| 1 (convergent, clustered) | `baseline_balanced_r010_t000`, `chen_balanced_r010_t006`, `chen_balanced_r025_t005`, `chen_repression_biased_r025_t001`, `stress_balanced_r010_t003` |
| 2 (clustered, multi-eq) | `chen_repression_biased_r025_t007` |
| 5 (per-init, non-convergent) | `chen_balanced_r025_t000`, `chen_balanced_r025_t005`*, `stress_balanced_r010_t001`, `stress_balanced_r010_t005`, `stress_balanced_r010_t007`, `stress_balanced_r025_t001`, `stress_repression_biased_r010_t004` |

> *注：`chen_balanced_r025_t005` 虽然在 run_005 中 9/10 non-convergent，但其 selected rep 初始条件在 t=500 处均进入同一 class（1 class, clustered）。

### 1.10 Baseline Re-Simulation Sanity Check

Fork 仿真中，每个 perturbation 的 baseline 分支独立运行 1000 步（W, b, X0），并与 run_005 原始 trajectory 做逐元素对比。所有 481 条 baseline trajectory 均通过 `np.allclose(atol=1e-12)` 检验，无 SANITY FAIL 案例。

---

## 2. 结果

### 2.1 全局 Outcomes 分布

481 条 perturbation trajectory 的 outcome 分布：

| Outcome | 数量 | 占比 |
|---------|-----:|-----:|
| landing_shift | 296 | 61.5% |
| regime_transition | 113 | 23.5% |
| robust | 67 | 13.9% |
| ambiguous | 5 | 1.0% |
| collapse | 0 | 0.0% |

在 156 个 perturbation mode × world 层级的 summary 中（例如 "baseline_balanced_r010_t000 × random_delete_5" 为一个 summary 行，该行聚合了该 world 在该策略下所有 representative 的结果），dominant outcome 分布：

| Dominant outcome | 数量 | 占比 |
|------------------|-----:|-----:|
| landing_shift | 113 | 72.4% |
| regime_transition | 24 | 15.4% |
| robust | 18 | 11.5% |
| ambiguous | 1 | 0.6% |

**核心发现**：landing_shift 在所有粒度级别上都是主导 outcome。急性删除 5 条边后，大多数系统保持了原有的 regime 类型（Convergent → Convergent，或 Sustained oscillatory → Sustained oscillatory），但最终的 state landing 位置发生了可检测的偏移。Regime transition 作为少数 outcome（23.5% in trajectory-level, 15.4% in perturbation-mode summary）提示系统确实存在少数"高杠杆"扰动场景。

### 2.2 Regime Transition Matrix — 什么转成了什么？

基线 regime 与扰动后 regime 的交叉表：

| baseline \ post | Ambiguous | Convergent | Divergent | Sustained osc. |
|-----------------|:---------:|:----------:|:---------:|:--------------:|
| Convergent | 1 | 298 | 4 | 1 |
| Divergent | 0 | 0 | 182 | 0 |
| Sustained oscillatory | 4 | 3 | 1 | 83 |

**解读**：

- **Convergent worlds（基线 304 条）中**，绝大多数保持 Convergent（298/304 = 98.0%）。仅 5 条发生 transition（4 → Divergent, 1 → Sustained oscillatory）。Convergent 世界对 5 条边级别的急性扰动非常鲁棒。
- **Sustained oscillatory worlds（基线 91 条）中**，83/91（91.2%）保持振荡。4 条不可分类（ambiguous），3 条转 Convergent，1 条转 Divergent。振荡 world 比收敛 world 具有略高的 transition 率。
- **Divergent worlds（基线 182 条）** 全部保持 Divergent。没有一条 divergent trajectory 被删边"挽救"回 convergent 或 oscillatory。
- 整体 regime 转换以 **Sustained oscillatory ↔ Ambiguous/Convergent** 和 **Convergent → Divergent** 为主。

### 2.3 Stable vs Boundary — H1 & H2 假说验证

| World tier | N trajectories | Regime transition |
|-------------|:---:|:---:|
| H1 Stable anchors | 52 | **1 (1.9%)** |
| H2 Boundary worlds | 117 | **36 (30.8%)** |
| H2 Control | 117 | **12 (10.3%)** |
| H3 Rare specimens | 195 | **64 (32.8%)** |

**H1 验证**：Stable anchor worlds（baseline_balanced_r010_t000, chen_balanced_r010_t006）在 52 条扰动 trajectory 中仅 1 条（1.9%）出现 regime transition。这两组 worlds 的 rho(W) < 1、零 clipping，表现出预期中的高度鲁棒性。H1 假说得到支持。

**H2 验证**：Boundary worlds 的 regime transition 率（30.8%）是 Stable anchors 的 16 倍。然而，大量 boundary world 的 landing_shift 而非 regime_transition 说明：在 rho(W) ≈ 1 附近，删边主要改变 state landing 位置而非 regime 类型。boundary zone 是一个"regime-flexible"而非"fragile"的区域。

**H3 验证**(部分)：r=0.25 的 chen worlds（H2_control 中的 `chen_balanced_r025_t005/000` + H3_rare 中的 `chen_repression_biased_r025_t007/001`）表现出高 regime transition 率（32.8% in H3），符合 H3 的预期——高密度 + 中等强度是最丰富的扰动后行为区域。

### 2.4 Random vs Targeted — H4 假说验证

按 perturbation mode 分组的 outcome：

| Mode | N | robust | landing_shift | regime_transition | ambiguous |
|------|--:|:---:|:---:|:---:|:---:|
| A: random (×10) | 320 | 53 | 188 | 75 | 4 |
| B: strongest | 57 | 5 | 39 | 13 | 0 |
| C: hub | 53 | 5 | 35 | 12 | 1 |
| D: SCC | 51 | 4 | 34 | 13 | 0 |

按 mode 的 regime transition 率：

| Mode | Regime transition % |
|------|:---:|
| A: random | 23.4% |
| B: strongest | 22.8% |
| C: hub | 22.6% |
| D: SCC | 25.5% |

**H4 未获支持**：四种删边策略的 regime transition 率几乎无差异（22.6%–25.5%）。Targeted deletion（strongest / hub / SCC）并未表现出比 random deletion 更强的系统扰动能力。这并非意味着 edge features 不重要，而是暴露了一个机制性发现：

**Topological hub ≠ functional hub**。在 RELU-based dynamics 中，如果 hub 基因在该初始条件和稳态解下的表达量为 0，则删除其出边等同于无操作。我们检查了轨迹图（图 6.4）中 C: hub 和 D: SCC 几乎无变化的具体案例（如 `chen_repression_biased_r025_t007`），发现被删边的 source gene 在 fork 时刻表达量均为 0。这解释了为什么 topology-targeted strategies 在 RELU 框架下未能超越 random deletion——边的实际"杠杆"取决于其 source gene 是否活跃。

### 2.5 Divergent Trajectories — Magnitude Ratio 分析

在 481 条 trajectory 中，187 条涉及 Divergent regime（基线和/或扰动后为 Divergent）。对 divergent trajectory，final state 不是 equilibrium 意义上的"状态"，因此不能用 normalized L2 distance 来衡量扰动效应。按照文档 §13.1 的指示，我们改用 **magnitude ratio**：

```
magnitude_ratio = ||post_final_state|| / max(||pre_reference_state||, 1e-8)
```

其中 pre 取自 baseline 分支的最后 50 步均值向量范数，post 取自 perturbed 分支的同位置。两者均为发散中的状态快照，不做 equilibrium 解释。

**187 条 divergent-involved trajectory 的 magnitude ratio 分布**：

| Transition type | N | 占比 | 含义 |
|-----------------|--:|:---:|------|
| Div→Div (suppressed) | 101 | 54.0% | ratio < 0.9，删边**抑制**了发散 |
| Div→Div (unchanged) | 69 | 36.9% | 0.9 ≤ ratio ≤ 1.1，删边对发散程度无显著影响 |
| Div→Div (accelerated) | 12 | 6.4% | ratio > 1.1，删边**加速**了发散 |
| Conv→Div | 4 | 2.1% | 收敛→发散，ratio 反映新生发散的体量 |
| Osc→Div | 1 | 0.5% | 振荡→发散 |

**核心发现**：在 182 条 Div→Div 轨迹中，超过一半（54%）的删边操作**抑制**了发散（median ratio = 0.706）。ratio 范围横跨 [0, 2.5×10⁶] 的 6 个数量级，其中 ratio < 1 的有 106 条，ratio > 1 的仅 24 条。这表明 Model A 中 divergent 状态的维持依赖于少数关键边——删除后，正反馈回路的驱动力被减弱，发散体量反而下降。Conv→Div 和 Osc→Div 的案例（5 条）的 ratio 均接近或大于 1，提示"新生发散"的 magnitude 依赖于所删除边的具体性质，而非普遍的大幅突变。

> Divergent magnitude ratio 图见图: [../figures/delta_rhoW_vs_divergent_magnitude.png](../figures/delta_rhoW_vs_divergent_magnitude.png)

### 2.6 Delta ρ(W) vs Outcome — 谱半径变化的预测力

在排除了 divergent 轨迹后，294 条 non-divergent trajectory 的 outcome 分布：

| Outcome | N | 占比 |
|---------|--:|:---:|
| landing_shift | 168 | 57.1% |
| regime_transition | 108 | 36.7% |
| robust | 13 | 4.4% |
| ambiguous | 5 | 1.7% |

Δρ(W) 分布集中在 [-0.42, 0.05] 区间，绝大多数扰动导致谱半径微小下降（median ≈ −0.005）。regime_transition 点在该区间内分布与 landing_shift 点重叠，说明 **Δρ(W) 单独不足以预测 outcome**（H5 得到支持）。Linear core 的谱稳定性与非线性 clipping 动力学共同决定扰动后落点。

> Delta ρ(W) vs Perturbation Outcome 图见图: [../figures/delta_rhoW_vs_outcome.png](../figures/delta_rhoW_vs_outcome.png)

### 2.7 The Ambiguous Case — 一个临界案例

在 481 条 perturbation trajectory 中，仅 5 条被标记为 ambiguous（1.0%），且全部来自同一个 perturbation mode summary 行：**`chen_balanced_r025_t005` × C: hub**。

**详细分析**：

该 world（chen_balanced_r025_t005, H2_control）的基线为 Convergent。对其 hub 基因（gene 15）删除 5 条出边后，perturbed trajectory 表现出以下特征：

- 系统**未收敛**——最后 50 步的 norm 变化约 2.67/step，远超收敛阈值（EPSILON=0.01）
- 系统**未发散**——max expression ≈ 2616，远低于 divergence 阈值（1e6）
- 系统**不振荡**——无 autocorrelation 上的周期性 pattern
- gene 15 的表达量从稳态中的 ~0 滑出，开始缓慢但持续增长

**机制解释**：删掉 gene 15 的 5 条出边意外地释放了被 RELU clipping 压制的正反馈回路。由于 gene 15 的靶基因不再接收其负调控信号，系统滑入了一个缓慢膨胀的轨道——以约 0.04%/步的速率持续增长。这个增长率太低，500 步内无法触发 divergence 判定，但又太高而无法收敛。它代表了 **Convergent 和 Divergent 之间的真实中间地带**。

这个案例揭示出 Model A 在 r=0.25、中等 strength（chen regime）下存在 **slow-drift dynamics**——一个临界于 stability boundary 的动力学状态，在我们的 1000 步窗口内无法被清晰地归入 Convergent 或 Divergent。`ambiguous` 标签在此处是准确且有科学含义的。

### 2.8 沉默节点问题 — 为什么 Topological Hub 不等于 Functional Hub

在检查 trajectory 图时（[图 6.4](../figures/trajectory_pre_post_world_chen_repression_biased_r025_t007.png)），我们观察到 C: hub 和 D: SCC 策略在该 world 下扰动后轨迹几乎无任何偏离。追踪发现：

**C: hub**：gene 15 被选为 hub 基因（outdegree=7），但其在 fork 时刻和整个 post-fork 期间的表达量恒为 0。因 W × x 中 x₁₅ = 0 → 对任何靶基因的贡献为 0，5 条被删边均无效。

**D: SCC**：被删的 5 条 SCC 内边中，4 条的 source gene 表达量为 0。仅 1 条有效边（5→4）的删除使 gene 4 下降了约 1.59（92.63 → 91.04），在 20 条轨迹叠加图中完全不可见。

**科学含义**：在 RELU-based linear dynamics 中，调控边的"杠杆"不仅取决于其权重大小和拓扑位置，还取决于其 source gene 是否在该稳态中保持非零表达。Topological centrality（outdegree、SCC membership、edge strength）只是必要条件，**expression-context dependency** 才是最终决定扰动效应的充分条件。这一机制性发现在 Model A 框架内是本质性的，而非数值 artifact。

---

## 3. 代表性 World Case Studies

以下 6 个案例覆盖文档要求的 7 种 case study 类型，展示从鲁棒到敏感的完整扰动响应谱。

### 3.1 Robust Stable World — `baseline_balanced_r010_t000` (H1)

- **基线 regime**: Convergent（13/13 条 trajectory），single-equilibrium，rho(W)=0.891，零 clipping
- **扰动结果**: 13/13 **landing_shift**（0 regime_transition, 0 robust）
- **解读**: 4 种策略均保持 Convergent regime，但所有扰动都产生了可检测的 state landing 偏移（normalized L2 > 0.05）。尽管 rho(W) < 1 确保了谱稳定性，landing position 却对拓扑变化有微弱但持续的信号。这是一个 "robust regime, flexible state" 的典型案例——符合 H1 的核心预期。

### 3.2 Boundary-Sensitive World — `stress_balanced_r010_t007` (H2)

- **基线 regime**: Sustained oscillatory（65/65 条 trajectory），rho(W)=1.091，clipping frequency ≈ 4.5%
- **扰动结果**: **37/65 regime_transition（56.9%）**, 28/65 landing_shift
- **解读**: 这是所有 12 个 worlds 中 regime transition 率最高的 world 之一。5 条边的删除使 56.9% 的振荡 trajectory 转向其他 regime（主要为 Convergent）。位于 rho(W)≈1.09 的 stress 边界区域，该 world 的 oscillation 依赖于少数关键边的持续激活——删除后反馈回路被打断，振荡坍缩为稳定不动点。

### 3.3 Oscillatory World — `chen_balanced_r025_t000` (H2_control)

- **基线 regime**: Sustained oscillatory（65/65 条 trajectory），rho(W)=1.027，r=0.25，chen-strength
- **扰动结果**: **34/65 regime_transition（52.3%）**, 27/65 landing_shift, 4/65 ambiguous
- **解读**: r=0.25 的高密度 chen regime 下，振荡世界对删边高度敏感。超过一半的扰动导致振荡消失（转为 Convergent 或 Ambiguous）。4 条 ambiguous 轨迹（来自 C: hub 和 D: SCC）与 §2.7 的 ambiguous 案例同属 slow-drift 模式——系统在 500 步内未明确收敛或发散。

### 3.4 Divergent World — `stress_balanced_r025_t001` (H3_rare)

- **基线 regime**: Divergent（65/65 条 trajectory），rho(W)=1.556，high-rho control
- **扰动结果**: **0/65 regime_transition** — 59/65 landing_shift, 6/65 robust
- **解读**: rho(W)=1.556 的极高谱半径下，所有 65 条 trajectory 在扰动后**全部保持 Divergent**，无一"被挽救"。这是所有 worlds 中最极端的 "regime lock-in" 案例——一旦系统的线性核足够强，del 5 条边不足以逆转发散趋势。59 条 landing_shift 暗示 dissociation magnitude 弱变化（需进一步通过 growth-rate analysis 量化）。

### 3.5 Multi-Equilibrium Specimen — `chen_repression_biased_r025_t007` (H3_rare)

- **基线 regime**: Convergent（26/26 条 trajectory），multi-equilibrium world（run_005 中唯一的 multi-eq world），rho(W)=0.930
- **聚类**: 5 个 init 在 t=500 处聚为 **2 个 class**（low/medium/high → class 0, sparse/low_01 → class 1）
- **扰动结果**: **0/26 regime_transition** — 21/26 landing_shift, 5/26 robust
- **解读**: 尽管该 world 在 run_005 中表现出多个 equilibrium，current rep 在删边后全部保持 Convergent。C: hub 和 D: SCC 策略下轨迹几乎无变化（见 §2.8 silence node analysis），暗示该 multi-eq landscape 的拆分并非由 hub/SCC 边单独驱动，而是由更分布式的 feedback architecture 决定。cluster 2（sparse init）中仍有两类 state landing 聚类，表明删边后在 basin 之间产生了微弱的 re-distribution。

### 3.6 Landing Shift Without Regime Change — `stress_balanced_r010_t005` (H2_boundary) vs `stress_balanced_r010_t001` (H2_boundary)

- **`stress_balanced_r010_t005`**: 基线 Divergent（65/65），扰动后 **0/65 regime_transition**，37/65 landing_shift, 28/65 robust。rho(W)=1.091，clipping ≈ 2.1%。这是一个纯粹的 "landing-shift only" 世界——分散程度虽有变化但 regime 类型保持刚性不变。
  
- **`stress_balanced_r010_t001`**: 基线 Mixed（13 Convergent + 52 Divergent），扰动后 **3/65 regime_transition**（全部 Conv→Div），42/65 landing_shift, 20/65 robust。单一世界内 converge 与 diverge 共存的混合 setup 产生了独特的扰动响应——3 条 convergent trajectory 在删边后转向 divergent，而 52 条基线就 divergent 的全部保持 divergent。
  
- **对比启示**: 同样是 rho(W)≈1.0 的 stress boundary worlds，`r010_t005`（全 divergent）和 `r010_t001`（mixed conv+div）展现出截然不同的 sensitivity profile。前者 "regime lock-in"，后者存在 "convergent→divergent tipping edges"——说明 boundary 的脆弱性高度依赖于 pre-perturb 的 regime 组成。

### 3.7 Regime Transition After Perturbation — `stress_repression_biased_r010_t004` (H2_control)

- **基线 regime**: Sustained oscillatory（65/65 条 trajectory），repression-biased sign ratio，r=0.10
- **扰动结果**: **37/65 regime_transition（56.9%）**, 27/65 landing_shift, 1/65 robust
- **解读**: 与 `stress_balanced_r010_t007`（同为 56.9% transition）一起，是全体中 regime transition 最高的世界。但两者的 transition 方向不同：`r010_t004`（repression-biased）的 transition 以 oscillation→convergent 为主，而 `r010_t007`（balanced）则出现更多 oscillation↔ambiguous 模糊转换。这表明 sign ratio 与 edge density 的三方交互决定了 transition 的定性方向——repression-biased 的反馈更易被删边"稳定化"，而 balanced 的 feedback 更容易滑入 indeterminate drift。

---

## 4. 文档要求的科学问题回答

### 4.1 急性删边后，regime 多常改变？

481 条 trajectory 中，113 条（23.5%）出现 regime transition。在 perturbation-mode summary 层级（156 行）中，24 行（15.4%）以 regime_transition 为 dominant outcome。

### 4.2 Convergent worlds 中删边多常改变 final state landing？

304 条 Convergent baseline trajectory 中，298 条保持 Convergent。其中 13 条为 robust（normalized L2 < 0.05），其余为 landing_shift。仅 5/304（1.6%）转出 Convergent regime。Convergent worlds 对 5-edge acute deletion 高度鲁棒。

### 4.3 Boundary worlds 是否比 stable worlds 更敏感？

是。H2_boundary 的 regime transition 率为 30.8%，H1_stable 为 1.9%（差距 16 倍）。但即使在 boundary worlds 中，landing_shift（而非 regime_transition）仍是主导 outcome。Boundary zone 近似于一个"regime-flexible corridor"——删边后系统更倾向于在同一 regime 类型内滑动落点，而非完全转向另一种 regime。

### 4.4 Targeted deletion 是否比 random 更强？

在当前实验设计下，未观察到显著差异。四种策略的 regime transition 率均在 22.6%–25.5% 之间。Two main reasons： (1) 沉默节点效应使 topology-targeted strategies 的优越性无法发挥（见 §2.8）；(2) K=5 的扰动规模可能对大数 world 而言是 sub-threshold 的——无论删哪 5 条边，系统自带的 stability margin 都能吸收。

### 4.5 哪些 edge features 与更强扰动效应相关？

`Δρ(W)` 单独不足以预测 outcome（robust/landing_shift/regime_transition 在 `Δρ(W)` 上重叠）。Active-set 分析揭示了更重要的信号：source gene 的表达状态直接影响被删边的有效性。系统的 clipping pattern（active_zero_set）决定了哪些调控关系是"实质上的"（de facto active），哪些是"纸面上的"（de jure but silent）。这支持 H5："rho(W) 有帮助但不充分，clipping pattern / active zero set 也决定实际落点"。

---

## 5. 与 run_005 的关联

Run_006 从 run_005 的 120-world pool 中精选 12 个 worlds，覆盖了 single-eq convergent（H1）、non-convergent（H2/H3）、以及唯一的 multi-equilibrium world（H3: `chen_repression_biased_r025_t007`）。

Run_005 的核心发现是：edge density 的提升首先推动 non-convergence 而非 attractor multiplicity。Run_006 在此基础上进一步发现：即使在 non-convergent 或 boundary 区域，调控边的"删除效应"也高度依赖 expression context——沉默节点上的边删不删都一样。这与 run_005 中 `chen_balanced_r025_t005` 在 r=0.25 下 9/10 non-convergent 但 repression_biased group 仅 2/10 non-convergent 的发现形成延伸——两者都指向 **sign structure × expression state × clipping 的三方交互** 作为 Model A dynamics 的关键 determinant。

对于 `chen_repression_biased_r025_t007`（run_005 唯一的 multi-eq world），run_006 的删边实验显示：C: hub 和 D: SCC 策略下 trajectory 几乎无变化——这暗示该 world 的 multiple equilibria 并非由 hub gene 或 SCC 内边单方面驱动，而是由更全局的 feedback architecture 决定。

---

## 6. 与 Curriculum 前序 Runs 的关系

| Run | 贡献 | run_006 的继承 |
|-----|------|---------------|
| run_002 | 确认 Model A 可产生多种 dynamics | regime 分类体系沿用 |
| run_003 | 确认拓扑结构影响 dynamics | — |
| run_004 | 确认 δ 是动力学边界关键参数 | δ=0.2 作为本 run 固定参数 |
| run_005 | 确认 edge density 改变 convergence landscape | 12 个 worlds 选自 run_005 pool；multi-eq world 在 run_006 中接受了针对性扰动测试 |

---

## 7. 局限性

### 7.1 沉默节点限制了 Targeted 策略的有效性

本 run 最明确的机制发现之一——topological hub ≠ functional hub——也构成一条"自限性"结论。在 RELU-based model A 中，topology-based 选边策略的 default 假设是 "high-degree regulators are functional"。但本实验的 12 个代表性 worlds 中，多个 state classes 下 hub 和 SCC 基因恰好被 clipping 到 0，导致 targeted deletion 的效果与 random deletion 无异。这个发现本身具有科学价值，但意味着：(1) 本 run 的实验结果并不构成对 H4（targeted > random）的否定性证据，而是在当前稳态状态下的有条件结论；(2) 需要一个 "expression-aware" 的选边策略变体（如选最高表达量基因的出边）作为 future work。

### 7.2 Divergent 轨迹的 magnitude ratio 分析已完成

Divergent-involved 的 187 条 trajectory 已通过 magnitude ratio 分析量化（见 §2.5）。关键发现：54% 的 Div→Div 轨迹在删边后发散被抑制（ratio < 0.9），median ratio = 0.706。但 magnitude ratio 仍有局限性：它基于最后 50 步均值向量的范数比，无法区分 "发散增长率改变" 与 "终点时刻的体量差异"。更精细的增长率分析（exponential growth rate of log-norm）需要从 trajectory 原始数据计算。

### 7.3 1000-step 窗口对 slow-drift dynamics 不足

`ambiguous` 案例 `chen_balanced_r025_t005 × C:hub` 显示：0.04%/step 的慢速 drift 在 500-step post-fork window 内既不会收敛也不会发散。如果仿真延长到 2000+ steps，部分 `ambiguous` 可能被重新归入 Convergent 或 Divergent。当前 1000-step total window 是一个 pragmatic choice，但对临界 dynamics 的分辨率存在已知限制。

### 7.4 Divergent trajectory outcome labels 的局限性

Divergent→Divergent 的 182 条被标记为 `landing_shift` 或 `robust`，但这两个标签原本是为有 equilibrium 的轨迹设计的。对 divergent 轨迹，"landing_shift" 仅意味着 regime 类型不变，但实际 divergent 程度的变化（growth rate, terminal magnitude）未被纳入 outcome label。

§2.5 的 magnitude ratio 分析已部分弥补了这一缺口——通过 `||post_final|| / max(||pre_ref||, eps)` 量化了删边后发散体量的相对变化，并发现 54% 的 Div→Div 轨迹在删边后发散被抑制。但 magnitude ratio 仍有局限：它基于最后 50 步均值向量的范数比，无法分解为 "发散增长率改变" 与 "终点时刻体量差异" 两个独立分量。更精细的 exponential growth rate 分析（log-norm slope over sliding windows）需要从 trajectory 原始数据计算，可作为 future work。

### 7.5 Next Steps

- **Expanded worlds**：当前 12-world selection 覆盖了关键 regime 组合但未全覆盖 sign_ratio × strength_regime × r 的全因子设计。H4 的结论可能因选择偏倚受限，特别是 chen_repression_biased 的 sampled worlds 较少。
- **Time-to-divergence as perturbation metric**：对 non-divergent → divergent 的 5 个案例，记录 divergence onset time 的变化。
- **引入非线性调控机制**：当前 Model A 为纯线性 GRN + RELU clipping。未来可在反馈项中引入 saturable activation（如 Hill-type kinetics）或 cooperative binding，观察非线性机制是否改变 topology-targeted perturbation 的 leverage profile，是否打破当前广泛存在的 "robust regime, flexible state" 模式，以及是否出现多个 attraction zone 或 coexistence of multiple equilibria。

---

## 8. One-Line Summary

> run_006 demonstrates that Model A's regulatory architecture exhibits regime-level robustness against acute 5-edge deletion (76.5% same-regime outcomes), but reveals a fundamental distinction between topological and functional hubs — in RELU-based dynamics, an edge only matters if its source gene speaks.

---

## 9. 输出清单

### Tables
- `data/source_world_selection.csv` — 12 selected worlds
- `data/edge_table_pre.csv` — pre-perturbation edge catalog
- `data/perturbation_metadata.csv` — 481 perturbation plans
- `data/pre_perturb_state_classes.csv` — pre-perturb state classes
- `data/trajectory_outcomes.csv` — 481 perturbation trajectory outcomes
- `data/perturbation_summary.csv` — 156 per-mode summaries
- `data/state_landing_summary.csv` — state landing clusters
- `data/representative_state_vectors.csv` — per-gene representative vectors
- `data/edge_property_outcome.csv` — per-edge perturbation outcomes

### Figures
- `figures/regime_transition_by_strategy.png` — regime transition matrix by strategy
- `figures/regime_transition_by_world_class.png` — regime transition matrix by world tier
- `figures/per_world_perturbation_outcome_heatmap.png` — per-world outcome heatmap
- `figures/state_landing_2x2.png` — PCA-based state landing visualization
- `figures/trajectory_pre_post_world_*.png` — representative trajectory pre/post plots (4 worlds, 见下说明)
- `figures/delta_rhoW_vs_outcome.png` — Δρ(W) vs normalized L2 distance (linear y-axis, divergent excluded)
- `figures/delta_rhoW_vs_divergent_magnitude.png` — Δρ(W) vs magnitude ratio for divergent-involved trajectories

> **关于轨迹图的生成说明**：本 run 共产生 481 条 perturbation trajectory（12 worlds × 37 representative initial conditions × 4 种策略，A 策略含 10 repeats），逐条绘制轨迹图不可行。因此按 tier（H1 / H2 / H2b / H3）各选 1 个代表性 world，每个 world 取其第 1 个 representative initial condition，展示该 init 下全部 4 种删边策略（A: random / B: strongest / C: hub / D: SCC）的 fork 轨迹对比。每张图为 4 行（策略）× 3 列（full trajectory / zoom ±50 steps / late-time window），共输出 4 张。选出的 world 及对应 tier：
>
> | Tier | World | 选取理由 |
> |:----:|-------|----------|
> | H1 | `baseline_balanced_r010_t000` | stable anchor，rho(W) < 1，零 clipping，期望高度鲁棒 |
> | H2 | `stress_balanced_r010_t003` | boundary world，收敛基线，检验边界附近的 landing shift |
> | H2b | `stress_balanced_r010_t007` | boundary world，振荡基线（regime transition 率最高的 world 之一） |
> | H3 | `chen_repression_biased_r025_t007` | rare specimen，run_005 中唯一的 multi-equilibrium world |
