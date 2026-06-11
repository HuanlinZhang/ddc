# Run 002 简要总结

**实验**: Classical Linear GRN 在 Chen-style Sparse Topology 下的动力学刻画
**日期**: 2026-05
**模型**: Model A（单层 abundance-compatible piecewise-linear GRN）

---

## 1. 动机

在研究 DDC 分层带来的 representation 变化之前，先建立一个干净的 baseline：**最简 GRN 在受控参数变化下会涌现出哪些动力学行为？** Run 002 就是这张参考坐标系。

---

## 2. 模型

单层离散时间 GRN，仅有 non-negativity clipping 作为非线性：

$$
x_i(t+1) = \max\bigl(0,\; (1-\delta_i)x_i(t) + b_i + \sum_{j} a_{ij}x_j(t)\bigr)
$$

- **G = 20** 个基因，**δ = 0.2**（self-retention 0.8），**b = 0.1**（basal transcription）
- **拓扑**：Chen-style sparse directed graph，edge density r = 0.1（20 基因中随机抽取 38 条有向边）
- **边强度**：|a_ij| ~ Uniform(a_min, a_max)，sign 独立分配为 activation (+1) 或 repression (-1)
- **初始表达量**：每个 cell 的 20 个基因表达量独立从 **Uniform(0, 1)** 中采样（每个 cell 用独立 seed 以保证可复现）

核心动力学是纯线性的，唯一的非线性来自 non-negativity clipping（abundance 不能为负）。

---

## 3. 实验设计

参数扫描沿**两个维度**展开：

| 维度 | 水平 | 设计意图 |
|------|------|----------|
| **interaction strength**（4 档） | baseline [0.02, 0.1] → chen_moderate [0.1, 0.2] → stress [0.1, 0.3] → chen_stress [0.2, 0.4] | interaction 从"远小于 decay"逐步增强到"超过 decay 2×"，观察动力学如何退化 |
| **sign ratio**（2 档） | balanced（activation:repression = 1:1）、repression_biased（1:2） | 增加 repression 比例是否能系统性提升稳定性？ |

- **8 个 combo × 10 个独立 topology = 80 个 world**，每个 world 5 个随机初始 cell → **400 条 trajectory**
- 每条 trajectory 运行 T = 1000 步
- 每个 world 归入以下 7 类 attractor 之一：

| Type | 名称 | 简要说明 |
|------|------|----------|
| A | stable equilibrium | 快速收敛到定点 |
| B | slow convergence | 收敛但超过 200 步 |
| C | damped oscillation | 振荡幅度逐步衰减 |
| D | sustained oscillation | 振荡幅度不衰减 |
| E | runaway divergence | 表达量无限增长 |
| F | numerical collapse | 所有表达量趋近 0 |
| G | unclassified | 防御性兜底 |

- 选取 9 个 world 做扰动分析：t=500 时对表达量最高的基因执行瞬时 knockdown（一个 timestep 内设为 0）

---

## 4. 核心结果

### 4.1 interaction strength 是 regime 切换的主控轴

| strength regime | Type A（稳定） | 发散（Type E） |
|----------------|:---:|:---:|
| baseline | **100%** | 0% |
| chen_moderate | 80% | 0% |
| stress | 50% | 20% |
| chen_stress | 55% | 35% |

随着 a_max 从 0.1 增长到 0.4，系统从全稳定退化到 balanced 条件下 60%+ 非稳定。仅靠 interaction strength 一个变量的变化，就能使系统依次经历全 Type A → 混合 B/C/D → 发散主导的全过程。

[图: figures/topo_strength_interaction.png]

### 4.2 repression 系统性增强稳定性

在同等 interaction strength 下，repression_biased（1:2）网络一致更稳定：

| 指标 | balanced（40 worlds） | repression_biased（40 worlds） |
|------|:---:|:---:|
| Type A 占比 | 55% | **87.5%** |
| 发散率 | 20% | **7.5%** |
| 中位 ρ(A) | 1.025 | 0.991 |

在高 stress regime 下差距最大：chen_stress 中 balanced 仅 30% Type A，而 repression_biased 仍保持 **80%** Type A。解释：负反馈充当了"动态阻尼器"——负调控回路将表达约束在 bounded region，防止 positive feedback runaway。

[图: figures/repression_ratio_vs_stability.png]

### 4.3 Linear GRN 的动力学比预期丰富

尽管没有非线性激活函数、没有 protein layer、没有延迟，系统仍然涌现出：

- **慢收敛**（Type B）：80 个 world 中 9 例
- **衰减振荡**（Type C）：2 例
- **持续振荡**（Type D）：1 例

这些行为完全来自 signed 调控回路 + non-negativity clipping 的组合，不需要 Hill-function 非线性。

[图: figures/canonical_trajectories.png, figures/convergence_time.png]

### 4.4 谱半径 ρ(A) 仅在不触发 clipping 时有效

对于离散时间线性系统 x(t+1) = A x(t)，稳定性的经典判据是转移矩阵 A 的谱半径 ρ(A) < 1（所有特征值落在单位圆内）。这对应 Chen (2019) 在连续时间下考察 Jacobian 特征值的做法——两者都试图从线性谱中提取稳定性信息。本文检验这一判据在 piecewise-linear（含 clipping）条件下是否仍然成立。

在 **non-clipping 区间**该判据完美成立：20 个非 clipping world 全部 ρ < 1 且全为 Type A，无一反例。但一旦进入 **clipping-dominated regime**（interaction 足够强，基因频繁被 non-negativity clipping 钳制到 0），ρ(A) 的预测力系统性失效：

- 37 个 clipping Type A 中 **54%（20 个）ρ > 1**——按线性理论应发散，实际却稳定
- 反过来，1 个发散 world（ρ = 0.995）ρ < 1 却发散

原因：non-negativity clipping 引入了非线性边界条件。ρ(A) 描述的是 unbounded linear core 的扩张/收缩倾向，一旦 clipping 频繁触发，系统实质上是 piecewise-linear，纯线性谱分析不再充分。

[图: figures/spectral_radius.png, figures/clipping_dominated_by_regime.png]

### 4.5 扰动恢复：收敛系统具有弹性

A/B 类型 world（已收敛者）均在瞬时 knockdown 后 ~40-60 步恢复，attractor 类别未变。振荡型和发散型 world 扰动后保持原有动力学行为。**没有任何 world 因 knockdown 改变了 attractor category**。

[图: figures/perturbation/perturbation_type_[a-e].png]

---

## 5. 核心结论

1. **interaction strength 是 regime 切换的主控轴**——固定其他参数，仅靠增大 interaction 强度就能驱动系统从全 Type A 逐步过渡到混合 B/C/D，再到发散为主
2. **repression 在 linear GRN 中充当 universal stabilizer**——同等强度下 repression_biased 发散风险约为 balanced 的 1/3（7.5% vs 20%）
3. **最简 linear GRN 能涌现多种动力学行为**——慢收敛、衰减振荡和持续振荡均可出现，这些行为来自 signed 调控回路 + non-negativity clipping，不需要 Hill-function 非线性
4. **线性谱理论在 clipping 活跃时不足**——ρ(A) < 1 仅在 non-clipping 区间有效，进入 clipping-dominated regime 后与 attractor 几乎脱钩，这是后续所有分析的 methodology caveat
