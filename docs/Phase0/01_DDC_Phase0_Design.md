# Designed Digital-Cell (DDC) Model

> Current Version: v1.0
> Last Updated: 02/16/2026
> Status: Frozen

## Change Log

v1.0 – Initial Phase 0 Freeze

---

## 1️. 变量体系（统一命名）

- $X_i(t)$ — mRNA
    
- $P_i(t)$ — functional protein
    
- $Z_i(t) \in [0,1]$ — chromatin gate
    
- $N(t)$ — population-level state (observation only)
    

**Phase0 不再使用 U / Y。**

**Rationale：**

- Phase 0 的核心目标是跑通“扰动→表达”的闭环；U/Y 的区分会引入额外隐变量与参数，

  增加不可辨识与实现复杂度。
    
- 用单一 protein 层 $P$ 继承原模型中 functional protein 的角色（≈Y），

  保持因果语义：$(P,Z)\to X$。
    
- 预留 Phase 1 恢复 U/Y 分层的接口（通过在 protein 模块内部扩展，不改主循环）。

---

## 2️. Regulatory Module

冻结为：

$$
\tilde P_j(t)
=
\frac{P_j(t)}
{\sum_{k=1}^{G} P_k(t) + \epsilon}
$$

$$
d_i = |P(i)|
$$
$$
\text{TFinput}_i(t)
=
\left[
\prod_{j \in P(i)}
\left(
\tilde P_j(t)
\right)^{a_{ij}}
\right]^{\frac{1}{d_i}}
$$



- 使用 multiplicative form
    
- $\epsilon = 10^{-8}$
    
- $a_{ij} > 0$

> $P(i)$ is defined by the transcription graph specified in the System Definition.

**Rationale：**

- $~\tilde P$ 把“绝对量”变成“组成比例”，与资源守恒/竞争机制一致，并显著降低尺度漂移导致的不稳定。
    
- 乘法形式天然体现 combinatorial requirement（AND-like），

  能在 Phase 0 就引入强非线性与高阶交互，而不需要显式中间层。
    
- 相比 min：可微、数值更平滑；相比线性：组合性更强，更符合项目的“机制挑战性”。

---

## 3️. Chromatin Gate

冻结为：

$$
Z_i(t+1)
=
\sigma\!\left(
\alpha_i
+
\sum_{j \in E(i)} \beta_{ij}\tilde P_j(t)
\right)
$$

$$
\sigma(x) = \frac{1}{1 + e^{-x}}
$$


- 使用连续 gate
    
- $δ_z = 1$（无记忆）
    
- 不启用 Markov gate
    
**Rationale：**

- gate $Z\in[0,1]$ 是 Phase 0 引入“潜在不可辨识性”的关键，同时仍保持计算稳定。
    
- $\delta_z=1$（无记忆）是为了冻结最小动态结构：避免时间尺度耦合与滞后参数爆炸。
    
- Markov ON/OFF 会引入离散随机过程与 burst 统计，属于 Phase 1+ 复杂度轴，

  不在 Phase 0 打开。

---

## 4️. Transcription

冻结为 Hill：

$$  
X_i(t+1)
=
(1 - \delta_{x,i}) X_i(t)  
+  
Z_i(t)  
\cdot  
\rho_i  
\cdot  
\frac{  
\text{TFinput}_i(t)^{n_i}  
}{  
K_i^{n_i}  
+  
\text{TFinput}_i(t)^{n_i}  
}  
$$

- $0 < \delta_{x,i} \le 1$
    
**Rationale：**

- Hill 函数有上界，能避免无界增长，是离散时间系统最重要的稳定性保障之一。
    
- $\delta_{x,i}>0$ 强制存在稳态/吸引域，便于 run_001 sanity checks（尤其是 perturbation 后回归平衡）。
    
- 这一步把“非线性响应”集中在一个地方，方便后续替换为更复杂的 transcription noise/burst。

---

## 5️. Translation

最小形式：

$$  
P_i^{\text{raw}}(t+1)
=
(1 - \delta_{p,i}) P_i(t)  
+  
\gamma_i X_i(t)  
$$

**Rationale：**

- Phase 0 的 protein 层主要是把 X 映射为能进入 regulation 的实体；

  线性形式最小参数、最易调试。
    
- 引入 $P^{raw}$ 是为了把“生成”和“资源约束”解耦：先生成，再投影守恒。
    
- 未来要加入 RBP/修饰 等复杂机制，只需把 $g(\cdot)$ 换掉，projection 与下游不变。

---

## 6️. Resource Projection


$$
P_i(t+1)
=
\begin{cases}
P_i^{\text{raw}}(t+1)
\cdot
\frac{R_{\text{total}}}
{\sum_k P_k^{\text{raw}}(t+1)}
& \text{if } \sum_k P_k^{\text{raw}}(t+1) > R_{\text{total}}
\\[8pt]
P_i^{\text{raw}}(t+1)
& \text{otherwise}
\end{cases}
$$


**Rationale：**

- 这是 Phase 0 的“全局耦合机制”，提供资源竞争与间接相互作用，不需要显式建模“其他基因”。
    
- 线性缩放投影是 deterministic、可复现、可解释的守恒实现方式，便于测试与审计。
    
- 更复杂的资源分配（soft constraint/目标函数投影）留到 Phase 1。

---

## 7️. Fate

$$
N(t+1)
=
N(t)
+
r N(t)
\left(
1 - \frac{N(t)}{K_{\text{pop}}}
\right)
$$

- 不进入 X/P 动力学
    
- 仅影响 observation

**Rationale：**

- Phase 0 的目标是 GRN world-building；

  Agent-based birth/death 会引入分支过程与状态集合管理，复杂度指数上升。
    
- 用 $N(t)$ 影响 observation 可以保留“population-level effect”的语义，

  同时不破坏 GRN 动力学闭环。
    
- Phase 1+ 可扩展为：A) fate→GRN gate（弱耦合）或 B) birth–death branching（强耦合），

  接口提前冻结避免未来重构。

