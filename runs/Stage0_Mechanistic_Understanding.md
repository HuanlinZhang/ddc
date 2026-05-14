# Stage 0 Mechanistic Understanding

> 基于 Phase 0 全部实验结果 (run_004 / run_005 / run_006 / run_007 / run_008) 的机制分析

---

## 核心因果链

```
模型结构（乘法TF调控 + Hill非线性 + 资源约束）
        ↓
动力学机制（TF衰减被放大 → 正反馈失活 → 稳态无法维持）
        ↓
系统行为（steady/collapse取决于初始TF激活是否足够强）
        ↓
扰动响应（瞬时扰动可恢复，持久扰动导致结构改变不可逆）
        ↓
剂量敏感性（系统在最大容量运行，无安全缓冲）
```

**一句话总结**: DDC系统的稳定性依赖于乘法TF调控产生的正反馈激活，而其脆弱性正源于此——任何持久性扰动都会破坏这个精心平衡的激活网络。

---

## 1️⃣ Viability & Regime Structure

### 核心问题：为什么有些world达到steady而有些collapse？

#### 因果链分析

```
初始TF激活（seed）
        ↓
乘法调控放大：TFinput = (∏P̃ⱼ^aᵢⱼ)^(1/dᵢ)
        ↓
Hill函数激活：hill = TFinputⁿ / (Kⁿ + TFinputⁿ)
        ↓
mRNA合成：X' = (1-δₓ)X + Z·ρ·hill
        ↓
蛋白合成：P' = (1-δₚ)P + γ·X
        ↓
循环往复（正反馈）
```

#### 机制解释

**为什么collapse发生？**

当TF蛋白略低于临界值时，乘法调控结构将其变化**放大**。例如：
- 假设TF下调10%，则TFinput ≈ (0.9)^n ≈ 0.69（假设n=4）
- 这导致Hill函数激活大幅下降
- 激活下降减少mRNA合成
- mRNA减少导致蛋白进一步降低
- 形成**正反馈循环**：低TF → 更低TF → 彻底失活

**为什么steady worlds是少数？**

需要同时满足：
1. 足够强的初始TF激活（跨过Hill函数阈值）
2. 平衡的TF网络拓扑
3. 足够的正反馈强度
4. 正确的参数组合（α, β, ρ, δ等）

Monte Carlo采样中，只有极少数随机参数组合能满足这些严苛条件。

**为什么增加TF support不能rescue collapse？**

Collapse是**self-reinforcing**的。关键原因：
- 当Hill激活降到最低点附近，即使大幅增加TF，Hill函数变化也很小
- 系统已经进入"TF耗竭"的吸引子
- 需要的是改变系统参数（结构），而非外部注入

**为什么初始条件变化不改变最终命运？**

系统存在单一的attractor landscape：
- **参数决定basin**：给定参数集合，系统最终趋向唯一的attractor（steady或collapse）
- **初始条件只影响速度**：在正确basin中，会更快收敛；在错误basin中，会更快崩溃

---

## 2️⃣ Transient Perturbation Response

### 核心问题：为什么系统在KD/OE后可以恢复？

#### 机制解释

**Attractor / Basin of Attraction 结构**

```
        ┌─────────────────────────────────────┐
        │         Steady Basin                │
        │                                     │
        │    ┌──────────────────────┐        │
        │    │   Steady Attractor   │        │
        │    │      (P* ≈ 0.02)      │        │
        │    └──────────────────────┘        │
        │                                     │
        │         ↑ 恢复                      │
        │   (瞬时KD后自然趋向)                │
        └─────────────────────────────────────┘

        ┌─────────────────────────────────────┐
        │        Collapse Basin                │
        │                                     │
        │    ┌──────────────────────┐        │
        │    │  Collapse Attractor   │        │
        │    │      (P* ≈ 0)         │        │
        │    └──────────────────────┘        │
        │                                     │
        └─────────────────────────────────────┘
```

**为什么steady systems的amplification更高？**

- Steady系统：TF调控网络处于活跃状态，能对扰动产生显著响应
- Collapse系统：网络已经"死寂"，响应微弱

**系统的robustness来自哪里？**

- **Transient robustness**：来自steady basin的结构——只要不持久改变系统参数，系统会自然恢复
- 关键点：KD/OE是**do-operator**（状态修改），不是**mechanism modification**（参数修改）

---

## 3️⃣ Structural Perturbation

### 核心问题：为什么KO TF会导致collapse？

#### 机制解释

**KO vs KD 的本质区别**

| 操作 | 类型 | 效果 |
|------|------|------|
| KD (Knockdown) | State modification | 瞬时状态改变，dynamics不变 |
| **KO (Knockout)** | Parameter modification | 永久性参数改变，dynamics被修改 |

**KO = ρᵢ = 0**：直接删除该基因的mRNA合成能力

```
乘法结构：TFinput = (∏P̃ⱼ^aᵢⱼ)^(1/dᵢ)
                ↓
        如果某个TF被KO → 缺少关键因子
                ↓
        TFinput显著下降
                ↓
        乘法放大效应
                ↓
        Hill函数激活降低 → 系统崩溃
```

**为什么基因0-3是100%崩溃，而基因4-5是75%？**

- **基因0-3**：位于网络核心，高度互联，KO后无法被其他TF补偿
- **基因4-5**：位于网络边缘，部分冗余，可能被其他TF替代

**关键洞察**：TF之间存在**层级依赖**——有些TF是"关键节点"，移除后整个网络失效。

---

## 4️⃣ Dosage Sensitivity

### 核心问题：为什么0.5× attenuation就已经显著影响系统？

#### 机制解释

**乘法结构的剂量放大效应**

```
λ = 0.5（50%剂量减少）
        ↓
蛋白合成：P' = (1-δₚ)·P + λ·γ·X
        ↓
稳态值：P* ≈ λ·γ/δₚ · X
        ↓
λ=0.5 → P* ≈ 0.5·P*_baseline
        ↓
TFinput = (∏P̃ⱼ^aᵢⱼ)^(1/dᵢ) ≈ (0.5)^n · TFinput_baseline
        ↓
对于n=2: TFinput ≈ 25% 原始值
对于n=4: TFinput ≈ 6.25% 原始值
        ↓
Hill函数激活急剧下降
```

**为什么collapse probability随dose单调上升？**

- 系统运行在**stability boundary边缘**
- λ越小，系统越远离稳定区
- 非线性响应：剂量-崩溃曲线不是线性，而是加速增长

**为什么不同TF的sensitivity曲线不同？**

- 不同TF在网络中的**拓扑位置**不同
- **Degree（连接数）**：连接越多，乘法结构中权重越大
- **位置（hub vs periphery）**：hub TF移除影响更大

**什么是haploinsufficiency-like行为？**

定义：50%基因剂量减少导致系统失败

在本模型中：
- TF是系统的**critical load-bearing component**
- 系统在最大容量运行，没有剂量缓冲（buffer = 0）
- 任何持久性减少都是"最后一根稻草"

---

## 5️⃣ 最终总结

### 这个系统的稳定性是如何产生的？

**稳定性的来源**：

1. **乘法TF调控产生的正反馈**：TF激活蛋白合成，蛋白激活TF表达，形成自维持循环
2. **Hill函数的阈值特性**：一旦激活超过阈值，系统能自我维持
3. **Attractor landscape**：steady state是一个robust attractor

**稳定性是如何被破坏的？**

1. **乘法结构的双刃剑**：放大激活的同时，也放大衰减
2. **Self-reinforcing collapse**：低TF → 更低TF，形成不可逆的下螺旋
3. **Attractor切换**：持久扰动改变dynamics，使系统进入collapse basin
4. **无剂量缓冲**：系统运行在临界点，任何减少都是致命的

**核心洞察**：

> 系统同时具有robustness（attractor保证恢复）和fragility（乘法结构放大扰动）。Transient perturbation可恢复，persistent structural perturbation不可逆。

---

## 附录：关键模型方程

### TF输入计算（乘法调控）

```python
TFinputᵢ = (∏ⱼ P̃ⱼ^aᵢⱼ)^(1/dᵢ)
```

### Hill激活函数

```python
hillᵢ = TFinputᵢⁿ / (Kᵢⁿ + TFinputᵢⁿ)
```

### mRNA更新

```python
X'ᵢ = (1-δₓ)·Xᵢ + Zᵢ·ρᵢ·hillᵢ
```

### 蛋白更新（含持久衰减）

```python
P'ᵢ = (1-δₚ)·Pᵢ + γᵢ·Xᵢ  # baseline
P'ᵢ = (1-δₚ)·Pᵢ + λ·γᵢ·Xᵢ  # run_008: 持久衰减模型
```
