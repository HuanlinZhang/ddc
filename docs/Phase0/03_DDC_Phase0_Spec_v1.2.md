# DDC Phase 0 – 执行规范（Execution Specification）

> Version: v1.2（统一 Perturbation + Intervention）  
> Status: Frozen  
> Last Updated: 2026-04-08

---

# 0. 目的（Purpose）

本文件定义 Phase 0 的 **执行契约（execution contract）**。

冻结内容：

- update order（更新顺序）
    
- module interfaces（模块接口）
    
- randomness rules（随机性规则）
    
- perturbation semantics（扰动语义）
    
- intervention semantics（干预语义）
    
- sanity checks（正确性检查）
    

本版本整合：

```text
v1.0（基础动力学）
+ v1.1（intervention 机制）
+ KO 语义修正
```

---

# 🔥 0.1 核心原则（Core Principle）

```text
Dynamics is frozen.
```

```text
动力学完全冻结，
只允许 perturbation 和 intervention 作用于系统。
```

---

# 1. 对象模型（Object Model）

---

## 1.1 World（不可变）

包含：

- regulatory graphs
    
- 参数（ρ, a_ij, α_i 等）
    
- 全局常数
    

```text
World MUST NOT be mutated during simulation
```

---

## 1.2 State

```text
S(t) = {X(t), P(t), Z(t), N(t)}
```

---

# 2. 更新顺序（Update Order）

```text
1. tilde_P
2. TFinput
3. Z(t+1)
4. X(t+1)（使用 Z(t)）
5. P_raw
6. resource projection
7. N(t+1)
```

```text
❌ 不允许修改
❌ intervention 之外不得插入新步骤
```

---

# 3. 模块接口（Module Interfaces）

完全继承 v1.0：

- normalize_protein
    
- compute_TFinput
    
- update_chromatin
    
- update_mRNA
    
- update_protein_raw
    
- apply_resource_projection
    
- update_fate
    

```text
❌ 不允许修改接口或数学形式
```

---

# 4. 扰动体系（Perturbation vs Intervention）

---

## 🔴 Phase 0 中存在两种完全不同的操作：

---

# 4.1 Type A — 参数扰动（Parameter Perturbation）

## 定义

```text
perturbation modifies World（参数层）
→ 改变 system dynamics
```

---

## 允许操作

- ρ_i = 0（基因 KO）
    
- 修改 a_ij
    
- 修改 α_i
    
- 修改 resource 参数
    

---

## 性质

```text
✔ 持续生效（persistent）
✔ 改变系统动力学
✔ 属于 system-internal
```

---

## 接口

```python
def apply_perturbation(world, state, config) -> (world_perturbed, state)
```

---

## 约束

```text
✔ 必须复制 world（immutable）
✔ 原始 world 不得修改
✔ 后续 dynamics 使用新参数
```

---

## 🔬 生物学解释

```text
ρ_i = 0 → transcription 被关闭
→ gene knockout（KO）
```

---

# 4.2 Type B — 状态干预（State Intervention）

---

## 定义

```text
intervention modifies State（状态层）
→ 不改变 dynamics
```

---

## 允许操作

- X_i = 0
    
- P_i scaling
    
- 多基因 override
    

---

## 性质

```text
✔ 单时间点（single-step）
✔ 不持续
✔ external do-operator
```

---

## 接口

```python
def apply_intervention(state, config) -> state
```

---

## 核心约束

```text
✔ 仅在指定时间点执行

✔ intervention step：
    ❌ 不执行 resource projection

✔ 不得修改 World

✔ intervention 之后：
    系统恢复正常 dynamics
```

---

## 🔬 生物学解释

```text
X_i(t) = 0 → 瞬时 knockdown
```

---

# 🔑 4.3 术语标准化（Terminology）

```text
KO（knockout）：
    parameter-level
    → ρ_i = 0

Knockdown：
    state-level
    → X_i(t) = 0 或 P_i(t) = 0
```

---

# 🚨 4.4 关键区分（必须理解）

```text
perturbation ≠ intervention
```

|类型|修改对象|是否持续|是否改变 dynamics|
|---|---|---|---|
|perturbation|World|✔|✔|
|intervention|State|❌|❌|

---

# 5. 主循环（Main Loop）

```python
for t in range(T):

    # Intervention Hook
    if t == intervention_time:
        state = apply_intervention(state, config)
        X, P, Z, N = state

    tilde_P = normalize_protein(P, world)
    TFinput = compute_TFinput(tilde_P, world)
    Z_next = update_chromatin(tilde_P, world)
    X_next = update_mRNA(X, Z, TFinput, world)
    P_raw = update_protein_raw(P, X, world)
    P_next = apply_resource_projection(P_raw, world)
    N_next = update_fate(N, world)

    X, P, Z, N = X_next, P_next, Z_next, N_next
```

---

# 6. Monte Carlo / Multi-cell

```text
与 v1.0 完全一致
```

---

# 7. Sanity Checks（新增）

---

## 7.1 Intervention Consistency

```text
✔ intervention 只执行一次

✔ intervention 后下一步：
    必须满足 resource constraint

✔ 无 intervention：
    行为必须与 v1.0 完全一致
```

---

# 8. 兼容性说明（Compatibility）

---

## v1.0 行为

```text
无 intervention → 完全等价 v1.0
```

---

## 🔧 KO 语义修正（关键）

```text
旧定义：
X_i(t) = 0 → KO

新定义：
X_i(t) = 0 → knockdown（intervention）
ρ_i = 0 → KO（perturbation）
```

---

# 🔒 最终状态（Final State）

```text
✔ dynamics 冻结
✔ world 生成冻结
✔ 扰动体系统一
✔ 接口语义清晰

Phase 0 execution contract 完全闭环
```