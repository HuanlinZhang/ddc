# Designed Digital-Cell (DDC) Model

## Architecture Master Plan (Phase 0 / MVP)

> Version: v1.2（统一 Perturbation + Intervention）
> Status: Frozen  
> Scope: Phase 0 World Construction
> Last_Updated: 2026-04-08

---

# 0. 项目定位

本文件定义 **DDC Phase 0 的系统架构蓝图**。

职责：

- 将 2–11 模型转化为可执行系统
    
- 统一 Design / System / Execution 三层语义
    
- 定义 Phase 0 的能力边界
    

目标：

```text
构建一个：
✔ 可运行（dynamical system）
✔ 可扰动（internal + external）
✔ 可观测（scRNA-like）
✔ 可复现（seed-deterministic）
的数字细胞世界
```

---

# 0.1 Phase 0 实现口径

核心策略：

```text
Monte Carlo world-building
+ discrete-time deterministic dynamics
```

---

## 0.1.1 单次 simulation run 结构

每个 run 分为两层：

### (A) World-level（run 内固定）

- Regulatory graphs
    
- Parameter sets
    

一次采样，整个 trajectory 固定

---

### (B) State dynamics（逐时间步）

状态变量：

```text
S(t) = {X(t), P(t), Z(t), N(t)}
```

每步更新（顺序冻结）：

```text
P → TFinput → Z → X → P_raw → projection → N
```

update order 不可改变

---

# 0.2 因果结构（Phase 0）

原始 2–11：

```text
X → U → Y → Z → X
```

Phase 0 简化：

```text
X → P → Z → X
```

核心闭环：

```
P(t),Z(t) → X(t+1)
X(t+1) → Praw(t+1) → projection → P(t+1)
```

---

# 🔥 0.3 扰动体系（关键升级）

Phase 0 明确支持两类操作：

---

## Type A：System-internal perturbation（run_005）

```text
✔ 修改参数 / 受约束 state
✔ 属于系统内部机制
✔ 必须 obey dynamics + resource constraint
```

---

## Type B：External intervention（run_006）

```text
✔ 在单一时间点直接修改 state（X / P）
✔ 属于外部 do-operator
✔ 可以暂时违反 resource constraint
✔ 下一步恢复正常 dynamics
```

---

## 🔒 核心区分

```text
run_005：机制内扰动（mechanism-level）
run_006：状态干预（state-level intervention）
```

---

## 🔒 原则（必须写清）

```text
intervention 改变 state，但不改变 dynamics
```

与 Execution Spec v1.1 完全对齐

---

# 1. MVP 四个不可删除机制

1. Combinatorial regulation
    
2. Strong nonlinearity
    
3. Global resource conservation
    
4. Fate / population effect
    

---

# 2. Phase 0 简化策略

---

## 2.1 Protein collapse（U/Y → P）

```text
U/Y 层合并为 P
P 直接承担 functional protein 角色
```

降维 + 保持调控语义

---

## 2.2 Chromatin

```text
连续 gate（sigmoid）
无时间记忆
```

---

## 2.3 Resource constraint

若 $ΣP>Rtotal：Pi ← Pi×(Rtotal/ΣP)$

deterministic global coupling

---

## 2.4 Fate

```text
N(t) logistic update
不进入 GRN 动力学
仅影响 observation
```

保留 population-level effect

---

# 3. 模块结构（Execution-ready）

---

## 3.1 Regulatory Module

- multiplicative combinatorial input
    
- 强非线性
    

---

## 3.2 Transcription

$$
X_i(t+1) = f(TFinput_i(t), Z_i(t))
$$
---

## 3.3 Translation + Projection

```text
X → P_raw → resource projection → P
```

---

## 3.4 Fate

```text
N(t+1) = logistic(N(t))
```

---

## 3.5 Perturbation / Intervention Interface

Phase 0 支持：

---

### Internal perturbation

```text
Parameter-level KO:
- ρ_i = 0（perturbation）
- 参数修改（ρ, a_ij, α_i）
- resource 参数变化
```

---

### External intervention

```text
State-level knockdown:
- X_i = 0 / P_i = 0（intervention）
- P_i ← scale × P_i
- X_i ← override
- 仅在指定时间点执行（intervention 必须是 single-step（t == T），禁止持续作用）
- 不执行 projection（当步）
```

---

## 🔒 接口原则

```text
✔ perturbation must NOT mutate original world (immutable)
✔ perturbation operates on a copied world (world_perturbed)
✔ intervention 不进入 update rule
✔ dynamics 始终由 Design 决定
```

---

# 3.6 Observation Model

默认：

$$
C_c(t) = X_c(t)
$$

输出包含：

- run_id
    
- seed
    
- perturbation / intervention
    
- gene_id
    
- time
    

---

# 4. Phase 0 验收标准

必须满足：

1. seed 可复现
    
2. perturbation 有方向性效应
    
3. resource constraint 成立
    
4. N(t) 影响观测
    
5. 可生成 dataset
    

---

# 5. Complexity Axes（冻结）

Phase 0 不引入：

- chromatin dynamics
    
- burst kinetics
    
- complex noise
    
- combinatorial perturbation
    

---

# 🔒 最终状态说明

```text
✔ 动力学系统：冻结
✔ world 生成：冻结
✔ 扰动接口：已升级（intervention-ready）
✔ 执行语义：完全闭环
```

---
