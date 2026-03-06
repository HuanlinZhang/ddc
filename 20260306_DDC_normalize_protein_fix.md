# DDC 模型最新进展与关键修改说明

**更新时间**: 2026-03-06
**项目**: Designed Digital Cell (DDC)

---

## 1. 原模型数据衰减问题的根因分析

在原始实现中，生成的轨迹数据显示所有基因的 X（mRNA）和 P（蛋白）值会随时间快速衰减趋近于 0。这并非代码 Bug，而是模型参数与数学公式相互作用的**物理必然结果**。

### 1.1 三个根本原因

#### 原因一：极高的 mRNA 衰减率

mRNA 更新公式前半部分为纯衰减项：

```
X_i(t+1) = (1 - δ_x,i) * X_i(t) + ...
```

参数初始化设定 `δ_x,i ~ Uniform(0.1, 0.5)`，即每个时间步会损失 10%~50% 的 mRNA。在没有任何新生产的情况下，起始值 `X_0 ~ Uniform(0, 1)` 会在几十步内归零。

#### 原因二：蛋白归一化导致信号稀释

蛋白归一化公式：

```
P̃_j = P_j / (ΣP_k + ε)
```

系统中共有 G=50 个基因，蛋白质资源被 50 个基因瓜分。平均每个转录因子的有效浓度仅为 **~0.02**（1/50）。

#### 原因三：希尔函数阈值过高

转录激活项核心是 Hill 函数：

```
TF_input^n / (K^n + TF_input^n)
```

- 输入信号仅 ~0.02
- 激活阈值 K ~ Uniform(0.1, 1.0)，平均 0.55
- 代入计算：0.02² / (0.55² + 0.02²) ≈ 0.0013

**结论**：由于信号远低于激活阈值，转录调控实际处于"关闭"状态。生产项远小于衰减项，导致数据归零。

### 1.2 解决方案

采用**平均场缩放（Mean-field scaling）**方法，在源头统一放大 tilde_P：

```python
def normalize_protein(P: Tensor, world: World) -> Tensor:
    # 核心修改：在归一化后乘以 G，抵消网络规模带来的稀释效应
    return (P / (torch.sum(P) + world.epsilon)) * G
```

此修改将平均输入信号从 0.02 拉回 ~1.0，使其能够跨越 Hill 函数的激活阈值，激活下游转录。

---

## 2. 初始状态增加 Protein Projection

### 2.1 原始实现

在 `sample_initial_state` 函数中，原始代码直接使用 γ * X0 作为初始蛋白 P0：

```python
# 原始实现
P0_raw: Tensor = world.gamma * X0
P0: Tensor = P0_raw  # 未做资源投影
```

### 2.2 当前实现

为保证资源守恒与物理一致性，在初始化时增加了 `apply_resource_projection`：

```python
# 当前实现
P0_raw: Tensor = world.gamma * X0
P0: Tensor = apply_resource_projection(P0_raw, world)
```

这一修改确保：
- 初始蛋白总量不超过 R_total
- 与时间循环中的资源投影逻辑保持一致

---

## 3. 当前代码状态

- **主文件**: `ddc.py`
- **测试文件**: `test_stages.py`
- **分析脚本**: `analyze_vector.ipynb`
- **数据目录**: `./data/`
- **测试目录**: `./test/`

所有修改已同步至 GitHub 仓库。

---

## 4. 下一步计划

- Stage D: Single-gene KO test（基因敲除扰动实验）
- 验证修改后的模型在不同扰动条件下的动力学表现
