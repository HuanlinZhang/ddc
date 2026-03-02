# DDC (Designed Digital Cell)

[![PyTorch](https://img.shields.io/badge/PyTorch-1.9+-ee4c2c?style=flat&logo=pytorch)](https://pytorch.org)
[![Python](https://img.shields.io/badge/Python-3.8+-3776ab?style=flat&logo=python)](https://www.python.org)
[![Status](https://img.shields.io/badge/Status-Frozen-blue?style=flat)]

基因调控网络动力学计算模拟框架

## 项目概述

DDC（Designed Digital Cell）是一个基于计算模拟的基因调控网络（Gene Regulatory Network, GRN）动力学模拟框架。该项目实现了单细胞水平的基因表达动态模拟，支持蒙特卡洛采样、多细胞数据集生成以及基因扰动实验。

### 核心功能

- **基因调控网络模拟**：模拟50个基因之间的转录调控关系
- **染色质动力学**：实现表观遗传修饰对基因表达的调控
- **单细胞模拟**：模拟单个细胞内mRNA、蛋白质和染色质状态的时间演化
- **多细胞数据集生成**：支持生成大规模单细胞RNA测序（scRNA-seq）风格的数据集
- **基因扰动实验**：支持基因敲除（knockout）、参数过表达等扰动操作

### 应用场景

- 研究基因调控网络的动态行为与稳定性
- 单细胞RNA测序数据分析与建模
- 细胞命运分化与分化轨迹研究
- 基因功能预测与验证
- 合成生物学回路设计

---

## 技术架构

### 核心概念

| 概念 | 描述 |
|------|------|
| **G** | 基因数量，默认50个 |
| **T** | 时间步长，默认200步 |
| **World** | 模拟世界对象，存储所有参数和调控图 |
| **X** | mRNA表达水平向量 |
| **P** | 蛋白质水平向量 |
| **Z** | 染色质状态向量（0-1之间） |
| **N** | 细胞数量 |

### 基因分类

| 类别 | 基因索引 | 描述 |
|------|----------|------|
| **转录因子基因 (TF)** | 0-5 | 调控其他基因表达的蛋白质 |
| **表观遗传基因 (EPI)** | 17-19 | 参与染色质修饰的基因 |
| **普通基因** | 其他 | 被调控的靶基因 |

### 核心模块

```
ddc.py
├── 全局常量定义
│   ├── G = 50              # 基因数量
│   ├── T = 200             # 时间步长
│   ├── R_TOTAL = 1.0       # 蛋白质资源总量
│   └── DTYPE = torch.float64
│
├── World 类
│   ├── 调控图结构
│   │   ├── P_graph         # 转录调控图
│   │   └── E_graph         # 染色质调控图
│   ├── 基因级参数
│   │   ├── alpha           # 染色质基础活性
│   │   ├── rho             # 转录速率
│   │   ├── K               # Hill系数常数
│   │   ├── n               # Hill系数
│   │   ├── delta_x         # mRNA降解率
│   │   ├── delta_p         # 蛋白质降解率
│   │   └── gamma           # 翻译效率
│   ├── 边级参数
│   │   ├── a_ij            # 转录调控权重
│   │   └── beta_ij         # 染色质调控权重
│   └── 运行参数
│       ├── r               # 细胞增长率
│       └── K_pop           # 环境容纳量
│
├── 蒙特卡洛采样器
│   └── sample_world()      # 随机生成基因网络世界
│
├── 模块接口函数
│   ├── normalize_protein()        # 蛋白质归一化
│   ├── compute_TFinput()          # 计算转录因子输入
│   ├── update_chromatin()         # 更新染色质状态
│   ├── update_mRNA()              # 更新mRNA水平
│   ├── update_protein_raw()       # 更新蛋白质（原始）
│   ├── apply_resource_projection() # 资源约束投影
│   └── update_fate()              # 更新细胞数量
│
├── 主模拟循环
│   └── simulate_single_cell()     # 单细胞时间演化
│
├── 数据生成器
│   ├── sample_initial_state()     # 采样初始状态
│   ├── run_simulation()           # 运行完整模拟
│   └── generate_dataset()         # 生成多细胞数据集
│
└── 扰动接口
    └── apply_perturbation()       # 应用基因扰动
```

### 数据流

```
┌─────────────────────────────────────────────────────────────────┐
│                        时间步 t → t+1                            │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  X(t), P(t), Z(t), N(t)                                        │
│         │                                                       │
│         ▼                                                       │
│  ┌─────────────────┐                                           │
│  │ normalize_protein│  P̃ = P / ΣP                              │
│  └────────┬────────┘                                           │
│           │                                                     │
│           ▼                                                     │
│  ┌─────────────────┐                                           │
│  │ compute_TFinput │  TF = (∏P̃ᵃʲ)^(1/dᵢ)                      │
│  └────────┬────────┘                                           │
│           │                                                     │
│           ▼                                                     │
│  ┌─────────────────┐                                           │
│  │update_chromatin│  Z = σ(α + ΣβᵢⱼP̃ⱼ)                        │
│  └────────┬────────┘                                           │
│           │                                                     │
│           ▼                                                     │
│  ┌─────────────────┐                                           │
│  │  update_mRNA   │  X' = (1-δₓ)X + Z·ρ·hill(TF)              │
│  └────────┬────────┘                                           │
│           │                                                     │
│           ▼                                                     │
│  ┌─────────────────┐                                           │
│  │update_protein  │  P_raw = (1-δₚ)P + γX                      │
│  └────────┬────────┘                                           │
│           │                                                     │
│           ▼                                                     │
│  ┌─────────────────┐                                           │
│  │resource_project │  若 ΣP_raw > R_total: P' = P_raw·R/ΣP   │
│  └────────┬────────┘              否则: P' = P_raw              │
│           │                                                     │
│           ▼                                                     │
│  ┌─────────────────┐                                           │
│  │  update_fate   │  N' = N + r·N·(1-N/K_pop)                  │
│  └────────┬────────┘                                           │
│           │                                                     │
│           ▼                                                     │
│  X(t+1), P(t+1), Z(t+1), N(t+1)                                │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## 环境要求

### 运行环境

- **Python**: 3.8 或更高版本
- **PyTorch**: 1.9.0 或更高版本

### 依赖项

| 包名 | 最低版本 | 用途 |
|------|----------|------|
| torch | 1.9.0 | 张量计算与神经网络支持 |
| (标准库) json | - | 数据序列化 |
| (标准库) copy | - | 对象深拷贝 |

---

## 安装部署

### 1. 安装Python

确保已安装Python 3.8+：

```bash
python --version
```

### 2. 安装PyTorch

```bash
# CPU版本
pip install torch

# GPU版本（需要CUDA）
pip install torch --index-url https://download.pytorch.org/whl/cu118
```

### 3. 克隆项目

```bash
git clone https://github.com/HuanlinZhang/ddc.git
cd ddc
```

### 4. 验证安装

```bash
python -c "import torch; print(f'PyTorch version: {torch.__version__}')"
```

---

## 使用指南

### 基础示例

#### 1. 运行烟雾测试

```python
from ddc import run_smoke_test

traj = run_smoke_test(seed=42, T=10)
```

#### 2. 运行完整模拟

```python
from ddc import run_simulation

result = run_simulation(seed=42)
# 返回包含 X_traj, P_traj, Z_traj, N_traj 的字典
```

#### 3. 生成多细胞数据集

```python
from ddc import generate_dataset

M = 100  # 细胞数量
dataset, world = generate_dataset(world_seed=42, M=M)
# dataset shape: (M, G) = (100, 50)
```

### 高级示例

#### 使用自定义World对象

```python
from ddc import sample_world, simulate_single_cell, sample_initial_state
import torch

# 创建世界
world = sample_world(seed=42)

# 采样初始状态
X0, P0, Z0, N0 = sample_initial_state(cell_seed=43, world=world)

# 运行模拟
result = simulate_single_cell(world, X0, P0, Z0, N0, t_steps=200)
```

#### 应用基因扰动

```python
from ddc import sample_world, sample_initial_state, simulate_single_cell, apply_perturbation

# 创建世界和初始状态
world = sample_world(seed=42)
X0, P0, Z0, N0 = sample_initial_state(cell_seed=43, world=world)

# 定义扰动配置
config = {
    'knockout': [0, 1],           # 敲除基因0和1
    'override_rho': {2: 5.0},     # 过表达基因2的转录速率
    'override_a_ij': [(3, 4, 0.0)], # 移除基因3对4的调控
    'override_alpha': {5: -2.0},  # 修改基因5的染色质基础活性
    'R_total': 0.8,               # 降低资源总量
}

# 应用扰动
world_pert, state_pert = apply_perturbation(world, 
    {'X': X0, 'P': P0, 'Z': Z0, 'N': N0}, config)

# 运行扰动后的模拟
result = simulate_single_cell(world_pert, 
    state_pert['X'], state_pert['P'], state_pert['Z'], state_pert['N'])
```

#### 保存和加载World对象

```python
import json
from ddc import sample_world

# 保存
world = sample_world(seed=42)
world_dict = world.to_dict()
with open('world.json', 'w') as f:
    json.dump(world_dict, f)

# 加载
with open('world.json', 'r') as f:
    loaded_dict = json.load(f)
world = World(0)
world.from_dict(loaded_dict)
```

---

## API文档

### 全局常量

| 常量 | 类型 | 默认值 | 描述 |
|------|------|--------|------|
| `G` | int | 50 | 基因数量 |
| `T` | int | 200 | 默认时间步长 |
| `R_TOTAL` | float | 1.0 | 蛋白质资源总量 |
| `EPSILON` | float | 1e-8 | 数值稳定性 epsilon |
| `K_POP` | float | 1.0 | 环境容纳量 |
| `DTYPE` | torch.dtype | torch.float64 | 浮点精度 |
| `TF_GENES` | List[int] | [0,1,2,3,4,5] | 转录因子基因索引 |
| `EPI_GENES` | List[int] | [17,18,19] | 表观遗传基因索引 |

### World 类

存储模拟世界的所有参数和结构。

#### 构造函数

```python
World(seed: int)
```

#### 主要属性

| 属性 | 类型 | 描述 |
|------|------|------|
| `seed` | int | 随机种子 |
| `P_graph` | Dict[int, List[int]] | 转录调控图 |
| `E_graph` | Dict[int, List[int]] | 染色质调控图 |
| `alpha` | Tensor | 染色质基础活性 (G,) |
| `rho` | Tensor | 转录速率 (G,) |
| `K` | Tensor | Hill常数 (G,) |
| `n` | Tensor | Hill系数 (G,) |
| `delta_x` | Tensor | mRNA降解率 (G,) |
| `delta_p` | Tensor | 蛋白质降解率 (G,) |
| `gamma` | Tensor | 翻译效率 (G,) |
| `a_ij` | Dict | 转录调控边权重 |
| `beta_ij` | Dict | 染色质调控边权重 |
| `r` | float | 细胞增长率 |
| `K_pop` | float | 种群容纳量 |

#### 方法

```python
def to_dict() -> Dict[str, Any]:
    """将World对象序列化为字典"""

def from_dict(cls, data: Dict[str, Any]) -> World:
    """从字典反序列化创建World对象"""
```

### 核心函数

#### sample_world

```python
def sample_world(seed: int) -> World
```

随机生成一个基因网络世界。

**参数：**
- `seed` (int): 随机种子

**返回：**
- `World`: 包含随机参数的模拟世界

**示例：**
```python
world = sample_world(seed=42)
```

#### simulate_single_cell

```python
def simulate_single_cell(
    world: World,
    X0: Tensor,
    P0: Tensor,
    Z0: Tensor,
    N0: float,
    t_steps: int = T
) -> Dict[str, Tensor]
```

模拟单个细胞的时间演化。

**参数：**
- `world` (World): 模拟世界对象
- `X0` (Tensor): 初始mRNA水平 (G,)
- `P0` (Tensor): 初始蛋白质水平 (G,)
- `Z0` (Tensor): 初始染色质状态 (G,)
- `N0` (float): 初始细胞数量
- `t_steps` (int): 模拟时间步数

**返回：**
```python
{
    'X_traj': Tensor,  # shape: (t_steps+1, G)
    'P_traj': Tensor,  # shape: (t_steps+1, G)
    'Z_traj': Tensor,  # shape: (t_steps+1, G)
    'N_traj': Tensor   # shape: (t_steps+1,)
}
```

**注意：** 索引0存储t=0时刻的初始状态。

#### generate_dataset

```python
def generate_dataset(world_seed: int, M: int) -> Tuple[Tensor, World]
```

生成多细胞数据集。

**参数：**
- `world_seed` (int): 世界随机种子
- `M` (int): 细胞数量

**返回：**
- `Tensor`: 细胞表达矩阵 (M, G)
- `World`: 对应的世界对象

**示例：**
```python
dataset, world = generate_dataset(world_seed=42, M=100)
# dataset.shape = (100, 50)
```

#### apply_perturbation

```python
def apply_perturbation(
    world: World,
    state: State,
    config: Dict[str, Any]
) -> Tuple[World, State]
```

应用基因扰动。

**参数：**
- `world` (World): 原始世界对象
- `state` (State): 当前状态字典 {'X': X, 'P': P, 'Z': Z, 'N': N}
- `config` (Dict): 扰动配置

**扰动配置选项：**
```python
config = {
    'knockout': [0, 1, 2],           # 基因敲除列表
    'override_rho': {0: 5.0},        # 转录速率覆盖 {gene: value}
    'override_a_ij': [(0, 1, 0.0)],  # 转录调控权重覆盖 (from, to, value)
    'override_alpha': {0: -2.0},     # 染色质基础活性覆盖 {gene: value}
    'R_total': 0.8,                  # 资源总量修改
}
```

**返回：**
- `Tuple[World, State]`: 扰动后的世界和状态

### 辅助函数

#### stable_sigmoid

```python
def stable_sigmoid(x: Tensor) -> Tensor
```

数值稳定的Sigmoid函数。

#### normalize_protein

```python
def normalize_protein(P: Tensor, world: World) -> Tensor
```

蛋白质归一化：`P̃ = P / (ΣP + ε)`

#### compute_TFinput

```python
def compute_TFinput(tilde_P: Tensor, world: World) -> Tensor
```

计算转录因子输入：`TFᵢ = (∏ⱼ∈P(i) P̃ⱼᵃⁱʲ)^(1/dᵢ)`

#### update_chromatin

```python
def update_chromatin(tilde_P: Tensor, world: World) -> Tensor
```

更新染色质状态：`Zᵢ = σ(αᵢ + Σⱼ∈E(i) βᵢⱼP̃ⱼ)`

#### update_mRNA

```python
def update_mRNA(X: Tensor, Z: Tensor, TFinput: Tensor, world: World) -> Tensor
```

更新mRNA水平：`X' = (1-δₓ)X + Z·ρ·hill(TF)`

#### update_protein_raw

```python
def update_protein_raw(P: Tensor, X: Tensor, world: World) -> Tensor
```

更新蛋白质（未应用资源约束）：`P' = (1-δₚ)P + γX`

#### apply_resource_projection

```python
def apply_resource_projection(P_raw: Tensor, world: World) -> Tensor
```

应用资源约束投影，确保 ΣP ≤ R_total

#### update_fate

```python
def update_fate(N: float, world: World) -> float
```

更新细胞数量：`N' = N + r·N·(1-N/K_pop)`

---

## 测试与验证

### 运行烟雾测试

```bash
python ddc.py
```

输出：
```
Running DDC Phase 0 Standard Pipeline...

--- T=10 Smoke Test ---
Running smoke test with T=10...
  X shape: torch.Size([11, 50])
  P shape: torch.Size([11, 50])
  Z shape: torch.Size([11, 50])
  N shape: torch.Size([11])
Smoke test passed!

--- T=200 Stability Test ---
Running sanity tests...
Reproducibility check passed.
Non-negativity and Resource bound checks passed.
Stability check passed.
All sanity tests passed!

--- Multi-cell Dataset ---
Dataset generated successfully: torch.Size([100, 50])
```

### 测试函数

#### run_smoke_test

```python
def run_smoke_test(seed: int, T: int = 10) -> Dict[str, Tensor]
```

快速验证基本功能是否正常。

#### run_sanity_tests

```python
def run_sanity_tests(seed: int) -> None
```

运行完整性测试，包括：
1. **可重复性**：相同种子产生相同结果
2. **非负性**：X, P ≥ 0，Z ∈ [0, 1]
3. **资源约束**：ΣP ≤ R_total
4. **数值稳定性**：所有值在T=200步内保持有限

---

## 参数配置

### 基因参数范围

| 参数 | 分布 | 范围 |
|------|------|------|
| alpha | 正态分布 | N(0, 1) |
| rho | 均匀分布 | [0.5, 2.0] |
| K | 均匀分布 | [0.1, 1.0] |
| n | 常数 | 2.0 |
| delta_x | 均匀分布 | [0.1, 0.5] |
| delta_p | 均匀分布 | [0.05, 0.3] |
| gamma | 常数 | 1.0 |
| r | 均匀分布 | [0.05, 0.2] |

### 调控图结构

- 每个基因受 **1-3个** 转录因子调控
- 每个基因受 **2个** 表观遗传基因调控
- 转录调控权重 a_ij ∈ [0.5, 2.0]
- 染色质调控权重 beta_ij ∈ N(0, 1.5)

---

## 版本历史

| 版本 | 日期 | 状态 | 描述 |
|------|------|------|------|
| v1.0 | 2026-03-02 | Frozen | 初始版本，Phase 0规范实现 |

---

## 贡献规范

### 提交问题

如果您发现bug或有功能请求，请提交Issue。请包含：

- 清晰的问题描述
- 复现步骤
- 环境信息（Python版本、PyTorch版本）
- 预期行为 vs 实际行为

### 代码贡献

1. Fork本仓库
2. 创建功能分支 (`git checkout -b feature/amazing-feature`)
3. 提交更改 (`git commit -m 'Add amazing feature'`)
4. 推送分支 (`git push origin feature/amazing-feature`)
5. 创建Pull Request

### 开发建议

- 遵循现有代码风格
- 添加单元测试覆盖新功能
- 更新文档以反映代码变更

---

## 许可证

本项目未指定明确的许可证。使用前请联系作者获取授权。

---

## 联系方式

- 作者：zhanghl
- 项目地址：https://github.com/HuanlinZhang/ddc

---

## 参考文献与延伸阅读

1. **基因调控网络 (GRN)**：描述基因之间相互作用的数学模型
2. **Hill函数**：用于建模转录调控的非线性动力学
3. **染色质动力学**：表观遗传修饰对基因表达的调控机制
4. **单细胞RNA测序 (scRNA-seq)**：高通量测量单个细胞基因表达的技术
5. **资源约束模型**：模拟细胞内有限资源对蛋白质合成的限制

---

*本项目基于 Phase 0 规范文档实现*
