# DDC 参数手册

本文档分类整理了 DDC (Designed Digital Cell) 模型中的全部参数及其分布情况。

---

## 1. 全局常量 (Global Constants)

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `G` | `int` | 50 | 基因总数 |
| `T` | `int` | 200 | 模拟时间步数 |
| `R_TOTAL` | `float` | 1.0 | 蛋白质总资源量（归一化基准） |
| `EPSILON` | `float` | 1e-8 | 数值稳定性常数（防止除零） |
| `K_POP` | `float` | 1.0 | 群体承载容量（Logistic 增长参数） |
| `DTYPE` | `torch.dtype` | `torch.float64` | 张量数据类型 |

---

## 2. 基因类别索引 (Gene Category Indices)

### 2.1 核心调控基因 (Core Regulatory Genes)

| 类别 | 索引范围 | 基因数量 |
|------|----------|----------|
| TF (转录因子) | 0-5 | 6 |
| RBP (RNA结合蛋白) | 6-10 | 5 |
| Kinase (激酶) | 11-13 | 3 |
| Phosphatase (磷酸酶) | 14-16 | 3 |
| Epigenetic (表观遗传) | 17-19 | 3 |

### 2.2 细胞命运基因 (Cell Fate Genes)

| 类别 | 索引范围 | 基因数量 |
|------|----------|----------|
| Cell Cycle (细胞周期) | 20-22 | 3 |
| Apoptosis (凋亡) | 23-25 | 3 |

### 2.3 背景基因 (Background Genes)

| 类别 | 索引范围 | 基因数量 |
|------|----------|----------|
| Background | 26-49 | 24 |

---

## 3. World 对象参数 (World Class Parameters)

### 3.1 基因级参数 (Gene-level Parameters)

| 参数 | 类型 | 维度 | 分布/取值 | 说明 |
|------|------|------|-----------|------|
| `alpha` | `Tensor` | (G,) | `N(0, 1)` | 表观遗传 baseline（Z 公式截距） |
| `rho` | `Tensor` | (G,) | `U(0.5, 2.0)` | 转录速率系数 |
| `K` | `Tensor` | (G,) | `U(0.1, 1.0)` | Hill 函数阈值 |
| `n` | `Tensor` | (G,) | `2.0` (固定) | Hill 系数 |
| `delta_x` | `Tensor` | (G,) | `U(0.1, 0.5)` | mRNA 衰减率 |
| `delta_p` | `Tensor` | (G,) | `U(0.05, 0.3)` | 蛋白衰减率 |
| `gamma` | `Tensor` | (G,) | `1.0` (固定) | 翻译速率系数 |

### 3.2 边级参数 (Edge-level Parameters)

| 参数 | 类型 | 维度 | 分布/取值 | 说明 |
|------|------|------|-----------|------|
| `P_graph` | `Dict[int, List[int]]` | - | 每个基因 i 连接 1-3 个 TF | 转录调控图：P_graph[i] = 受基因 i 调控的上游 TF 列表 |
| `a_ij` | `Dict[int, Dict[int, float]]` | - | `U(0.5, 2.0)` | 转录调控强度：基因 j 对基因 i 的调控系数 |
| `E_graph` | `Dict[int, List[int]]` | - | 每个基因 i 连接 2 个表观基因 | 染色质调控图：E_graph[i] = 影响基因 i 的表观基因列表 |
| `beta_ij` | `Dict[int, Dict[int, float]]` | - | `N(0, 1.5)` | 染色质调控强度：表观基因 j 对基因 i 的调控系数 |

### 3.3 运行级参数 (Run-level Parameters)

| 参数 | 类型 | 分布/取值 | 说明 |
|------|------|-----------|------|
| `r` | `float` | `U(0.05, 0.2)` | 群体增长率（Logistic） |
| `K_pop` | `float` | `1.0` | 群体承载容量 |

### 3.4 全局常数引用

| 参数 | 类型 | 说明 |
|------|------|------|
| `R_total` | `float` | 蛋白质总资源量（默认 1.0） |
| `epsilon` | `float` | 数值稳定性常数（默认 1e-8） |

### 3.5 基因注释

| 参数 | 类型 | 说明 |
|------|------|------|
| `gene_categories` | `Dict` | 基因类别层级结构 |
| `gene_to_macro` | `Dict[int, str]` | 基因索引 → 大类映射 |
| `gene_to_micro` | `Dict[int, str]` | 基因索引 → 细分类映射 |

---

## 4. 状态变量 (State Variables)

在模拟过程中，每个时间步包含以下状态：

| 变量 | 类型 | 形状 | 说明 |
|------|------|------|------|
| `X` | `Tensor` | (G,) | mRNA 表达量 |
| `P` | `Tensor` | (G,) | 蛋白水平 |
| `Z` | `Tensor` | (G,) | 染色质开放状态 (0-1) |
| `N` | `float` | scalar | 细胞数量/群体状态 |

---

## 5. 输出轨迹 (Output Trajectories)

| 变量 | 类型 | 形状 | 说明 |
|------|------|------|------|
| `X_traj` | `Tensor` | (T+1, G) | mRNA 轨迹，index 0 为初始状态 |
| `P_traj` | `Tensor` | (T+1, G) | 蛋白轨迹 |
| `Z_traj` | `Tensor` | (T+1, G) | 染色质状态轨迹 |
| `N_traj` | `Tensor` | (T+1,) | 群体状态轨迹 |

---

*文档版本: v1.1*
*更新时间: 2026-03-06*
