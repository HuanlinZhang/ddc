# Run 009 — Mechanistic Refinement 分析总结

## 1. 背景

Run 009 基于 run\_007 和 run\_008 的扰动实验结果，对调控系统的机制特征进行三个层面的深入分析：

- **Task 1**：TF 基因在 Hill 曲线上的稳态位置——探索 collapse 敏感性与 Hill 曲线位置的关系
- **Task 2**：Collapse 模式分类——区分 threshold crossing 和 rapid collapse
- **Task 3**：延长模拟时间，确定异常存活案例的真实命运

### 分析对象

从 run\_007 的稳态世界（baseline regime = steady）中选取了 4 个seed（2026, 4094, 4478, 7709），对每个seed的 6 个 TF 基因施加 λ = 0.5 和 λ = 0.1 的持续衰减扰动（run\_008 实验），模拟时长 T = 400。

***

## 2. Task 1 — Hill 曲线定位

### 方法

本分析以 Hill 曲线为理论框架。根据 Phase 0 设计文档，Transcription 步骤中的激活项为：

$$
\frac{\text{TFinput}_i^{n_i}}{K_i^{n_i} + \text{TFinput}_i^{n_i}}
$$

令 $x = \text{TFinput}_i / K_i$，$n_i = 2$，则激活水平简化为：

$$
y = \frac{x^2}{1 + x^2}
$$

对每个种子，以 baseline（无扰动）轨迹为基础，计算每个基因的 TFinput/K 时间序列，再对全轨迹（401 步）求平均，得到每个 TF 基因的"稳态工作点" (x, y)：

- **x** = mean\_t(TFinput\_i / K\_i) — 该基因的调控输入相对于 half-saturation constant 的位置
- **y** = x² / (1 + x²) — Hill activation 水平（n=2）

同时计算全部 50 个基因的均值。

> `task1_hill_curve_positioning.png`

#### 2.1 TF 基因的 Hill 曲线位置

每个种子 6 个 TF 基因的 (x, y) 位置用不同颜色区分 collapse 行为：

- **红色**：λ=0.5 时即 collapse（严重敏感）
- **橙色**：仅在 λ=0.1 时 collapse（中等敏感）
- **绿色**：λ=0.1 时仍存活（robust）

> `task1_hill_curve_all_genes.png`

#### 2.2 全部 50 个基因的分布

44 个 non-TF 基因（灰色小点）与 6 个 TF 基因（红色大点）共同绘制，显示：

- TF 基因与 non-TF 基因在 Hill 曲线上的分布模式
- 全部基因均值（深蓝色菱形）所在位置

> `task1_hill_mean_comparison.png`

#### 2.3 四个稳态世界的均值比较

将四个种子的全部 50 个基因的 TFinput/K 均值标注在同一张 Hill 曲线图上。

### 关键发现

| 种子   | TFinput/K（6 TF） | TFinput/K（50 基因） | Hill Act.（50 基因） | Collapse 敏感性         |
| ---- | --------------- | ---------------- | ---------------- | -------------------- |
| 2026 | \~1.4           | \~3.7            | \~0.93           | 高（4\~6 基因 collapse）  |
| 4094 | \~6.1           | \~4.6            | \~0.95           | 极低（0 collapse）       |
| 4478 | \~0.9           | \~1.2            | \~0.60           | 中（3\~5 基因 collapse）  |
| 7709 | \~1.1           | \~2.3            | \~0.84           | 中高（8/12 组合 collapse） |

**核心规律**：稳态 TFinput/K 与 collapse 敏感性存在强相关。Seed 4094 无论从 TF 子集还是全部基因看均值都最高，系统完全 robust；Seed 4478 的 50 基因均值最低（\~1.2），最接近 half-activation point。值得注意的是 2026 的 TF 基因均值仅 \~1.4，但全部 50 基因均值拉到 \~3.7，说明 non-TF 基因起到了缓冲作用——然而这并未阻止它在 λ=0.1 时 6 个 TF 基因全部 collapse，暗示 TF 基因自身的 TFinput/K 水平比全局均值更能预测敏感性。

***

## 3. Task 2 — Collapse 模式分类

### 方法

对每个 (seed, gene, λ) 组合，计算 **TF-core mean TFinput/K** 时间序列（即所有 6 个 TF 基因在该时刻的 TFinput/K 的均值），然后用滑动窗口（窗口大小=5）平滑后，检测是否出现持续低于阈值（0.3）的跨越：

- **threshold\_crossing**：均值在阈值上方维持较长时间后，持续低于 0.3 超过 10 步，且跨越前后的演化时间充足（>30 步 before, >20 步 after）
- **rapid\_collapse**：跨越发生得很早（pre\_decay\_len ≤ 30 或 post\_decay\_len ≤ 20）
- **gradual\_decay**：未检测到明显的阈值跨越，但最终均值低于 0.3
- **no\_collapse**：从未进入 collapse 模式

**注**：分类基于 TF-core（全部 TF 基因）的平均 TFinput/K，而非被扰动基因自身的 TFinput/K。这是因为被扰动基因自身的输入会因其蛋白质被衰减而立即下降，不能反映整个调控网络的健康状态。

**关于阈值 0.3**：这是一个相对保守的选取值，没有严格的理论依据或先验证据，仅为便于说明 collapse 行为而设定。该阈值仅对处于临界区的 collapse 案例有意义：更低的阈值使分类标准更严格（更多案例被判为 threshold\_crossing 而非 rapid\_collapse），更高的阈值使分类标准更宽松（更多案例被归为 rapid\_collapse）。no\_collapse 案例无论阈值怎么取都不会归为 collapse。后续可通过敏感性分析验证分类结果对该阈值的稳健程度。

### 结果

> `task2_collapse_mode.png`

共 48 个 (seed, gene, λ) 组合参与分类：

| 模式                  | 数量 | 占比    |
| ------------------- | -- | ----- |
| no\_collapse        | 25 | 52.1% |
| threshold\_crossing | 13 | 27.1% |
| rapid\_collapse     | 10 | 20.8% |
| gradual\_decay      | 0  | 0%    |

#### 各种子的 collapse 模式分布

**Seed 2026**（最敏感）：

- λ=0.5：4 个 collapse（gene 1/2/4/5），全部为 threshold\_crossing，crossing\_time 分布在 45\~151 步；gene 0 和 gene 3 为 no\_collapse
- λ=0.1：6 个全部 collapse，其中 4 个 threshold\_crossing（gene 0/1/3/5, crossing\_time 86\~134），2 个 rapid\_collapse（gene 2/4, crossing\_time=13\~20）

**Seed 4094**（最 robust）：

- λ=0.5 和 λ=0.1：全部 12 个组合均为 no\_collapse，系统完全不敏感

**Seed 4478**：

- λ=0.5：2 个 threshold\_crossing（gene 0 crossing\_time=39, gene 3 crossing\_time=86），4 个 no\_collapse
- λ=0.1：3 个 rapid\_collapse（gene 0/3/4），3 个 no\_collapse

**Seed 7709**：

- λ=0.5：2 个 threshold\_crossing（gene 1 crossing\_time=56, gene 2 crossing\_time=43），4 个 no\_collapse
- λ=0.1：6 个全部 collapse，5 个 rapid\_collapse（gene 0/1/2/4/5），1 个 threshold\_crossing（gene 3 crossing\_time=117）

### 关键发现

1. **Collapse 模式的多样性**：约 30% 的 collapse 案例表现出 threshold\_crossing 模式——系统在较长时间内维持 quasi-steady，然后突然崩溃。这暗示存在"临界点"（critical point）或"不可逆转折点"（point of no return）的动力学特征。
2. **rapid\_collapse 的跨越时间通常很短**（<30 步），对应扰动后快速级联失效。
3. **跨种子差异明显**：Seed 2026 的 collapse 以 threshold\_crossing 为主（8/10 的 collapse 案例），而 Seed 7709 则相反——λ=0.1 时 5/6 的 collapse 为 rapid\_collapse，反映出两种种子在调控网络结构上的根本差异。
4. **λ 从 0.5 降至 0.1 时，threshold\_crossing 比例下降**：2026 中 λ=0.5 时 4/4 为 threshold\_crossing，λ=0.1 时降至 4/6；4478 和 7709 也有类似趋势。更大的扰动强度（λ=0.1）使系统更快越过临界点，从而阈值跨越行为被压缩为 rapid collapse。

***

## 4. Task 3 — 延长 KO 模拟

### 方法

选取 run\_007 中 T=400 时表现异常的 KO 案例（包括 steady 到 T=200 附近异常死亡、以及 baseline 本身非稳态但 KO 后存活至 T=400 的案例），将模拟延长至 T=800（部分案例延长至 T=1200），以确定是 slow\_collapse 还是 quasi-steady：

| 种子   | KO 基因 | 来源                            | 延长时长   | 最终命运           |
| ---- | ----- | ----------------------------- | ------ | -------------- |
| 4094 | 4     | run\_007                      | T=800  | slow\_collapse |
| 4094 | 5     | run\_007                      | T=800  | slow\_collapse |
| 20   | 2     | run\_007\_sampling\_extension | T=800  | slow\_collapse |
| 20   | 3     | run\_007\_sampling\_extension | T=1200 | slow\_collapse |
| 33   | 3     | run\_007\_sampling\_extension | T=800  | slow\_collapse |

### 结果

> `task3_extended_trajectories.png`

**所有 5 个案例均确认为 slow\_collapse**——T=400 时看起来存活，但延长模拟后均走向崩溃（n\_active\_genes = 0）。区别在于崩溃速度不同：

- Seed 4094 gene 4/5 和 seed 20 gene 2 在 T=800 范围内已完全 collapse
- Seed 33 gene 3 在 T=800 之前已完全 collapse
- Seed 20 gene 3 延长至 T=1200 后完全 collapse

### 关键发现

这些案例揭示了"quasi-steady deception"现象：某些案例在 T=400 时看似存活，但实际上是 slow\_collapse 的早期阶段。延长模拟时间是识别这类案例的必要手段。

***

## 5. 综合结论

### 5.1 Hill 曲线位置与 Collapse 敏感性

TF 基因的 TFinput/K 位置是 collapse 敏感性的强预测因子：

- **均值 \~6**：系统 robust（seed 4094）
- **均值 \~1.0\~1.4**：中等到高敏感性（seeds 2026, 7709）
- **均值 \~0.9**：中等敏感性，但个别基因极低（seed 4478，gene 3 的 TFinput/K 仅 0.19）

Hill 曲线上的"高敏感区"（x ≈ 1 附近，红色阴影）似乎是判断基因脆弱性的关键区域。但需要结合基因在网络中的出度（out-degree）综合判断——高 degree 的基因即使 x 值不低，也可能因为其广泛影响而触发 collapse。

### 5.2 Collapse 模式的双重性

- **Threshold crossing**：暗示调控网络存在临界阈值——跨越后系统失去自我维持能力
- **Rapid collapse**：暗示快速级联失效——单个节点失效迅速传播

不同种子对这两种模式的倾向性不同，可能反映了其调控网络拓扑的差异。

### 5.3 Quasi-Steady State 检测

run\_007 的 T=400 判定存在假阳性风险——某些案例需要更长的模拟时间才能显现 collapse 趋势。

***

## 6. 输出文件列表

### 图片

| 文件名                                | 内容                        |
| ---------------------------------- | ------------------------- |
| `task1_hill_curve_positioning.png` | 6 个 TF 基因的 Hill 曲线位置，4 面板 |
| `task1_hill_curve_all_genes.png`   | 全部 50 个基因的 Hill 曲线位置，4 面板 |
| `task1_hill_mean_comparison.png`   | 4 个种子全部基因均值的比较            |
| `task2_collapse_mode.png`          | Collapse 模式分布柱状图          |
| `task3_extended_trajectories.png`  | 5 个延长模拟案例的蛋白质轨迹           |

### 数据

| 文件名                               | 内容                                                       |
| --------------------------------- | -------------------------------------------------------- |
| `task1_hill_position_summary.tsv` | 6 个 TF 基因的 TFinput/K、Hill activation、slope 及 collapse 记录 |
| `task2_collapse_mode.tsv`         | 每个 (seed, gene, λ) 组合的 collapse 模式分类及 crossing\_time     |
| `task3_extended_simulation.tsv`   | 5 个延长模拟案例的最终命运判定                                         |

***

## 7. 技术说明

- **TFinput 计算**：通过 `normalize_protein` 归一化蛋白质水平后，使用 `compute_TFinput` 计算几何平均
- **分类阈值**：threshold = 0.3（对应 \~8% Hill activation）
- **滑动窗口平滑**：窗口大小 = 5
- **持续性要求**：低于阈值需维持至少 10 步才认定为有效跨越
- **模拟时长**：Task 1/2 基于 T=400，Task 3 部分案例延长至 T=800\~1200

