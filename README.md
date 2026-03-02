# DDC (Designed Digital Cell)

A computational simulation framework for gene regulatory network dynamics.

## Overview

DDC (Designed Digital Cell) is a computational framework for simulating gene regulatory network (GRN) dynamics at the single-cell level. It supports Monte Carlo sampling, multi-cell dataset generation, and gene perturbation experiments.

### Core Features

- Gene regulatory network simulation (50 genes)
- Chromatin dynamics modeling
- Single-cell mRNA/protein trajectory simulation
- Multi-cell dataset generation (scRNA-seq style)
- Gene perturbation experiments (knockout, overexpression)

## Technical Architecture

### Core Concepts

| Concept | Description |
|---------|-------------|
| **G** | Number of genes (default: 50) |
| **T** | Time steps (default: 200) |
| **World** | Simulation world object |
| **X** | mRNA expression levels |
| **P** | Protein levels |
| **Z** | Chromatin states (0-1) |
| **N** | Cell count |

### Gene Categories

| Category | Indices | Description |
|----------|---------|-------------|
| **TF (Transcription Factors)** | 0-5 | Regulatory proteins |
| **EPI (Epigenetic)** | 17-19 | Chromatin modifiers |
| **Target Genes** | Others | Regulated targets |

### Data Flow

```
t → t+1:

X(t), P(t), Z(t), N(t)
    │
    ▼
normalize_protein: P̃ = P / ΣP
    │
    ▼
compute_TFinput: TF = (∏P̃ᵃʲ)^(1/dᵢ)
    │
    ▼
update_chromatin: Z = σ(α + ΣβᵢⱼP̃ⱼ)
    │
    ▼
update_mRNA: X' = (1-δₓ)X + Z·ρ·hill(TF)
    │
    ▼
update_protein: P_raw = (1-δₚ)P + γX
    │
    ▼
resource_projection: if ΣP > R_total, scale proportionally
    │
    ▼
update_fate: N' = N + r·N·(1-N/K_pop)
    │
    ▼
X(t+1), P(t+1), Z(t+1), N(t+1)
```

## Model Parameters

### Gene-level Parameters

| Parameter | Symbol | Distribution/Value | Description |
|-----------|--------|-------------------|-------------|
| **alpha** | α | N(0, 1) | Basal chromatin activation |
| **rho** | ρ | U(0.5, 2.0) | mRNA production rate |
| **K** | K | U(0.1, 1.0) | Hill coefficient denominator |
| **n** | n | 2.0 (fixed) | Hill coefficient |
| **delta_x** | δₓ | U(0.1, 0.5) | mRNA degradation rate |
| **delta_p** | δₚ | U(0.05, 0.3) | Protein degradation rate |
| **gamma** | γ | 1.0 (fixed) | Translation rate |

### Edge-level Parameters

| Parameter | Symbol | Distribution | Description |
|-----------|--------|--------------|-------------|
| **a_ij** | aᵢⱼ | U(0.5, 2.0) | TF regulatory strength |
| **beta_ij** | βᵢⱼ | N(0, 1.5) | Epigenetic regulatory strength |

### Run-level Parameters

| Parameter | Symbol | Distribution/Value | Description |
|-----------|--------|-------------------|-------------|
| **r** | r | U(0.05, 0.2) | Cell growth rate |
| **K_pop** | K_pop | 1.0 (fixed) | Carrying capacity |
| **R_total** | R_total | 1.0 (fixed) | Total protein resource |

### Network Topology

| Graph | Source | Target | Degree |
|-------|--------|--------|--------|
| **P_graph** | TF (0-5) | All genes (≠ itself) | 1-3 per gene |
| **E_graph** | EPI (17-19) | All genes | 2 per gene |

## Core Functions

### sample_world

```python
def sample_world(seed: int) -> World
```

Generate a random gene network world with all parameters sampled from the distributions defined above. The world contains:
- Gene-level parameters (alpha, rho, K, n, delta_x, delta_p, gamma)
- Edge-level parameters (a_ij, beta_ij)
- Network topology (P_graph, E_graph)

### to_dict / from_dict

```python
def to_dict(self) -> Dict[str, Any]
def from_dict(self, data: Dict[str, Any]) -> None
```

Serialize/deserialize World object for reproducibility. Includes all graph adjacency lists, parameter arrays, and random seed.

### simulate_single_cell

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

Simulate single cell time evolution from initial state to t_steps. Returns:
- `X_traj`: Tensor of shape (t_steps+1, G) - mRNA trajectories
- `P_traj`: Tensor of shape (t_steps+1, G) - protein trajectories
- `Z_traj`: Tensor of shape (t_steps+1, G) - chromatin state trajectories
- `N_traj`: Tensor of shape (t_steps+1,) - cell count trajectory

Convention: index 0 stores initial state at t=0.

### sample_initial_state

```python
def sample_initial_state(cell_seed: int, world: World) -> Tuple[Tensor, Tensor, Tensor, float]
```

Sample initial state for a single cell:
- X0 ~ U(0, 1) - random initial mRNA
- P0 = gamma * X0 with resource projection
- Z0 = sigmoid(alpha) - chromatin initial state
- N0 = 1.0 - initial cell count

### generate_dataset

```python
def generate_dataset(
    world_seed: int,
    M: int,
    save_path: str = None
) -> Tuple[Tensor, World]
```

Generate multi-cell dataset. Returns:
- Expression matrix: Tensor of shape (M, G)
- World object

When `save_path` is provided, data is automatically saved to `save_path` containing:
- `expression`: Expression matrix (M, G)
- `world`: Serialized World object

### run_simulation

```python
def run_simulation(
    seed: int,
    save_path: str = None
) -> Dict[str, Tensor]
```

Run single-cell simulation. Returns:
- `X_traj`: Tensor of shape (T+1, G) - mRNA trajectory
- `P_traj`: Tensor of shape (T+1, G) - protein trajectory
- `Z_traj`: Tensor of shape (T+1, G) - chromatin state trajectory
- `N_traj`: Tensor of shape (T+1,) - cell count trajectory

When `save_path` is provided, data is automatically saved containing:
- `X_traj`, `P_traj`, `Z_traj`, `N_traj`: Trajectory tensors
- `world`: Serialized World object

### apply_perturbation

```python
def apply_perturbation(
    world: World,
    state: State,
    config: Dict[str, Any]
) -> Tuple[World, State]
```

Apply gene perturbation. Config options:
- `knockout`: List of gene indices
- `override_rho`: Dict {gene: value}
- `override_a_ij`: List of (from, to, value)
- `override_alpha`: Dict {gene: value}
- `R_total`: float

### run_smoke_test

```python
def run_smoke_test(seed: int, T: int = 10) -> Dict[str, Tensor]
```

Quick sanity check with T time steps. Tests basic functionality.

### run_sanity_tests

```python
def run_sanity_tests(seed: int) -> None
```

Comprehensive sanity tests:
1. Reproducibility: same seed produces identical results
2. Non-negativity: X, P >= 0, Z in [0, 1]
3. Resource bound: sum(P) <= R_total
4. Stability: all values finite for T=200

## Quick Start

```python
from ddc import run_simulation

result = run_simulation(seed=42)
```

## Version History

| Version | Date | Description |
|---------|------|-------------|
| v1.0 | 2026-03-02 | Initial release, Phase 0 spec |
