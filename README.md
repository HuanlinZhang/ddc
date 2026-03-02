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

## Environment Requirements

- **Python**: 3.8+
- **PyTorch**: 1.9.0+ (tested with 2.10.0)
- **NumPy**: 1.20.0+ (tested with 1.23.5)

## Installation

```bash
pip install torch numpy
```

Or clone and install in development mode:

```bash
git clone https://github.com/HuanlinZhang/ddc.git
cd ddc
pip install -e .
```

## Core Functions

### sample_world

```python
def sample_world(seed: int) -> World
```

Generate a random gene network world.

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

Simulate single cell time evolution. Returns:
- `X_traj`: (t_steps+1, G)
- `P_traj`: (t_steps+1, G)
- `Z_traj`: (t_steps+1, G)
- `N_traj`: (t_steps+1,)

### generate_dataset

```python
def generate_dataset(world_seed: int, M: int) -> Tuple[Tensor, World]
```

Generate multi-cell dataset. Returns:
- Expression matrix: (M, G)
- World object

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

## Quick Start

```python
from ddc import run_simulation

result = run_simulation(seed=42)
```

## Version History

| Version | Date | Description |
|---------|------|-------------|
| v1.0 | 2026-03-02 | Initial release, Phase 0 spec |
