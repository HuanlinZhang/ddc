# DDC (Designed Digital Cell)

A computational simulation framework for gene regulatory network (GRN) dynamics.

## Overview

DDC simulates gene regulatory network dynamics at the single-cell level. It supports Monte Carlo sampling, multi-cell dataset generation, and perturbation experiments.

## Installation

```bash
pip install -e .
```

Or from GitHub:

```bash
pip install git+https://github.com/EasyPiPi/ddc.git
```

Requires: Python 3.8+, PyTorch 2.0+, NumPy 1.20+

## CLI

```bash
# Run simulation
ddc run --seed 42 --output ./traj.pt

# Generate multi-cell dataset
ddc dataset --world-seed 0 --M 100 --output dataset.pt

# Sample and save a world
ddc world --seed 42 --output world.json

# Smoke test
ddc smoke-test

# Sanity tests
ddc sanity
```

## Python API

```python
from ddc import run_simulation, generate_dataset, sample_world

# Single-cell simulation
traj = run_simulation(seed=42, save_path='./traj.pt')
# traj: X_traj, P_traj, Z_traj, N_traj

# Multi-cell dataset
dataset, world = generate_dataset(world_seed=0, M=100, save_path='./dataset.pt')
# dataset: (M, 50), world: World object
```

## Architecture

### State Variables
- `X(t)`: mRNA expression (50 genes)
- `P(t)`: Protein levels
- `Z(t)`: Chromatin state [0, 1]
- `N(t)`: Cell population count

### Update Order (per timestep)
```
normalize_protein → TFinput → Z → X → P_raw → projection → N
```

### Gene Categories (50 total)
| Category | Indices |
|----------|---------|
| TF | 0-5 |
| RBP | 6-10 |
| Kinase | 11-13 |
| Phosphatase | 14-16 |
| Epigenetic | 17-19 |
| Cell Cycle | 20-22 |
| Apoptosis | 23-25 |
| Background | 26-49 |

### Perturbation vs Intervention
- **Perturbation**: parameter-level (e.g., `rho=0` for knockout), persistent
- **Intervention**: state-level (e.g., `X_i=0`), single-step

## Output Format

Saved `.pt` files:
```python
{
    'X_traj': Tensor,  # (T+1, G)
    'P_traj': Tensor,  # (T+1, G)
    'Z_traj': Tensor,  # (T+1, G)
    'N_traj': Tensor,  # (T+1,)
    'world': Dict       # Serialized World
}
```

## Documentation

- `docs/Phase0/`: Phase 0 specifications (Architecture, Design, System Definition, Spec v1.2)
- `docs/mics/`: Parameter reference
- `ENV_SETUP.md`: Environment setup notes (WSL/mamba/pip)
