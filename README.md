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
pip install git+https://github.com/HuanlinZhang/ddc.git
```

Requires: Python 3.8+, PyTorch 2.0+, NumPy 1.20+

## CLI

```bash
# Run simulation (all adjustable parameters shown)
ddc run --seed 42 --output ./traj.pt
ddc run --seed 42 --T 500 --output ./traj.pt                    # custom timesteps
ddc run --seed 42 --cell-seed 100 --output ./traj.pt           # separate cell seed
ddc run --seed 42 --knockout-gene 0 1 --output ./traj.pt       # gene knockout (persistent)
ddc run --seed 42 --knockdown-gene 0 --intervention-time 50 --output ./traj.pt  # knockdown (single-step)
ddc run --seed 42 --disable-resource-projection --output ./traj.pt  # disable ΣP ≤ R_total

# Generate multi-cell dataset
ddc dataset --world-seed 0 --M 100 --output dataset.pt

# Sample and save a world
ddc world --seed 42 --output world.json

# Smoke test
ddc smoke-test

# Sanity tests
ddc sanity
```

### `ddc run` — Parameter Reference

| CLI Parameter | Default | Description |
|--------------|---------|-------------|
| `--seed` | *required* | Random seed for gene network world |
| `--output` | *required* | Output `.pt` file path |
| `--T` | 200 | Simulation timesteps |
| `--cell-seed` | `--seed` + 1 | Random seed for cell initial state |
| `--intervention-time` | None | Apply state-level intervention at this timestep |
| `--knockdown-gene` | None | Gene indices for knockdown (X_i=0, single-step) |
| `--knockout-gene` | None | Gene indices for knockout (rho_i=0, persistent) |
| `--disable-resource-projection` | False | Disable ΣP ≤ R_total resource constraint |

## Python API

```python
from ddc import run_simulation, generate_dataset, sample_world

# Single-cell simulation (world_seed=42, cell_seed defaults to 43)
traj = run_simulation(world_seed=42, save_path='./traj.pt')
# traj: X_traj, P_traj, Z_traj, N_traj

# With explicit cell seed
traj = run_simulation(world_seed=42, cell_seed=100, save_path='./traj.pt')

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

### Simulation Parameters (Hardcoded Defaults)

The following parameters are currently hardcoded in `src/ddc/core.py`. They are documented here for reference; adjust them in code or via the Python API.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `G` | 50 | Number of genes |
| `T` | 200 | Default simulation timesteps |
| `R_TOTAL` | 1.0 | Total protein resource capacity |
| `K_POP` | 1.0 | Cell population carrying capacity |
| `EPSILON` | 1e-8 | Numerical stability constant |
| `DTYPE` | float64 | Tensor data type |
| `ENABLE_RESOURCE_PROJECTION` | True | Resource constraint toggle |
| `A_IJ_RANGE` | (0.5, 2.0) | Regulatory interaction strength range |

**World sampling ranges** (set in `sample_world()`):

| Parameter | Range | Description |
|-----------|-------|-------------|
| `alpha` | N(0, 1) | Basal transcription rate |
| `rho` | U(0.5, 2.0) | Transcription efficiency |
| `K` | U(0.1, 1.0) / G | Hill function half-saturation |
| `n` | 2.0 (fixed) | Hill coefficient |
| `delta_x` | U(0.1, 0.5) | mRNA degradation rate |
| `delta_p` | U(0.05, 0.3) | Protein degradation rate |
| `gamma` | 1.0 (fixed) | Translation efficiency |
| `r` | U(0.05, 0.2) | Cell growth rate |
| `a_ij` | U(0.5, 2.0) | TF-mediated regulatory strength |
| `beta_ij` | N(0, 1.5) | Chromatin interaction strength |

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
- `ENV_SETUP.md`: Environment setup notes (WSL/mamba/pip)
- `CHANGELOG.md`: Version history
- `LICENSE`: MIT License

## Testing

```bash
pytest tests/
```

## Quick Start

```python
from ddc import run_simulation, generate_dataset

# Run simulation
traj = run_simulation(seed=42, save_path='./traj.pt')

# Generate dataset
dataset, world = generate_dataset(world_seed=0, M=100)
```
