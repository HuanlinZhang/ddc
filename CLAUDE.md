# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DDC (Designed Digital Cell) is a computational framework for simulating gene regulatory network (GRN) dynamics at the single-cell level. It supports Monte Carlo sampling, multi-cell dataset generation, and gene perturbation experiments.

## Environment Setup

```bash
pip install torch numpy
# or
pip install -e .
```

Requirements: Python 3.8+, PyTorch 1.9.0+, NumPy 1.20.0+

## Running

**As module:**
```python
from ddc import run_simulation
result = run_simulation(seed=42, save_path='./data/traj.pt')
```

**Direct execution:**
```bash
python src/ddc.py  # runs smoke test with T=10
```

**Run tests:**
```bash
python test_stages.py        # Stage A/B/C tests
python convergence_test.py   # Convergence testing
```

## Architecture

### Core State Variables
- `X` (G,): mRNA expression levels
- `P` (G,): Protein levels
- `Z` (G,): Chromatin states (0-1)
- `N`: Cell count

### Gene Categories (50 genes total)
| Category | Indices | Count |
|----------|---------|-------|
| TF (Transcription Factors) | 0-5 | 6 |
| RBP | 6-10 | 5 |
| Kinase | 11-13 | 3 |
| Phosphatase | 14-16 | 3 |
| Epigenetic | 17-19 | 3 |
| Cell Cycle | 20-22 | 3 |
| Apoptosis | 23-25 | 3 |
| Background | 26-49 | 24 |

### Key Functions (`src/ddc.py`)
| Function | Purpose |
|----------|---------|
| `sample_world(seed)` | Generate random gene network world |
| `simulate_single_cell(world, X0, P0, Z0, N0, t_steps)` | Simulate single cell trajectory |
| `generate_dataset(world_seed, M, save_path)` | Generate multi-cell dataset (M cells) |
| `run_simulation(seed, save_path)` | Full simulation pipeline |
| `apply_perturbation(world, state, config)` | Parameter-level perturbation |
| `apply_intervention(state, config)` | State-level intervention |

### Data Flow (per timestep)
```
X(t), P(t), Z(t), N(t)
    → normalize_protein: P̃ = P / ΣP
    → compute_TFinput: TF = (∏P̃ᵃʲ)^(1/dᵢ)
    → update_chromatin: Z = σ(α + ΣβᵢⱼP̃ⱼ)
    → update_mRNA: X' = (1-δₓ)X + Z·ρ·hill(TF)
    → update_protein: P_raw = (1-δₚ)P + γX
    → resource_projection: if ΣP > R_total, scale proportionally
    → update_fate: N' = N + r·N·(1-N/K_pop)
    → X(t+1), P(t+1), Z(t+1), N(t+1)
```

### Global Flags
`ENABLE_RESOURCE_PROJECTION` (default: True) in `src/ddc.py` - toggle resource constraint

### Output Format
Saved `.pt` files contain:
```python
{
    'X_traj': Tensor,  # (T+1, G)
    'P_traj': Tensor,  # (T+1, G)
    'Z_traj': Tensor,  # (T+1, G)
    'N_traj': Tensor,  # (T+1,)
    'world': Dict      # Serialized World object
}
```
