# CLAUDE.md

DDC (Designed Digital Cell) - Gene regulatory network simulator.

## Project Structure

```
src/ddc/
├── __init__.py    # Public API: run_simulation, generate_dataset, sample_world, etc.
├── core.py        # Core simulation logic (World, simulate_single_cell, update functions)
├── cli.py         # CLI entry point (argparse)
pyproject.toml
docs/Phase0/       # Phase 0 specifications
docs/mics/         # Parameter reference
```

## Key Functions

| Function | Location | Purpose |
|----------|----------|---------|
| `sample_world(seed)` | core.py | Generate random gene network world |
| `simulate_single_cell(world, X0, P0, Z0, N0, t_steps)` | core.py | Simulate single cell trajectory |
| `run_simulation(seed, save_path)` | core.py | Full pipeline, returns trajectory dict |
| `generate_dataset(world_seed, M, save_path)` | core.py | Multi-cell dataset (M, 50) |
| `apply_perturbation(world, state, config)` | core.py | Parameter-level perturbation |
| `apply_intervention(state, config)` | core.py | State-level intervention |

## Running

**CLI:**
```bash
ddc run --seed 42 --output ./traj.pt
ddc dataset --world-seed 0 --M 100 --output dataset.pt
ddc smoke-test --seed 42
```

**Python:**
```python
from ddc import run_simulation, generate_dataset
traj = run_simulation(seed=42)
dataset, world = generate_dataset(world_seed=0, M=100)
```

**Editable install:**
```bash
pip install -e .
```

## Architecture Notes

- `src/ddc.py` is a backward-compatible shim; all logic is in `src/ddc/core.py`
- `src/ddc/__init__.py` exposes the public API
- Simulation uses `torch.float64` throughout
- `ENABLE_RESOURCE_PROJECTION` global flag in core.py controls resource constraint

## Data Flow

```
normalize_protein → compute_TFinput → update_chromatin → update_mRNA
    → update_protein_raw → apply_resource_projection → update_fate
```

## Perturbation vs Intervention

- **Perturbation** (`apply_perturbation`): modifies World parameters (persistent, changes dynamics)
  - KO: `rho_i = 0` via `config={'knockout': [i]}`
- **Intervention** (`apply_intervention`): modifies State (single-step, does not change dynamics)
  - Knockdown: `config={'knockdown_X': [i]}` sets X_i = 0 at intervention time only
