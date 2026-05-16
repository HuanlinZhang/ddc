# Performance & Optimization Log

## v1.2 Benchmark Results

Single-cell simulation, CPU single-thread:

| T (timesteps) | Time | Throughput |
|----------------|------|------------|
| 100 | 0.023s | — |
| 200 | 0.037s | — |
| 500 | 0.084s | — |
| 1000 | 0.174s | — |

Batch simulation (CPU, M cells, T=400):

| M | Time | Speed |
|---|------|-------|
| 100 | 0.53s | 189 cells/sec |
| 1000 | 5.2s | 192 cells/sec |

Batch simulation (GPU, RTX 3060 Laptop, M cells, T=400):

| M | Time | Speed |
|---|------|-------|
| 100 | 0.38s | 263 cells/sec |
| 1000 | 2.9s | 345 cells/sec |

**Note**: At current gene count (G=50) and batch size (M≤1000), GPU speedup over CPU is modest (~1.5-1.8x). GPU acceleration becomes significant with larger G or M.

Multiprocessing CPU (M=1000, T=400):

| n_jobs | Time | Speed |
|---------|------|-------|
| 1 | 33s | 30 cells/sec |
| 4 | 13s | 76 cells/sec |

**Recommendation**: Use batch simulation (`M>1`) over multiprocessing for M≥100. Batch CPU is 2.5x faster than multiprocessing n_jobs=4.

## Optimization History

### v1.2 — GPU Batch Simulation (2026-05-16)
- `simulate_batch()`: vectorized batch simulation accepting `(M, G)` inputs, returns `(M, T+1, G)` trajectories
- `run_simulation()`: add `M` param (default 1) and `device` param (`'auto'|'cuda'|'cpu'`)
- Auto-detects CUDA, gracefully falls back to CPU if unavailable
- CPU vs GPU produce bit-identical results (max_diff=0)
- CLI: add `--M` and `--device` flags

### v1.1.0 — Vectorized Updates (2026-05-16)
- `compute_TFinput`: Python for-loop → tensor ops (`**` + `masked_fill` + `prod(dim=1)`)
- `update_chromatin`: Python for-loop → matrix multiply (`beta_matrix @ tilde_P`)
- Added `P_mask`, `P_degree`, `a_ij_matrix`, `beta_matrix` tensors to `World`
- **Result**: ~12-14x speedup on single-cell simulation

### v1.0.0 → v1.1.0
- Single-cell T=200: 0.47s → 0.037s (**13x faster**)
- Added multiprocessing via `--n-jobs` (up to 2.6x on 4 cores)

---

## Batch Simulation — Design

### API

```python
# Single-cell (backward compatible, M=1, CPU)
traj = run_simulation(world_seed=42, T=200)
# X_traj: (T+1, G)

# Batch (CPU or GPU auto-detected)
traj = run_simulation(world_seed=42, T=200, M=100)
# X_traj: (M, T+1, G)

# Explicit device
traj = run_simulation(world_seed=42, T=200, M=100, device='cuda')   # GPU
traj = run_simulation(world_seed=42, T=200, M=100, device='cpu')    # CPU
traj = run_simulation(world_seed=42, T=200, M=100, device='auto')  # auto-detect (default)
```

### CLI

```bash
# Auto-detect GPU (default)
ddc run --seed 42 --M 1000 --output traj.pt

# Force CPU
ddc run --seed 42 --M 1000 --device cpu --output traj.pt

# Force GPU
ddc run --seed 42 --M 1000 --device cuda --output traj.pt
```

### Design Principles
1. **Auto-detection by default**: `device='auto'` uses CUDA if available, else CPU
2. **Backward compatible**: `M=1` (default) uses scalar path, result shapes unchanged
3. **Identical results**: CPU and GPU produce bit-identical trajectories (max_diff=0)
4. **Graceful fallback**: if CUDA requested but unavailable, falls back to CPU silently

## Next Steps (Ideas)

### High Priority
1. **snakemake pipeline**: integrate ddc into a reproducible analysis workflow in `ddc-analysis` repo.

2. **Larger G or M for GPU saturation**: current G=50 is too small for GPU to shine. GPU acceleration becomes significant at G≥500 or M≥10000.

### Medium Priority
3. **World graph stats**: expose network statistics (in-degree distribution, edge density, etc.) for analysis.

4. **Checkpointing**: save intermediate trajectory snapshots for long T simulations.

5. **Single-cell GPU**: accelerate M=1 path on GPU (currently CPU only).

### Low Priority / Experimental
6. **JIT compilation**: `torch.compile()` on the update functions for additional speedup.

7. **Mixed precision**: `float32` for intermediate computations, `float64` only for final output — faster on GPU.

8. **Adaptive T**: detect convergence and stop early instead of fixed T steps.
