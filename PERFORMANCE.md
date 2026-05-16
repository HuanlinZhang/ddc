# Performance & Optimization Log

## v1.1.0 Benchmark Results

Single-cell simulation, CPU single-thread:

| T (timesteps) | Time | Throughput |
|----------------|------|------------|
| 100 | 0.023s | — |
| 200 | 0.037s | — |
| 500 | 0.084s | — |
| 1000 | 0.174s | — |

Multi-cell dataset generation (M=1000, T=400, 8-core server):

| n_jobs | Time | Speed |
|---------|------|-------|
| 1 | 33s | 30 cells/sec |
| 4 | 13s | 76 cells/sec |
| 8 | 14s | 73 cells/sec |
| 12 | 17s | 60 cells/sec |

**Recommendation**: use `n_jobs=4` — 2.6x speedup, beyond that yields no gain or slowdown.

## Optimization History

### v1.1.0 — Vectorized Updates (2026-05-16)
- `compute_TFinput`: Python for-loop → tensor ops (`**` + `masked_fill` + `prod(dim=1)`)
- `update_chromatin`: Python for-loop → matrix multiply (`beta_matrix @ tilde_P`)
- Added `P_mask`, `P_degree`, `a_ij_matrix`, `beta_matrix` tensors to `World`
- **Result**: ~12-14x speedup on single-cell simulation

### v1.0.0 → v1.1.0
- Single-cell T=200: 0.47s → 0.037s (**13x faster**)
- Added multiprocessing via `--n-jobs` (up to 2.6x on 4 cores)

## Next Steps (Ideas)

### High Priority
1. **GPU batch simulation**: modify `simulate_single_cell` to accept batched initial states, output `(M, T+1, G)` tensor — single instruction multiple data on GPU. Potential 10-100x speedup for large M.

2. **snakemake pipeline**: integrate ddc into a reproducible analysis workflow in `ddc-analysis` repo.

### Medium Priority
3. **Batch `simulate_single_cell`**: refactor main loop to accept `X0, P0, Z0` as `(M, G)` tensors — GPU acceleration without changing the scalar simulation function.

4. **World graph stats**: expose network statistics (in-degree distribution, edge density, etc.) for analysis.

5. **Checkpointing**: save intermediate trajectory snapshots for long T simulations.

### Low Priority / Experimental
6. **JIT compilation**: `torch.compile()` on the update functions for additional speedup.

7. **Mixed precision**: `float32` for intermediate computations, `float64` only for final output — faster on GPU.

8. **Adaptive T**: detect convergence and stop early instead of fixed T steps.

---

## GPU Batch Simulation — Design (v1.2)

### Goal
Support batch simulation on GPU while remaining backward-compatible with CPU-only servers.

### Design Principles
1. **Auto-detection**: `torch.cuda.is_available()` decides whether to use GPU
2. **Zero API change**: existing code (`run_simulation(world_seed=42)`) continues to work unchanged
3. **Opt-in batch**: new `M` parameter enables batch mode for multi-cell simulations
4. **Graceful fallback**: if GPU unavailable in batch mode, fallback to CPU multiprocessing

### API Changes

```python
# Existing API — unchanged, works on CPU or GPU
traj = run_simulation(world_seed=42, T=500)
# traj: {'X_traj': (T+1, G), ...}

# New batch API — GPU-accelerated when available
traj = run_simulation(world_seed=42, T=500, M=100)
# traj: {'X_traj': (M, T+1, G), 'P_traj': (M, T+1, G), ...}
```

### Implementation Plan

#### Layer 1: Vectorized per-timestep functions (already done in v1.1.0)
- `compute_TFinput`: works on `(M, G)` batched tensors natively
- `update_chromatin`: `beta_matrix @ tilde_P.T` → `.T` gives `(M, G)`
- All per-step functions are already broadcast-compatible

#### Layer 2: `simulate_batch(world, X0, P0, Z0, N0, t_steps)` → `(M, T+1, G)`
- Refactor main loop to handle `(M, G)` tensors per step
- Each step processes all M cells simultaneously
- GPU tensor operations = M-way parallelism per step

#### Layer 3: `run_simulation()` adds `M` parameter
- `M=1` (default): calls scalar `simulate_single_cell`, result shapes match current API
- `M>1`: calls `simulate_batch`
  - GPU available → GPU batch (fastest)
  - GPU unavailable → CPU multiprocessing (existing path)

#### Backward Compatibility
- `ddc run --seed 42` works exactly as before
- `generate_dataset()` continues to work, may use batch internally in future
