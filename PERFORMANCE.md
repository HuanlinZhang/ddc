# Performance & Optimization Log

## Benchmark Environment

**Test date**: 2026-05-16

**Machine**: Dell G15-5511 (Windows 11, WSL2)
**CPU**: Intel Core i7-11800H @ 2.30GHz (8 cores, 16 threads)
**RAM**: 64.0 GB (63.7 GB usable)
**GPU**: NVIDIA GeForce RTX 3060 Laptop GPU, 6 GB GDDR6, compute 8.6, 30 multiprocessors
**PyTorch**: 2.4.1+cu121
**Python**: 3.10

## v1.2 Benchmark Results

Single-cell simulation (M=1, CPU, T varying):

| T (timesteps) | Time | Throughput |
|----------------|------|------------|
| 100 | 0.023s | — |
| 200 | 0.037s | — |
| 500 | 0.084s | — |
| 1000 | 0.174s | — |

### Batch Simulation: CPU vs GPU (T=400)

| M | CPU Time | CPU Speed | GPU Time | GPU Speed | GPU Speedup |
|---|---------|-----------|---------|-----------|-------------|
| 100 | 0.59s | 170 cells/sec | 1.05s | 95 cells/sec | 0.6x (GPU slower) |
| 1,000 | 4.95s | 202 cells/sec | 1.85s | 540 cells/sec | 2.7x |
| 10,000 | 54.56s | 183 cells/sec | 16.46s | 608 cells/sec | **3.3x** |
| 20,000 | 106.7s | 187 cells/sec | 32.6s | 614 cells/sec | **3.3x** |
| 50,000 | 267.4s | 187 cells/sec | 112.2s | 446 cells/sec | **2.4x** |
| 100,000 | — | — | OOM killed | — | — |

### Key Findings
- **M < ~200**: CPU faster (GPU overhead exceeds compute time)
- **M ≥ 1,000**: GPU 2.4-3.3x faster
- **M=10,000-20,000**: GPU 3.3x speedup (peak, RTX 3060 VRAM ~6.4GB)
- **M=50,000**: GPU speedup drops to 2.4x (GPU RAM exhausted, spills to system RAM)
- **M=100,000**: OOM killed on this machine (16GB RAM insufficient)
- CPU throughput is constant ~185-200 cells/sec regardless of M (bandwidth-bound)

### Multiprocessing CPU (M=1000, T=400, for comparison):

| n_jobs | Time | Speed |
|---------|------|-------|
| 1 | 33s | 30 cells/sec |
| 4 | 13s | 76 cells/sec |

**Recommendation**: Use batch simulation (`M>1`) over multiprocessing for M≥100. Batch CPU (202 cells/sec) is 2.7x faster than multiprocessing n_jobs=4 (76 cells/sec).

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
