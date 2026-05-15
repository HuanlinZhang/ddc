# Changelog

All notable changes to DDC are documented here.

## [1.0.0] — 2026-05-16

### Added
- Python package structure (`src/ddc/`) with installable layout
- CLI entry point (`ddc run/dataset/world/smoke-test/sanity`)
- `sample_world(seed)` — random gene regulatory network generation
- `simulate_single_cell(world, X0, P0, Z0, N0, t_steps)` — single-cell trajectory
- `run_simulation(seed, save_path)` — full simulation pipeline
- `generate_dataset(world_seed, M, save_path)` — multi-cell dataset (scRNA-like)
- `apply_perturbation(world, state, config)` — parameter-level perturbation (KO)
- `apply_intervention(state, config)` — state-level intervention (knockdown)
- Gene category annotations (TF, RBP, Kinase, Phosphatase, EPI, CellCycle, Apoptosis, Background)
- `tests/` with pytest-based unit tests
- MIT License
- Phase 0 specifications in `docs/Phase0/`

### Changed
- Refactored from single `ddc.py` into structured package
- CLI replaces direct script execution
- K parameter now sampled as `U(0.1, 1.0) / G` (scaled by gene count)

### Fixed
- Perturbation vs Intervention semantics corrected (KO = parameter-level, knockdown = state-level)
- Intervention records post-intervention state at correct trajectory index

### Architecture
- 50 genes: TF(0-5), RBP(6-10), Kinase(11-13), Phosphatase(14-16), EPI(17-19), CellCycle(20-22), Apoptosis(23-25), Background(26-49)
- Update order: normalize_protein → TFinput → Z → X → P_raw → projection → N
- Resource constraint: ΣP ≤ R_total (toggleable via `ENABLE_RESOURCE_PROJECTION`)
