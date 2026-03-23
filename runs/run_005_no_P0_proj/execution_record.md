# Execution Record

## 1. Analysis Overview

| Item                        | Content                                                   |
| --------------------------- | --------------------------------------------------------- |
| **Analysis Name**           | Phase 0 Regime Analysis — Mechanism Testing               |
| **Run ID**                  | run\_005\_no\_P0\_proj                                    |
| **Author**                  | zhanghl                                                   |
| **Analysis Date**           | 2026-03-23                                                |
| **Script**                  | run\_005\_analysis.py                                     |
| **Simulation Code Version** | v1.1 (no P0 projection); git commit: 53d9fbe (src/ddc.py) |

<br />

***

## 2. Analysis Objective

Test mechanistic hypotheses about the viability boundary of Phase 0 regulatory worlds:

```
What determines whether a regulatory world collapses or maintains steady-state?
```

Four experimental tasks were designed to probe different mechanistic hypotheses through targeted perturbations.

***

## 3. Task Overview

| Task   | Experiment                    | Key Question                                                   |
| ------ | ----------------------------- | -------------------------------------------------------------- |
| Task 1 | Resource Projection Ablation  | Does global resource constraint determine viability?           |
| Task 2 | Effective Gain Proxy Analysis | What dynamical metrics distinguish collapse from steady-state? |
| Task 3 | TF Support / Rescue Test      | Can increased TF support rescue collapsed worlds?              |
| Task 4 | Initial Condition Scaling     | Does viability depend on initial expression levels?            |

***

## 4. Input Data

### 4.1 Source Directories

| Directory                                           | Description                              |
| --------------------------------------------------- | ---------------------------------------- |
| `test_convergence/v1_1_enable_resource_projection`  | Worlds with resource projection enabled  |
| `test_convergence/v1_1_disable_resource_projection` | Worlds with resource projection disabled |

### 4.2 World Seeds

| Total | Collapse Seeds                          | Steady Seeds           |
| ----- | --------------------------------------- | ---------------------- |
| 11    | 546, 2177, 3441, 3492, 4578, 6798, 7483 | 2026, 4094, 4478, 7709 |

***

## 5. Analysis Parameters

### 5.1 Common Parameters

| Parameter         | Value | Description                               |
| ----------------- | ----- | ----------------------------------------- |
| `THRESHOLD`       | 0.1   | Active gene classification threshold      |
| Simulation Length | T     | Full trajectory for regime classification |

### 5.2 Task 3 Parameters (TF Rescue)

| Parameter          | Values                 | Description                               |
| ------------------ | ---------------------- | ----------------------------------------- |
| `zeta_values`      | 2, 3, 4, 5             | Overexpression multiplier for TF mRNA     |
| `clamp_window`     | 50                     | Time window for TF clamping (t=0 to t=50) |
| `clamp_thresholds` | 0.02, 0.04, 0.06, 0.08 | Minimum TF protein level after clamp      |

### 5.3 Task 4 Parameters (Initial Condition Scaling)

| Parameter      | Values             | Description                  |
| -------------- | ------------------ | ---------------------------- |
| `gamma_values` | 0.5, 1.0, 2.0, 5.0 | Scaling factor for X0 and P0 |

***

## 6. Implementation Constraints

### P0 Projection Disabled (v1.1 modification)

The simulation code (v1.1) has been modified to **not perform resource projection on initial protein state P0**:

```python
# In sample_initial_state():
P0_raw: Tensor = world.gamma * X0
# P0: Tensor = apply_resource_projection(P0_raw, world)  # Disabled
P0: Tensor = P0_raw
```

### Task 3 (TF Support) — Critical Rules

All interventions must be implemented via **system-internal variables only**:

- **Allowed**: Modify X\_TF\_initial (mRNA), apply clamp to P\_raw then re-apply resource projection
- **Forbidden**: Direct modification of protein state P, skipping resource projection in simulation loop, creating ΣP > R\_total states

This ensures we test the system's intrinsic capability to support TF expression, not external perturbations.

***

## 7. Output Files

### 7.1 Task 1 — Resource Projection Ablation

| File                                                      | Description             |
| --------------------------------------------------------- | ----------------------- |
| `results/resource_ablation/regime_comparison.tsv`         | Regime comparison table |
| `results/resource_ablation/resource_ablation_analysis.md` | Analysis report         |
| `data/trajectories/resource_ablation_trajectories.tsv`    | Full trajectory data    |

### 7.2 Task 2 — Effective Gain Proxy Analysis

| File                                                           | Description         |
| -------------------------------------------------------------- | ------------------- |
| `results/gain_proxy_analysis/proxy_metrics.tsv`                | Proxy metrics table |
| `results/gain_proxy_analysis/viability_proxy_summary.tsv`      | Summary statistics  |
| `results/gain_proxy_analysis/effective_gain_proxy_analysis.md` | Analysis report     |

### 7.3 Task 3 — TF Rescue

| File                                           | Description               |
| ---------------------------------------------- | ------------------------- |
| `results/TF_rescue/rescue_summary.tsv`         | Rescue experiment summary |
| `results/TF_rescue/TF_rescue_trajectories.tsv` | Full trajectory data      |
| `results/TF_rescue/TF_rescue_analysis.md`      | Analysis report           |

### 7.4 Task 4 — Initial Condition Scaling

| File                                                                  | Description                      |
| --------------------------------------------------------------------- | -------------------------------- |
| `results/initial_condition_scaling/scaling_summary.tsv`               | Scaling experiment summary       |
| `data/trajectories/scaling_trajectories.tsv`                          | Full trajectory data             |
| `plots/initial_condition_scaling/seed{seed}_scaling_trajectories.png` | Per-seed visualization (7 files) |

### 7.5 Final Synthesis

| File                                            | Description               |
| ----------------------------------------------- | ------------------------- |
| `results/Phase0_viability_mechanism_summary.md` | Overall mechanism summary |

***

## 8. Key Findings Summary

| Finding                                          | Evidence                                              |
| ------------------------------------------------ | ----------------------------------------------------- |
| Resource projection does not determine viability | No regime transitions after ablation                  |
| TF decay drives collapse                         | TF Retention Ratio: 4.5% (collapse) vs 44% (steady)   |
| TF support cannot rescue collapse                | 32/32 experiments showed 0% rescue                    |
| Collapse is initial-condition independent        | No regime transitions across γ ∈ {0.5, 1.0, 2.0, 5.0} |

***

## 9. Reproduction

To reproduce this analysis:

```bash
python scripts/run_005_analysis.py
```

Requires:

- Source data in `test_convergence/v1_1_enable_resource_projection` and `test_convergence/v1_1_disable_resource_projection`
- DDC module at `/home/zhanghl/projects/ddc_github/src`

