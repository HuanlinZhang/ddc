# Task 4: Initial Condition Scaling Analysis

## 1. Objective

Test whether collapse depends on initial expression levels, and whether the system exhibits multiple attractors.

## 2. Methods

Initial conditions were scaled by factor γ:

```
X0 → γ × X0
P0 → γ × P0
```

γ values tested: [0.5, 1.0, 2.0, 5.0]

## 3. Results

### 3.1 Scaling Summary

| Seed | γ   | Original Regime | Scaled Regime | Mean X (terminal) | Mean P (terminal) |
| ---- | --- | --------------- | ------------- | ----------------- | ----------------- |
| 546  | 0.5 | collapse        | collapse      | 2.94e-08          | 3.00e-07          |
| 546  | 1.0 | collapse        | collapse      | 3.06e-08          | 3.21e-07          |
| 546  | 2.0 | collapse        | collapse      | 3.17e-08          | 3.46e-07          |
| 546  | 5.0 | collapse        | collapse      | 3.23e-08          | 3.85e-07          |
| 6798 | 0.5 | collapse        | collapse      | 1.48e-08          | 2.81e-06          |
| 6798 | 1.0 | collapse        | collapse      | 7.11e-09          | 3.11e-06          |
| 6798 | 2.0 | collapse        | collapse      | 2.81e-09          | 3.54e-06          |
| 6798 | 5.0 | collapse        | collapse      | 6.91e-10          | 4.41e-06          |
| 4578 | 0.5 | collapse        | collapse      | 3.27e-06          | 3.91e-05          |
| 4578 | 1.0 | collapse        | collapse      | 3.45e-06          | 4.13e-05          |
| 4578 | 2.0 | collapse        | collapse      | 3.67e-06          | 4.38e-05          |
| 4578 | 5.0 | collapse        | collapse      | 3.96e-06          | 4.73e-05          |
| 2177 | 0.5 | collapse        | collapse      | 1.92e-07          | 1.31e-06          |
| 2177 | 1.0 | collapse        | collapse      | 2.08e-07          | 1.43e-06          |
| 2177 | 2.0 | collapse        | collapse      | 2.33e-07          | 1.61e-06          |
| 2177 | 5.0 | collapse        | collapse      | 2.76e-07          | 1.92e-06          |
| 3492 | 0.5 | collapse        | collapse      | 3.39e-07          | 3.67e-06          |
| 3492 | 1.0 | collapse        | collapse      | 3.26e-07          | 4.15e-06          |
| 3492 | 2.0 | collapse        | collapse      | 2.77e-07          | 4.85e-06          |
| 3492 | 5.0 | collapse        | collapse      | 1.41e-07          | 6.94e-06          |
| 7483 | 0.5 | collapse        | collapse      | 1.06e-06          | 7.98e-06          |
| 7483 | 1.0 | collapse        | collapse      | 9.99e-07          | 7.49e-06          |
| 7483 | 2.0 | collapse        | collapse      | 9.10e-07          | 6.79e-06          |
| 7483 | 5.0 | collapse        | collapse      | 7.54e-07          | 5.59e-06          |
| 3441 | 0.5 | collapse        | collapse      | 4.35e-06          | 1.13e-03          |
| 3441 | 1.0 | collapse        | collapse      | 5.12e-06          | 1.27e-03          |
| 3441 | 2.0 | collapse        | collapse      | 5.90e-06          | 1.40e-03          |
| 3441 | 5.0 | collapse        | collapse      | 6.68e-06          | 1.52e-03          |

### 3.2 Regime Transitions

No regime transitions were observed. All scaled initial conditions remained in their original regime.

| γ Value | Regime Transitions |
| ------- | ----------------- |
| 0.5     | 0 / 7             |
| 1.0     | 0 / 7             |
| 2.0     | 0 / 7             |
| 5.0     | 0 / 7             |

## 4. Interpretation

The absence of regime transitions suggests that:

1. **Collapse is not initial-condition dependent**: The collapse regime appears to be a stable attractor that is not influenced by the initial expression levels.
2. **No multiple attractors**: The system does not appear to have multiple attractors based on initial conditions.
3. **Collapse is driven by intrinsic dynamics**: Rather than being triggered by initialization, collapse is likely determined by the underlying regulatory network dynamics and parameters.

## 5. Output Files

| File | Description |
| ---- | ----------- |
| `results/initial_condition_scaling/scaling_summary.tsv` | Summary table with all scaling experiments |
| `data/trajectories/scaling_trajectories.tsv` | Full trajectory data |
| `plots/initial_condition_scaling/seed{seed}_scaling_trajectories.png` | Per-seed visualization (1×5 grid), 7 files |
