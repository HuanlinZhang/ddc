# Task 4: Initial Condition Scaling Analysis

## 1. Objective

Test whether collapse depends on initial expression levels, and whether the system exhibits multiple attractors.

## 2. Methods

Initial conditions were scaled by factor γ:

```
X0 → γ × X0
P0 → γ × P0
```

γ values tested: \[0.5, 1.0, 2.0, 5.0]

## 3. Results

### 3.1 Scaling Summary

| Seed | γ   | Original Regime | Scaled Regime | Mean X (terminal) | Mean P (terminal) |
| ---- | --- | --------------- | ------------- | ----------------- | ----------------- |
| 546  | 0.5 | collapse        | collapse      | 0.000000          | 0.000000          |
| 546  | 1.0 | collapse        | collapse      | 0.000000          | 0.000000          |
| 546  | 2.0 | collapse        | collapse      | 0.000000          | 0.000000          |
| 546  | 5.0 | collapse        | collapse      | 0.000000          | 0.000000          |
| 6798 | 0.5 | collapse        | collapse      | 0.000000          | 0.000003          |
| 6798 | 1.0 | collapse        | collapse      | 0.000000          | 0.000003          |
| 6798 | 2.0 | collapse        | collapse      | 0.000000          | 0.000004          |
| 6798 | 5.0 | collapse        | collapse      | 0.000000          | 0.000004          |
| 4578 | 0.5 | collapse        | collapse      | 0.000003          | 0.000039          |
| 4578 | 1.0 | collapse        | collapse      | 0.000003          | 0.000041          |
| 4578 | 2.0 | collapse        | collapse      | 0.000004          | 0.000044          |
| 4578 | 5.0 | collapse        | collapse      | 0.000004          | 0.000047          |
| 2177 | 0.5 | collapse        | collapse      | 0.000000          | 0.000001          |
| 2177 | 1.0 | collapse        | collapse      | 0.000000          | 0.000001          |
| 2177 | 2.0 | collapse        | collapse      | 0.000000          | 0.000002          |
| 2177 | 5.0 | collapse        | collapse      | 0.000000          | 0.000002          |
| 3492 | 0.5 | collapse        | collapse      | 0.000000          | 0.000004          |
| 3492 | 1.0 | collapse        | collapse      | 0.000000          | 0.000004          |
| 3492 | 2.0 | collapse        | collapse      | 0.000000          | 0.000005          |
| 3492 | 5.0 | collapse        | collapse      | 0.000000          | 0.000007          |
| 7483 | 0.5 | collapse        | collapse      | 0.000001          | 0.000008          |
| 7483 | 1.0 | collapse        | collapse      | 0.000001          | 0.000007          |
| 7483 | 2.0 | collapse        | collapse      | 0.000001          | 0.000007          |
| 7483 | 5.0 | collapse        | collapse      | 0.000001          | 0.000006          |
| 3441 | 0.5 | collapse        | collapse      | 0.000004          | 0.001142          |
| 3441 | 1.0 | collapse        | collapse      | 0.000005          | 0.001274          |
| 3441 | 2.0 | collapse        | collapse      | 0.000006          | 0.001403          |
| 3441 | 5.0 | collapse        | collapse      | 0.000007          | 0.001524          |

### 3.2 Regime Transitions

No regime transitions were observed. All scaled initial conditions remained in their original regime.

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

