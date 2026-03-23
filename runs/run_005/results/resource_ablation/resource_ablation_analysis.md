# Task 1: Resource Projection Ablation Analysis

## 1. Summary

This analysis tests whether the global resource constraint (`Σ P_i ≤ R_total`) influences the viability boundary of Phase0 regulatory worlds.

**Key Finding:** Removing resource projection has **no significant effect** on regime classification. This suggests that collapse vs steady-state behavior is determined primarily by regulatory dynamics, not by the global resource constraint.

***

## 2. Change in Collapse Ratio

| Condition                     | Collapse | Steady | Collapse Ratio |
| ----------------------------- | -------- | ------ | -------------- |
| Original (with projection)    | 7        | 4      | 63.6%          |
| Ablation (projection removed) | 7        | 4      | 63.6%          |
| **Change**                    | 0        | 0      | **0%**         |

### Regime Transition Matrix

| Original \ Ablation | Collapse | Steady |
| ------------------- | -------- | ------ |
| **Collapse**        | 7        | 0      |
| **Steady**          | 0        | 4      |

**Interpretation:** No regime transitions occurred when resource projection was removed. All collapse worlds remained collapsed, and all steady-state worlds remained in steady state.

***

## 3. Effect on TF Activity

### Collapse Worlds (n=7)

| Metric           | Original | Ablation | Change          |
| ---------------- | -------- | -------- | --------------- |
| N\_active\_TF    | 0        | 0        | None            |
| Mean X\_terminal | \~10⁻⁷   | \~10⁻³   | Slight increase |

TF activity remained at zero for all collapse worlds under both conditions. The slight increase in mean X\_terminal after ablation is negligible for biological function.

### Steady-State Worlds (n=4)

| Seed | N\_active\_TF (Original) | N\_active\_TF (Ablation) | Change |
| ---- | ------------------------ | ------------------------ | ------ |
| 2026 | 6                        | 4                        | -2     |
| 4094 | 6                        | 6                        | 0      |
| 4478 | 5                        | 4                        | -1     |
| 7709 | 6                        | 6                        | 0      |

| Metric             | Original | Ablation | Change |
| ------------------ | -------- | -------- | ------ |
| Mean N\_active\_TF | 5.75     | 5.0      | -0.75  |
| Mean X\_terminal   | 1.19     | 0.82     | -31%   |

TF activity decreased slightly in 2 out of 4 steady-state worlds after ablation, but all remained in steady-state regime.

***

## 4. Qualitative Trajectory Differences

### Collapse Worlds

- **Original:** Rapid decay to near-zero expression, TF expression drops quickly
- **Ablation:** Similar rapid decay pattern, slightly higher terminal expression but still collapsed
- **Conclusion:** Resource projection does not prevent collapse

### Steady-State Worlds

- **Original:** Sustained expression levels, stable TF activity
- **Ablation:** Slightly reduced expression levels, minor reduction in TF count, but regime maintained
- **Conclusion:** Resource projection provides modest stabilization but is not essential for viability

***

## 5. Interpretation

### Answer to Core Question

**Does resource projection influence the viability boundary?**

**Answer: No.** The resource projection module has minimal impact on regime classification.

### Case Classification

This result corresponds to **Case B** from the analysis plan:

> **Projection removal shows no significant effect**
>
> → Collapse is determined by regulatory dynamics

### Mechanistic Implications

1. **Resource constraint is not the bottleneck**: The global protein limit (`Σ P_i ≤ R_total`) does not determine whether a regulatory world collapses or survives.
2. **Regulatory dynamics dominate**: Collapse vs steady-state behavior may be primarily driven by:
   - TF network connectivity
   - Regulatory feedback strength
   - Initial condition sensitivity
   - Gene circuit topology
3. **Resource projection provides marginal benefit**: In steady-state worlds, resource projection slightly increases TF activity (by \~13% on average), suggesting a modest stabilizing effect, but this is not sufficient to change regime classification.

***

## 6. Next Steps

Based on these findings, subsequent tasks should focus on:

1. **Task 2**: Identify dynamical metrics that capture effective regulatory gain
2. **Task 3**: Test whether TF support can rescue collapse worlds
3. **Task 4**: Evaluate sensitivity to initial conditions

These directions address the regulatory dynamics that appear to be the primary determinant of viability.

***

## 7. Data Files

- Regime comparison: `results/resource_ablation/regime_comparison.tsv`
- Trajectory data: `data/trajectories/resource_ablation_trajectories.tsv`
- Visualizations: `plots/resource_ablation/seed*_resource_ablation_trajectories.png`
- Summary heatmap: `plots/resource_ablation/resource_ablation_summary.png`

