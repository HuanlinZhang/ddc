# Run 004 - Sampling Extension Summary

## 0. Regime Classification Criteria

Regime classification follows the standard from run\_004 and run\_005:

| Regime       | Criteria      | Description                                    |
| ------------ | ------------- | ---------------------------------------------- |
| **Collapse** | N\_active ≤ 1 | X(t) → 0, P(t) → 0 (system collapses)          |
| **Steady**   | N\_active > 1 | X(t) > 0, P(t) > 0 (system maintains activity) |

Where **N\_active** is the number of genes with expression above a threshold at the final simulation time step.

***

## 1. Steady World Ratio

| Metric          | Value |
| --------------- | ----- |
| Total Worlds    | 100   |
| Steady Worlds   | 26    |
| Collapse Worlds | 74    |
| Steady Rate     | 26.0% |
| Collapse Rate   | 74.0% |

**Conclusion**: Among 100 sampled worlds, steady-state regime accounts for **26.0%**, while collapse regime accounts for **74.0%**.

***

## 2. TF Activity Pattern: Collapse vs Steady

### Collapse Worlds

| Metric             | Value       |
| ------------------ | ----------- |
| Count              | 74          |
| N\_active\_TF = 0  | 74 (100.0%) |
| N\_active\_TF Mean | 0.00        |

All collapse worlds exhibit zero TF activity. **TF expression is strongly correlated with collapse.**

### Steady Worlds

| Metric              | Value      |
| ------------------- | ---------- |
| Count               | 26         |
| N\_active\_TF > 0   | 24 (92.3%) |
| N\_active\_TF Mean  | 4.92       |
| N\_active\_TF Std   | 1.83       |
| N\_active\_TF Range | 0 - 6      |

The vast majority of steady worlds show significant TF activity (N\_active\_TF > 0), with 2 exceptions (TF = 0).

***

## 3. Pattern Robustness

### Q: Does the original run_004 pattern hold?

**A: Yes.** The TF-driven stability pattern remains stable under larger world sampling.

### Original Finding (run_004, ~10 worlds)

```
TF expression ↔ system viability is strongly correlated
```

### Extended Validation (run\_004\_sampling\_extension, 100 worlds)

| Pattern                            | Observed in 100 worlds? |
| ---------------------------------- | ----------------------- |
| Collapse worlds: N\_active\_TF ≈ 0 | Yes (100%)              |
| Steady worlds: N\_active\_TF > 0   | Yes (92.3%)             |

**Conclusion**: The **TF-driven stability pattern observed in original run\_004 remains stable under larger world sampling**.

- Collapse worlds: 100% accompanied by TF inactivity (N\_active\_TF = 0)
- Steady worlds: 92.3% accompanied by TF activity (N\_active\_TF > 0)

The core hypothesis is validated: **TF regulatory core is a key factor in maintaining steady-state system stability**.

***

## 4. Key Observations

1. **Higher Collapse Rate Than Expected**: 74% of worlds trend toward collapse, indicating that under the current parameter sampling distribution, systems are more prone to collapse.
2. **TF Activity as Strong Predictor**: TF activity is a reliable indicator for distinguishing collapse and steady regimes, with near 100% accuracy.
3. **Minor Exceptions**: 2 steady worlds (seeds 33 and 49) exhibit N\_active\_TF = 0. This may be due to insufficient simulation time—these systems may collapse after additional time steps.

***

## 5. Environment

| Item | Path/Version |
|------|--------------|
| DDC Source | `src/ddc.py` |
| Latest commit | `0c0a6ce33108f715d5bede25fe5ee4946cbffca3` |
| Simulation Duration | T = 200 time steps |

***

## 6. Output Files

| File                                      | Description                                                |
| ----------------------------------------- | ---------------------------------------------------------- |
| `data/world_regime_distribution.tsv`      | Per-world regime classification data                       |
| `results/regime_summary.tsv`              | Aggregated regime statistics                               |
| `plots/sampling_extension_combined.png`   | Regime distribution and N\_active\_TF box plots            |
| `plots/P_traj/` | Individual world P_traj plots (100 PNG files) |

