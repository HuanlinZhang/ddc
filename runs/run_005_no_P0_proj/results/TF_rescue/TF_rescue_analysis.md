# Task 3: TF Support / Rescue Analysis

## 1. Experimental Setup

### 1.1 Objective

Test whether increased TF (Transcription Factor) support can rescue collapsed worlds from collapse to steady-state regime, using **system-internal variables** as required by the analysis plan.

### 1.2 Implementation Constraints (per Analysis Plan)

- All system update rules must remain unchanged
- Resource projection must be preserved (Σ P_i ≤ R_total)
- Interventions must be via system-internal variables (mRNA or parameters)
- **Method A**: Modify X_TF_initial (mRNA), not P (protein)
- **Method B**: Apply clamp to P_raw, then re-apply resource projection

### 1.3 Selected Collapse Seeds

Four collapse seeds were selected based on TF_retention_ratio distribution (see task 2):

| Seed | TF_retention_ratio |
| ---- | ------------------ |
| 546  | 0.001491           |
| 4578 | 0.010444           |
| 3492 | 0.013712           |
| 3441 | 0.155999           |

### 1.4 Intervention Methods

**Method A: TF Overexpression (via mRNA)**

- Initial TF mRNA levels multiplied by ζ: `X_TF_initial ← ζ × X_TF_initial`
- ζ values tested: 2, 3, 4, 5
- This tests whether increased TF transcription can rescue collapse

**Method B: TF Clamp**

- Apply clamp to P_raw before resource projection
- TF protein levels clamped to minimum threshold in early time window (t=0 to t=50)
- Thresholds tested: 0.02, 0.04, 0.06, 0.08
  > Threshold selection rationale: By observing P_traj of steady worlds, the stable protein levels are approximately in the range of 0.02~0.05. Setting threshold too high would not be biologically reasonable.

## 2. Results

### 2.1 Rescue Summary

All 32 experiments across 4 collapse seeds showed **0% rescue success**.

| Method                  | Experiments | Rescued | Success Rate |
| ----------------------- | ----------- | ------- | ------------ |
| TF Overexpression (ζ)  | 16          | 0       | 0%           |
| TF Clamp (threshold)   | 16          | 0       | 0%           |
| **Total**              | **32**      | **0**   | **0%**       |

### 2.2 Per-Seed Results

| Seed | Overexpression (ζ=2~5) | Clamp (thresh=0.02~0.08) |
| ---- | ---------------------- | ------------------------- |
| 546  | 0/4 rescued            | 0/4 rescued              |
| 4578 | 0/4 rescued           | 0/4 rescued              |
| 3492 | 0/4 rescued           | 0/4 rescued              |
| 3441 | 0/4 rescued           | 0/4 rescued              |

## 3. Interpretation

### 3.1 Collapse Is Robust and Self-Reinforcing

The complete absence of rescue success indicates that:

1. **Collapse is not due to insufficient TF levels at initialization**: Even when TF mRNA is artificially overexpressed (up to 5×), the system cannot maintain steady-state expression
2. **Collapse is not due to early-stage TF depletion**: Even when TF protein levels are clamped to minimum thresholds during early development (t=0-50), collapse persists
3. **Collapse is driven by insufficient regulatory gain**: The TF regulatory network cannot sustain expression through its own feedback mechanisms

### 3.2 Mechanistic Implications

These results support the hypothesis that collapse is determined by **intrinsic regulatory dynamics**, not by:

- Initial expression levels (tested in Task 4)
- Resource constraints (tested in Task 1)
- Early TF support (tested in Task 3)

The system appears to have a **low-gain attractor** from which it cannot escape through internal variables alone.

## 4. Output Files

| File | Description |
| ---- | ----------- |
| `results/TF_rescue/rescue_summary.tsv` | Rescue experiment summary (32 rows) |
| `data/trajectories/TF_rescue_trajectories.tsv` | Full trajectory data |
| `plots/TF_rescue/seed{seed}_TF_rescue_trajectories.png` | Per-seed visualization (7 files) |
