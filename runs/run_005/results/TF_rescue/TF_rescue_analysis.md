# Task 3: TF Support / Rescue Analysis

## 1. Experimental Setup

### 1.1 Objective

Test whether increased TF (Transcription Factor) support can rescue collapsed worlds from collapse to steady-state regime, using **system-internal variables** as required by the analysis plan.

### 1.2 Implementation Constraints (per Analysis Plan)

- All system update rules must remain unchanged
- Resource projection must be preserved (Σ P\_i ≤ R\_total)
- Interventions must be via system-internal variables (mRNA or parameters)
- **Method A**: Modify X\_TF\_initial (mRNA), not P (protein)
- **Method B**: Apply clamp to P\_raw, then re-apply resource projection

### 1.3 Selected Collapse Seeds

Four collapse seeds were selected based on TF\_retention\_ratio distribution (see task 2):

| Seed | TF\_retention\_ratio |
| ---- | -------------------- |
| 546  | 0.001491             |
| 4578 | 0.010444             |
| 3492 | 0.013712             |
| 3441 | 0.155999             |

### 1.4 Intervention Methods

**Method A: TF Overexpression (via mRNA)**

- Initial TF mRNA levels multiplied by ζ: `X_TF_initial ← ζ × X_TF_initial`
- ζ values tested: 2, 3, 4, 5
- This tests whether increased TF transcription can rescue collapse

**Method B: TF Clamp**

- Apply clamp to P\_raw before resource projection
- TF protein levels clamped to minimum threshold in early time window (t=0 to t=50)
- Thresholds tested: 0.02, 0.04, 0.06, 0.08
  > Threshold selection rationale: By observing P_traj of steady worlds, the stable protein levels are approximately in the range of 0.02~0.05. Setting threshold too high would not be biologically reasonable.
  >
  > see `runs/run_004/plots/TF_dynamics.png`

## 2. Results

### 2.1 Rescue Summary

| Seed | Method         | Support Strength | Original Regime | Rescued Regime |
| ---- | -------------- | ---------------- | --------------- | -------------- |
| 546  | overexpression | ζ=2              | collapse        | collapse       |
| 546  | overexpression | ζ=3              | collapse        | collapse       |
| 546  | overexpression | ζ=4              | collapse        | collapse       |
| 546  | overexpression | ζ=5              | collapse        | collapse       |
| 546  | clamp          | thresh=0.02      | collapse        | collapse       |
| 546  | clamp          | thresh=0.04      | collapse        | collapse       |
| 546  | clamp          | thresh=0.06      | collapse        | collapse       |
| 546  | clamp          | thresh=0.08      | collapse        | collapse       |
| 4578 | overexpression | ζ=2              | collapse        | collapse       |
| 4578 | overexpression | ζ=3              | collapse        | collapse       |
| 4578 | overexpression | ζ=4              | collapse        | collapse       |
| 4578 | overexpression | ζ=5              | collapse        | collapse       |
| 4578 | clamp          | thresh=0.02      | collapse        | collapse       |
| 4578 | clamp          | thresh=0.04      | collapse        | collapse       |
| 4578 | clamp          | thresh=0.06      | collapse        | collapse       |
| 4578 | clamp          | thresh=0.08      | collapse        | collapse       |
| 3492 | overexpression | ζ=2              | collapse        | collapse       |
| 3492 | overexpression | ζ=3              | collapse        | collapse       |
| 3492 | overexpression | ζ=4              | collapse        | collapse       |
| 3492 | overexpression | ζ=5              | collapse        | collapse       |
| 3492 | clamp          | thresh=0.02      | collapse        | collapse       |
| 3492 | clamp          | thresh=0.04      | collapse        | collapse       |
| 3492 | clamp          | thresh=0.06      | collapse        | collapse       |
| 3492 | clamp          | thresh=0.08      | collapse        | collapse       |
| 3441 | overexpression | ζ=2              | collapse        | collapse       |
| 3441 | overexpression | ζ=3              | collapse        | collapse       |
| 3441 | overexpression | ζ=4              | collapse        | collapse       |
| 3441 | overexpression | ζ=5              | collapse        | collapse       |
| 3441 | clamp          | thresh=0.02      | collapse        | collapse       |
| 3441 | clamp          | thresh=0.04      | collapse        | collapse       |
| 3441 | clamp          | thresh=0.06      | collapse        | collapse       |
| 3441 | clamp          | thresh=0.08      | collapse        | collapse       |

### 2.2 Summary Statistics

- **Total experiments**: 32 (4 seeds × 2 methods × 4 strength levels)
- **Rescued**: 0 (0%)
- **Not rescued**: 32 (100%)

## 3. Analysis

### 3.1 Key Findings

1. **No Rescue Achieved**: Neither TF overexpression (via mRNA) nor TF clamp interventions succeeded in rescuing any of the four collapse worlds to steady-state regime.
2. **Seed Variation**: The selected seeds span a wide range of TF\_retention\_ratio (0.0015 to 0.156), representing different levels of TF retention in collapsed states. Despite this variation, all interventions failed uniformly.
3. **Method Comparison**: Both intervention methods showed identical results (no rescue), suggesting that:
   - The collapse mechanism may be self-reinforcing and not easily reversed by early TF intervention
   - The resource projection dynamics may dominate over TF levels in determining system fate

### 3.2 Interpretation

The inability to rescue collapsed worlds through TF intervention suggests that **collapse is inherently robust**. Combined with Task 2 results which show that collapse worlds have extremely low TF retention ratio (~4.5% vs ~44% in steady worlds), this indicates that collapse is likely driven by insufficient regulatory gain in the system—where even elevated TF levels cannot effectively sustain downstream gene activation due to inherently weak regulatory connections.

## 4. Output Files

| File                                                    | Description                                                 |
| ------------------------------------------------------- | ----------------------------------------------------------- |
| `TF_rescue/rescue_summary.tsv`                          | Summary table with all experimental results                 |
| `TF_rescue/TF_rescue_trajectories.tsv`                  | Full trajectory data (X, P, Z for all genes and time steps) |
| `plots/TF_rescue/seed{seed}_TF_rescue_trajectories.png` | Per-seed visualization (3×4 grid), 4 files for 4 seeds      |
| `TF_rescue/TF_rescue_analysis.md`                       | This analysis report                                        |

## 5. Conclusions

- **TF support interventions are insufficient to rescue collapsed worlds** in this model, even when implemented correctly via system-internal variables (mRNA modification)
- The collapse phenotype appears **robust across different seed types and intervention strengths**
- Results suggest collapse is driven by **deeper regulatory mechanisms** beyond simple TF expression levels
- Future work may explore alternative rescue strategies or investigate the fundamental limits of TF-based rescue in this system

