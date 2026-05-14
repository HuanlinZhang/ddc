# Run 007 — Execution Record

## 1. Seeds Used

4 steady worlds from run\_005:

| Seed ID | Baseline Regime |
| ------- | --------------- |
| 2026    | steady          |
| 4094    | steady          |
| 4478    | steady          |
| 7709    | steady          |

## 2. Gene List

**Target Genes**: All TF genes + all epigenetic modifier genes (9 genes total)

| Gene ID | Gene Type                  |
| ------- | -------------------------- |
| 0-5     | TF (Transcription Factors) |
| 17-19   | Epigenetic Modifiers       |

- TF\_GENES: \[0, 1, 2, 3, 4, 5]
- EPI\_GENES: \[17, 18, 19]

## 3. Simulation Parameters

| Parameter            | Value | Description                                |
| -------------------- | ----- | ------------------------------------------ |
| T                    | 400   | Total simulation time                      |
| TERMINAL\_WINDOW     | 20    | Window size for terminal state calculation |
| N\_ACTIVE\_THRESHOLD | 1     | Threshold for active gene count            |
| THRESHOLD            | 0.1   | Protein concentration threshold            |
| N                    | 1     | Number of simulation runs per condition    |

## 4. Code Version / Commit Hash

- Latest commit: `0c0a6ce33108f715d5bede25fe5ee4946cbffca3`
- Commit message: "Add run\_006: External Perturbation Response Analysis"
- Source code: `src/ddc.py`

## 5. Deviation Metric Definition

**deviation_score**: L1 norm of (P_ko - P_baseline), averaged over time and genes.

### Concept:

```
deviation_score = (1/T) * Σ_t ( Σ_g |P_ko(t,g) - P_baseline(t,g)| )
```

Where:
- T: Number of simulation time steps (400)
- g: Gene index
- P_ko(t,g): Protein concentration of gene g at time t under KO condition
- P_baseline(t,g): Protein concentration of gene g at time t under Baseline condition

### Physical Interpretation:

- Measures the impact of KO perturbation on the entire system's protein concentration
- Averages over all genes and time points within the terminal window
- Larger values indicate more severe disruption of the steady state system by the KO

### Related Metrics:

| Metric | Definition |
|--------|------------|
| terminal_mean_P_baseline | Mean protein concentration under Baseline within terminal window |
| terminal_mean_P_ko | Mean protein concentration under KO within terminal window |
| delta_mean_P_terminal | terminal_mean_P_ko - terminal_mean_P_baseline |

## 6. Terminal Window Definition

- **Window Size**: 20 time points
- **Definition**: The last 20 time points of the simulation (t = 380 to 399)
- **Purpose**: Used for calculating terminal state metrics (terminal_mean_P, delta_mean_P_terminal)

