# Run 007 — Summary

## 1. Experiment Overview

### Objective

Analyze how gene-specific knockout (KO) induces regime transitions in steady DDC systems, and quantify system fragility.

### Design

- **Perturbation**: Knockout (rho_i = 0) at the transcription level via parameter perturbation
- **Target Genes**: 6 TF genes (0-5) + 3 epigenetic modifier genes (17-19) = 9 genes total
- **Seeds**: 4 steady worlds from run_005 (seeds: 2026, 4094, 4478, 7709)
- **Simulation**: T=400, terminal window=20

## 2. Key Findings

### 2.1 TF Genes Show High Collapse Susceptibility

- **Genes 0-3**: 100% collapse rate (4/4 seeds) — most fragile
- **Genes 4-5**: 75% collapse rate (3/4 seeds)
- **All TF genes**: average collapse rate = 91.7%
- **Deviation score**: ~1.06-1.23, indicating severe protein level disruption

### 2.2 Epigenetic Genes Are Robust

- **All EPI genes**: 0% collapse rate — maintain steady regime under KO
- **Deviation score**: ~0.0002-0.03, nearly unchanged from baseline
- This suggests EPI genes are not essential regulators in the steady-state GRN

### 2.3 Regime Transition Matrix

| Gene | steady→collapse-like | steady→steady | Total |
|------|---------------------|---------------|-------|
| 0-3 (TF) | 4 | 0 | 4 |
| 4-5 (TF) | 3 | 1 | 4 |
| 17-19 (EPI) | 0 | 4 | 4 |

### 2.4 Class-Level Comparison

| Gene Type | Avg Collapse Rate | Avg ΔP_terminal | Avg Deviation |
|-----------|-------------------|-----------------|---------------|
| TF | 91.7% | -0.018 | 1.06 |
| EPI | 0% | ~0 | 0.10 |

### 2.5 Cross-Run Comparison with run_006

| Gene Type | run_006 Recovery | run_007 Collapse |
|-----------|------------------|-----------------|
| TF | High (1.0) | High (91.7%) |
| EPI | N/A | None (0%) |

**Key Insight**: Genes showing high recovery after transient intervention (run_006) may still cause regime collapse under persistent KO (run_007). This reveals that:
1. Transient vs persistent perturbations yield fundamentally different system responses
2. TF genes are both recoverable under transient perturbation and critical under persistent perturbation
3. EPI genes play a minimal role in maintaining steady-state dynamics

## 3. Output Files

### 3.1 Data Files

#### Trajectories (`data/trajectories/`)

```
per seed (2026, 4094, 4478, 7709):
  ├── baseline_X.tsv    (mRNA trajectory)
  ├── baseline_P.tsv    (protein trajectory)
  ├── baseline_Z.tsv    (protein half-life trajectory)
  ├── KO_gene_i_X.tsv   (KO mRNA trajectory, i=0,1,2,3,4,5,17,18,19)
  ├── KO_gene_i_P.tsv   (KO protein trajectory)
  └── KO_gene_i_Z.tsv   (KO protein half-life trajectory)
```

#### Response Metrics (`data/response_metrics/`)

```
response_summary.tsv
  Columns: seed, gene_id, gene_type, terminal_mean_P_baseline,
           terminal_mean_P_ko, Δmean_P_terminal, deviation_score,
           baseline_regime, ko_regime, transition_type
```

### 3.2 Results Files

#### Per-Seed Results (`results/per_seed/`)

```
seed_2026_results.tsv
seed_4094_results.tsv
seed_4478_results.tsv
seed_7709_results.tsv

Columns: seed_id, gene_id, gene_type, baseline_regime, ko_regime,
         transition_type, Δmean_P_terminal, deviation_score
```

#### Aggregated Results (`results/aggregated/`)

```
gene_level_summary.tsv
  Columns: gene_id, gene_type, collapse_count, collapse_rate,
           mean_Δmean_P_terminal, mean_deviation_score

class_level_summary.tsv
  Columns: gene_type, average_collapse_rate, average_Δmean_P_terminal,
           average_deviation_score

cross_run_comparison.tsv
  Columns: gene_id, gene_type, recovery_score(run_006),
           collapse_rate(run_007), mean_Δmean_P_terminal(run_007)

regime_transition_matrix.tsv
  Rows: gene_id
  Columns: steady→steady, steady→collapse-like
```

### 3.3 Visualization Files

#### Trajectory Plots (`plots/trajectories/`)

```
traj_comparison_seed_2026.png
traj_comparison_seed_4094.png
traj_comparison_seed_4478.png
traj_comparison_seed_7709.png

Layout: 4×3 subplots
  Row 0: Baseline, run_006 reference, empty
  Rows 1-3: 9 KO conditions per seed
```

#### Summary Plots (`plots/`)

```
class_comparison/          - TF vs EPI comparison plots
regime_transition/
  transition_summary.png   - Regime transition bar plot
global_summary.png         - mean_P(t) baseline vs KO comparison
delta_P_over_time.png      - Δmean_P(t) auxiliary plot
```

## 4. Metrics Definitions

| Metric | Definition |
|--------|------------|
| deviation_score | L1 norm of (P_ko - P_baseline), averaged over time and genes |
| Δmean_P_terminal | terminal_mean_P_ko - terminal_mean_P_baseline |
| collapse_rate | Proportion of KO experiments resulting in steady→collapse-like transition |
| transition_type | steady→steady or steady→collapse-like |
| terminal_window | Last 20 time points (t=380-399) |
