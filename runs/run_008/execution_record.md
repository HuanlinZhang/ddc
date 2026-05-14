# Run 008 — Execution Record

## 1. Seeds Used

4 steady worlds from run\_007 (which originally selected from run\_005):

| Seed ID | Baseline Regime |
| ------- | --------------- |
| 2026    | steady          |
| 4094    | steady          |
| 4478    | steady          |
| 7709    | steady          |

## 2. Gene List

**Target Genes**: All TF genes (6 genes total)

| Gene ID | Gene Type                  |
| ------- | -------------------------- |
| 0-5     | TF (Transcription Factors) |

- TF\_GENES: \[0, 1, 2, 3, 4, 5]
- EPI\_GENES: \[17, 18, 19] (not perturbed in this run)

## 3. Simulation Parameters

| Parameter        | Value       | Description                                |
| ---------------- | ----------- | ------------------------------------------ |
| T                | 400         | Total simulation time                      |
| TERMINAL\_WINDOW | 20          | Window size for terminal state calculation |
| THRESHOLD        | 0.1         | Protein concentration threshold            |
| LAMBDA\_VALUES   | \[0.5, 0.1] | Dosage attenuation factors                 |

## 4. Code Version / Commit Hash

- Source script: `runs/run_008/scripts/run_008_dosage_sensitivity.py`
- Source code: `src/ddc.py`
- Git commit: `0c0a6ce33108f715d5bede25fe5ee4946cbffca3` (2026-04-08)

## 5. Persistent Attenuation Model

**Formula**: P\_i(t+1) = (1 - δ\_p) \* P\_i(t) + λ \* γ\_i \* X\_i(t)

Where:

- δ\_p = 0.01 (protein decay rate)
- λ ∈ {0.5, 0.1} (dosage factor)
- γ\_i = 1.0 (gene-specific factor for target genes)
- X\_i(t) = TF input at time t

**Interpretation**:

- λ = 1.0: Baseline (no attenuation)
- λ = 0.5: 50% dosage reduction
- λ = 0.1: 90% dosage reduction

**Note**: KO (Knockout) results are sourced from run\_007, where KO was implemented by setting γ\_i = 0 (completely removing gene influence from the network), which is different from the continuous attenuation approach used here.

## 6. Collapse Classification

| Classification                  | Criteria                                            |
| ------------------------------- | --------------------------------------------------- |
| robust                          | collapse\_rate\_0.1 = 0                             |
| haploinsufficient\_like         | collapse\_rate\_0.5 > 0                             |
| intermediate\_dosage\_sensitive | collapse\_rate\_0.5 = 0 AND collapse\_rate\_0.1 > 0 |
| KO\_only\_sensitive             | collapse\_rate\_KO > 0 AND collapse\_rate\_0.1 = 0  |
| non\_essential                  | collapse\_rate\_KO = 0                              |

## 7. Collapse Rate Results

| Gene ID | λ=0.5 | λ=0.1 | KO   | Classification          |
| ------- | ----- | ----- | ---- | ----------------------- |
| 0       | 0.25  | 0.75  | 1.00 | haploinsufficient\_like |
| 1       | 0.50  | 0.50  | 1.00 | haploinsufficient\_like |
| 2       | 0.50  | 0.50  | 1.00 | haploinsufficient\_like |
| 3       | 0.25  | 0.75  | 1.00 | haploinsufficient\_like |
| 4       | 0.25  | 0.75  | 0.75 | haploinsufficient\_like |
| 5       | 0.25  | 0.50  | 0.75 | haploinsufficient\_like |

**Average Collapse Rates**:

- λ=0.5: 0.33
- λ=0.1: 0.63
- KO: 0.92

## 8. Output Files

### Trajectories

```
data/trajectories/
    seed_xxxx/
        baseline_X_traj.tsv
        baseline_P_traj.tsv
        baseline_Z_traj.tsv
        gene_i_lambda_0.5_X_traj.tsv
        gene_i_lambda_0.5_P_traj.tsv
        gene_i_lambda_0.5_Z_traj.tsv
        gene_i_lambda_0.1_X_traj.tsv
        gene_i_lambda_0.1_P_traj.tsv
        gene_i_lambda_0.1_Z_traj.tsv
```

### Response Metrics

```
data/response_metrics/
    attenuation_summary.tsv
```

### Aggregated Results

```
results/aggregated/
    TF_dosage_sensitivity.tsv
```

### Plots

```
plots/
    collapse_rate_vs_lambda.png
    trajectories/
        P_traj_seed_xxxx_lambda_0.5.png
        P_traj_seed_xxxx_lambda_0.1.png
```

