# Execution Record - Run 007 Sampling Extension

**Generated:** 2026-04-17
**Purpose:** Engineering management document for 100% reproducibility and traceability

---

## 1. Selected Seeds (World Seeds)

**Sampling Method:** Stratified sampling based on N_active_TF levels
**Total Seeds:** 10

| Seed ID | N_active_TF Level | Baseline Regime |
| ------- | ------------------- | --------------- |
| 20      | low                 | steady          |
| 23      | low                 | steady          |
| 29      | mid                 | steady          |
| 33      | mid                 | collapse-like   |
| 49      | high                | collapse-like   |
| 52      | high                | steady          |
| 55      | high                | steady          |
| 62      | mid                 | steady          |
| 69      | low                 | steady          |
| 83      | mid                 | steady          |

**Data Source:** Seeds selected from combined pool of run_004/03_sampling_extension and run_005 steady worlds

---

## 2. Gene List (KO Test Pool)

### 2.1 Transcription Factor (TF) Genes
- **Gene IDs:** 0, 1, 2, 3, 4, 5
- **Total:** 6 TF genes

### 2.2 Epigenetic Modifier Genes
- **Gene IDs:** 17, 18, 19
- **Total:** 3 EPI genes

### 2.3 Total KO Test Pool
- **Total Genes:** 9 (6 TF + 3 EPI)

---

## 3. Simulation Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| T | 400 | Total simulation time steps |
| Terminal Window | t ∈ [T-20, T] | Final 20 time steps for terminal state assessment |
| N_active_TF Threshold | 1 | Threshold for classifying active TFs |
| X Activity Threshold | 0.1 | Threshold for gene activity in regime classification |

---

## 4. Code Version

| Item | Value |
|------|-------|
| Git Commit | 0c0a6ce33108f715d5bede25fe5ee4946cbffca3 |
| Commit Message | Add run_006: External Perturbation Response Analysis |
| Script | run_007_sampling_extension.py |
| Script Path | `/home/zhanghl/projects/ddc_github/runs/run_007_sampling_extension/scripts/` |

---

## 5. Deviation Metric Definition

| Metric | Definition |
|--------|------------|
| **Metric Used** | Mean L1 Norm per Time Step |
| **Formula** | `mean(|P_ko - P_baseline|.sum(dim=1))` |
| **Calculation** | Computed over full trajectory (all time steps) |
| **Output Column** | `deviation_score` |

**Note:** The formula computes the L1 norm (sum of absolute differences across all genes) at each time step, then takes the mean across all time steps.

---

## 6. Terminal Window Definition

| Definition | Value |
|------------|-------|
| Window Range | t ∈ [T-20, T] |
| Window Size | 20 time steps |
| Calculation Method | Mean over last 20 time steps |
| Usage | Terminal state assessment, regime classification, deviation scoring |

---

## 7. Output Structure

```
run_007_sampling_extension/
├── data/
│   ├── trajectories/          # Raw trajectory data (.tsv)
│   ├── regime_results/        # Regime classification results
│   ├── response_metrics/      # KO response metrics
│   └── surviving_TF/          # Surviving TF analysis
├── results/
│   ├── aggregated/            # Aggregated summary tables
│   ├── cross_run/             # Cross-run comparison results
│   ├── per_seed/              # Per-seed detailed results
│   ├── surviving_TF/           # Surviving TF comparison
│   └── surviving_TF_summary.md # Scientific summary report
├── plots/
│   ├── global_summary.png
│   ├── delta_P_over_time.png
│   ├── regime_transition/
│   ├── class_comparison/
│   ├── surviving_TF/
│   ├── trajectories/
│   └── TF_network_topology.png
├── scripts/
│   └── run_007_sampling_extension.py
├── summary.md                  # This summary
├── execution_record.md          # This document
└── surviving_TF_summary.md     # Scientific report
```

---

## 8. Notes

- Seeds 33 and 49 were classified as `collapse-like` based on T=400 baseline trajectory analysis
- All other seeds (20, 23, 29, 52, 55, 62, 69, 83) were classified as `steady`
- The 8 steady seeds underwent full KO perturbation analysis
- The 2 collapse-like seeds (33, 49) provided additional insight into gene effects in already-compromised systems
