# Surviving TF Analysis Summary Report

**Generated:** 2026-04-17
**Purpose:** Scientific discovery report explaining "Why can the system survive when some core genes are knocked out?"

---

## 1. Network Centrality Assessment

### 1.1 Surviving TF Network Metrics

| Gene ID | Collapse Rate | Out Degree | In Degree | SCC ID |
|---------|---------------|------------|-----------|--------|
| 0 | 0.500 | 2.0 | 1.8 | 0 |
| 1 | 0.750 | 1.8 | 2.2 | 0 |
| 2 | 0.875 | 2.2 | 2.5 | 0 |
| 3 | 0.500 | 2.0 | 1.2 | 0 |
| 4 | 0.750 | 2.3 | 2.3 | 0 |
| 5 | 0.875 | 2.1 | 2.4 | 0 |

### 1.2 Centrality Statistics
- **Mean Out Degree:** 2.07 ± 0.18
- **Mean In Degree:** 2.07 ± 0.49
- **Degree Range:** 1.8 - 2.5

### 1.3 Key Observation
Unlike typical expectations where surviving nodes would have **lower** network centrality, the surviving TFs in this analysis show **moderate to high degree connectivity**. Notably:
- Gene 2 (out_degree=2.2, in_degree=2.5) and Gene 5 (out_degree=2.1, in_degree=2.4) have the highest connectivity yet still survive in some worlds
- Gene 3 has the **lowest in_degree (1.2)** among surviving TFs, which is more consistent with the "low centrality" hypothesis

---

## 2. Strongly Connected Component (SCC) Membership

### 2.1 SCC Analysis

| Metric | Value |
|--------|-------|
| All TF SCC ID | 0 |
| Largest SCC Size | 6 (all TFs) |
| TFs Outside Largest SCC | 0 |

### 2.2 Key Finding
**All surviving TFs belong to the core strongly connected component (SCC).** This directly contradicts the hypothesis that "surviving TFs are non-core nodes outside the SCC."

The network topology analysis reveals that the TF regulatory network forms a single fully connected SCC where all 6 TF genes are mutually reachable, indicating a highly interconnected regulatory architecture.

---

## 3. Hypothesis Validation Conclusions

### 3.1 Original Hypothesis
> "Non-critical nodes / redundant nodes" hypothesis: Genes that survive KO are peripheral, low-centrality nodes not essential for network cohesion.

### 3.2 Evidence Assessment

| Hypothesis Aspect | Evidence | Supports Hypothesis? |
|-------------------|----------|---------------------|
| Low network centrality | Mixed: some surviving TFs (Gene 3) have low in_degree=1.2, but others have high connectivity (Gene 2, 5) | **Partially** |
| Outside SCC membership | All TFs are within the core SCC (size=6) | **No** |
| Redundant/Non-critical | All 6 TFs survive in at least one world, suggesting network has built-in redundancy | **Yes** |

### 3.3 Revised Interpretation

Based on the observed data, the "redundancy" hypothesis is better supported than the "non-critical node" hypothesis:

1. **Network Redundancy:** The fact that ALL tested TFs can survive in at least some worlds suggests the network has **functional redundancy** - multiple genes can compensate for the loss of any single TF.

2. **Collapse Rate Variation:** Even within the same SCC, different TFs show varying collapse rates (0.5 - 0.875), indicating that **position within the SCC matters** - genes with higher connectivity (Gene 2, 5) tend to have higher collapse rates.

3. **Context-Dependent Criticality:** The survival of TFs depends on both the gene's role AND the specific world configuration, suggesting **condition-dependent essentiality** rather than fixed criticality.

### 3.4 Conclusion

> **The observed phenomenon supports a modified "Redundant + Position-Dependent" hypothesis**: The TF network exhibits functional redundancy where no single TF is absolutely essential. However, genes with higher connectivity within the SCC tend to be more important across multiple world configurations. This explains why some core genes can be knocked out while the system survives - the network's interconnected architecture provides backup pathways for critical functions.

---

## 4. Data Sources

- **KO Response Data:** `results/per_seed/response_summary.tsv`
- **Surviving TF Metrics:** `data/surviving_TF/surviving_TF_network_metrics.tsv`
- **Network Topology:** `plots/TF_network_topology.png`
