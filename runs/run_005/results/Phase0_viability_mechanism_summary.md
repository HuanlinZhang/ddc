# Phase 0 Viability Mechanism Summary

## Executive Summary

This document synthesizes findings from four experimental tasks designed to understand the **viability boundary** of Phase 0 regulatory worlds — specifically, why some worlds collapse while others maintain steady-state expression.

**Central Question:** What determines whether a regulatory world collapses or survives?

***

## Task Overview

| Task   | Experiment                    | Key Question                                                           |
| ------ | ----------------------------- | ---------------------------------------------------------------------- |
| Task 1 | Resource Projection Ablation  | Does global resource constraint (`ΣPᵢ ≤ R_total`) determine viability? |
| Task 2 | Effective Gain Proxy Analysis | What dynamical metrics best distinguish collapse from steady-state?    |
| Task 3 | TF Support / Rescue Test      | Can increased TF support rescue collapsed worlds?                      |
| Task 4 | Initial Condition Scaling     | Does viability depend on initial expression levels?                    |

***

## Key Findings

### 1. Resource Projection Has Minimal Impact on Viability

**Finding:** Removing resource projection produced **no regime transitions**. All 7 collapse worlds remained collapsed, and all 4 steady-state worlds remained stable.

| Condition                     | Collapse | Steady |
| ----------------------------- | -------- | ------ |
| Original (with projection)    | 7        | 4      |
| Ablation (projection removed) | 7        | 4      |

**Interpretation:** The global resource constraint is **not the bottleneck** determining viability. Collapse is driven by regulatory dynamics, not resource limitation.

***

### 2. TF Decay Is the Primary Collapse Mechanism

**Finding:** All proxy metrics **effectively distinguish** collapse from steady-state worlds, in both original and projection_removed conditions:

| Condition | Metric | Collapse (n=7) | Steady (n=4) | Interpretation |
| --------- | ------ | -------------- | ------------ | -------------- |
| **Original** | mean(TFinput/K) | 0.245 ± 0.217 | 2.955 ± 1.470 | Collapse worlds have near-zero TF input |
| | TF Retention Ratio | 0.045 ± 0.063 | 0.444 ± 0.130 | Steady worlds retain ~44% of peak TF |
| | Terminal TFinput | ~1.0e-6 | 0.019 ± 0.009 | Collapse worlds have complete TF input collapse |
| | mean Hill Activation | 0.040 ± 0.025 | 0.479 ± 0.139 | TF regulatory function is lost in collapse |
| **Projection Removed** | mean(TFinput/K) | 0.432 ± 0.354 | 2.841 ± 1.137 | Same pattern as original |
| | TF Retention Ratio | 0.076 ± 0.070 | 0.363 ± 0.084 | Slightly higher but still discriminative |
| | Terminal TFinput | ~1.6e-5 | 0.014 ± 0.006 | Same pattern as original |
| | mean Hill Activation | 0.068 ± 0.048 | 0.421 ± 0.095 | Same pattern as original |

**Interpretation:** Collapse is driven by **TF decay** — after initial activation, the system cannot maintain TF expression levels. The consistency across both conditions confirms that **resource projection is not the determining factor**.

***

### 3. TF Support Cannot Rescue Collapse Worlds

**Finding:** Neither TF overexpression (ζ = 2–5×) nor TF clamping (threshold = 0.02–0.08) rescued any of 4 tested collapse worlds.

| Method                       | Experiments | Rescued | Success Rate |
| ---------------------------- | ----------- | ------- | ------------ |
| TF Overexpression (via mRNA) | 16          | 0       | 0%           |
| TF Clamp                     | 16          | 0       | 0%           |

**Interpretation:** Collapse is **robust and self-reinforcing**. Even when TF levels are artificially elevated through system-internal variables, the collapse phenotype persists. This suggests collapse is driven by **insufficient regulatory gain** — the regulatory system cannot sustain TF expression regardless of initial TF levels.

***

### 4. Collapse Is Not Initial-Condition Dependent

**Finding:** Scaling initial conditions (X₀, P₀) by γ ∈ {0.5, 1.0, 2.0, 5.0} produced **no regime transitions** across all 7 collapse seeds.

| γ Value | Regime Transitions |
| ------- | ------------------ |
| 0.5     | 0 / 7              |
| 1.0     | 0 / 7              |
| 2.0     | 0 / 7              |
| 5.0     | 0 / 7              |

**Interpretation:** The collapse attractor is **stable and independent of initial conditions**. The system does not exhibit multiple attractors based on initialization. Collapse is determined by intrinsic regulatory dynamics, not starting state.

***

## Answers to Core Questions

| Question                                            | Answer                                                                                                       |
| --------------------------------------------------- | ------------------------------------------------------------------------------------------------------------ |
| Does resource projection influence viability?       | **No** — removal of resource projection has no effect on regime classification                               |
| What proxies best capture regulatory gain?          | **TF Retention Ratio** is most intuitive; mean(TFinput/K) and Hill activation are also highly discriminative |
| Can TF support rescue collapse regimes?             | **No** — 32 experiments across 4 collapse seeds showed 0% rescue success                                     |
| Does viability depend on initial expression levels? | **No** — scaling initial conditions by 0.5× to 5× produces no regime transitions                             |

***

## Output Files

| Category           | Files                                                                                                            |
| ------------------ | ---------------------------------------------------------------------------------------------------------------- |
| **Task 1**         | `results/resource_ablation/regime_comparison.tsv`, `resource_ablation_analysis.md`                               |
| **Task 2**         | `results/gain_proxy_analysis/proxy_metrics.tsv`, `viability_proxy_summary.tsv`, `effective_gain_proxy_analysis.md` |
| **Task 3**         | `results/TF_rescue/rescue_summary.tsv`, `TF_rescue_analysis.md`                                                  |
| **Task 4**         | `results/initial_condition_scaling/scaling_summary.tsv`, `scaling_analysis.md`                                   |
| **Visualizations** | `plots/resource_ablation/`, `plots/gain_proxy_analysis/`, `plots/TF_rescue/`, `plots/initial_condition_scaling/` |

***

## Conclusions

Phase 0 regulatory viability is primarily determined by **intrinsic regulatory dynamics**, specifically the effective gain of the TF regulatory network. Key findings:

1. **Resource constraints are not limiting** — the ΣPᵢ ≤ R_total constraint does not determine whether a world collapses
2. **TF decay drives collapse** — collapse worlds fail to maintain TF expression after initial activation
3. **Collapse is robust and irreversible** — TF support interventions cannot rescue collapsed worlds
4. **Viability is an emergent property of regulatory dynamics** — initial conditions and resource allocation are not determining factors

This suggests that designing viable Phase 0 regulatory worlds requires attention to **regulatory parameters and feedback strength**, rather than resource allocation or initial expression levels.
