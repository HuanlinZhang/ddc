#!/usr/bin/env python3
"""
Shared utility functions for DDC Lite analysis runs.

These functions provide standardised implementations of the core
analysis primitives used across run_002, run_004, run_005,
run_grn_perturbation, and their variants.

All trajectory-analysis functions accept torch.Tensor (matching
run_002 / 004 / 005).  Matrix-building functions return numpy arrays.

Default thresholds assume b = 10.
"""

import numpy as np
import torch
import networkx as nx
from typing import Any, Dict, List, Optional, Tuple, Set

# ── ddc_model_a constants ─────────────────────────────────────────
# Imported lazily to avoid circular imports when __name__ == '__main__'.
# Callers that already have dma can pass G explicitly to functions
# that need it (e.g. build_A_from_world_dict).
G_DEFAULT: int = 20


# ═══════════════════════════════════════════════════════════════════
# Module-level threshold defaults (b = 10)
# ═══════════════════════════════════════════════════════════════════

EPSILON: float                     = 1e-2
CONVERGENCE_WINDOW: int            = 50
COLLAPSE_THRESHOLD: float          = 1e-1
DIVERGENCE_THRESHOLD: float        = 1e5
CLIPPING_FRAC_THRESHOLD: float     = 0.1
SLOW_CONVERGENCE_THRESHOLD: int    = 200


# ═══════════════════════════════════════════════════════════════════
# Tier 3 — Naming helpers
# ═══════════════════════════════════════════════════════════════════

def make_world_id(strength_key: str, sign_key: str, topo_idx: int) -> str:
    """Canonical world-id string: {strength_key}_{sign_key}_t{topo_idx:03d}."""
    return f"{strength_key}_{sign_key}_t{topo_idx:03d}"


def make_combo_key(strength_key: str, sign_key: str) -> str:
    """Combo key used for grouping: {strength_key}_{sign_key}."""
    return f"{strength_key}_{sign_key}"


# ═══════════════════════════════════════════════════════════════════
# Tier 2 — Matrix construction from world metadata
# ═══════════════════════════════════════════════════════════════════

def build_A_from_world_dict(world_dict: Dict[str, Any],
                            G: int = G_DEFAULT) -> Tuple[np.ndarray, np.ndarray]:
    """Build full transition matrix A and bias vector b from world_dict.

    P_graph convention (from ddc_model_a):
        P_graph[i] = list of genes that regulate gene i
        key   (i) = target (regulated gene)
        value (j) = regulators (genes that regulate i), i.e. j → i

    A[i,i] = 1 - delta[i]                    (self-decay)
    A[i,j] = sign_ij * strength_ij  (j≠i)    (regulation)

    Returns (A, b) where A is (G×G) and b is (G×1).

    Canonical implementation from run_005.
    """
    A = np.zeros((G, G))
    delta = world_dict['parameters']['delta']
    for i in range(G):
        A[i, i] = 1.0 - delta[i]

    P_graph = {int(k): v for k, v in world_dict['P_graph'].items()}
    edge_signs = {int(k): {int(sk): sv for sk, sv in v.items()}
                  for k, v in world_dict['edge_signs'].items()}
    edge_strengths_dict = world_dict.get('edge_strengths', {})
    edge_strengths = {}
    if edge_strengths_dict:
        edge_strengths = {int(k): {int(sk): sv for sk, sv in v.items()}
                          for k, v in edge_strengths_dict.items()}

    for i in range(G):
        for j in P_graph.get(i, []):
            sign = edge_signs[i][j]
            strength = edge_strengths.get(i, {}).get(j, 0.0)
            A[i, j] = sign * strength

    b_arr = np.array(world_dict['parameters']['b']).reshape(-1, 1)
    return A, b_arr


# ═══════════════════════════════════════════════════════════════════
# Tier 1 — Core analysis primitives
# ═══════════════════════════════════════════════════════════════════

def detect_equilibrium(
    X_traj: torch.Tensor,
    epsilon: float = EPSILON,
    window: int = CONVERGENCE_WINDOW,
) -> Dict[str, Any]:
    """Detect convergence in a trajectory.

    Scans forward from t=0 for the first window of *window* consecutive
    steps where ‖X(t+1)−X(t)‖₂ < *epsilon*.

    Returns dict with:
      converged         bool
      convergence_time  int   (–1 if not converged)
      equilibrium_magnitude  ‖X[-1]‖₂
      equilibrium_sparsity   fraction of genes with X[-1] < epsilon
    """
    t_steps = X_traj.shape[0] - 1
    G = X_traj.shape[1]
    converged = False
    conv_time = -1

    for t in range(t_steps - window + 1):
        all_consecutive = True
        for w in range(window):
            diff = float(torch.norm(X_traj[t + w + 1] - X_traj[t + w]).item())
            if diff >= epsilon:
                all_consecutive = False
                break
        if all_consecutive:
            converged = True
            conv_time = t
            break

    X_eq = X_traj[-1]
    eq_magnitude = float(torch.norm(X_eq).item())
    eq_sparsity = float((X_eq < epsilon).sum().item()) / G

    return {
        'converged': converged,
        'convergence_time': conv_time,
        'equilibrium_magnitude': eq_magnitude,
        'equilibrium_sparsity': eq_sparsity,
    }


def analyze_stability(
    X_traj: torch.Tensor,
    clip_count: torch.Tensor,
    total_clips: int,
    divergence_threshold: float = DIVERGENCE_THRESHOLD,
    collapse_threshold: float = COLLAPSE_THRESHOLD,
    clipping_frac_threshold: float = CLIPPING_FRAC_THRESHOLD,
    G: int = G_DEFAULT,
) -> Dict[str, Any]:
    """Analyze stability: divergence, numerical collapse, clipping.

    Returns dict with:
      bounded               bool
      max_expression         float
      final_mean_expression  float
      collapsed              bool
      divergence_existence   bool
      divergence_time        int | None
      numerical_collapse     bool
      clipping_dominated     bool
      total_clips            int
    """
    max_expression = float(X_traj.max().item())
    bounded = max_expression < divergence_threshold
    final_mean = float(X_traj[-1].mean().item())
    collapsed = final_mean < collapse_threshold
    total_steps = X_traj.shape[0]
    clipping_dominated = total_clips > total_steps * G * clipping_frac_threshold

    divergence_time = None
    if not bounded:
        for t in range(X_traj.shape[0]):
            if float(X_traj[t].max().item()) >= divergence_threshold:
                divergence_time = t
                break

    return {
        'bounded': bounded,
        'max_expression': max_expression,
        'final_mean_expression': final_mean,
        'collapsed': collapsed,
        'divergence_existence': not bounded,
        'divergence_time': divergence_time,
        'numerical_collapse': collapsed,
        'clipping_dominated': clipping_dominated,
        'total_clips': total_clips,
    }


def analyze_oscillation(
    X_traj: torch.Tensor,
    converged: bool,
    conv_time: int,
    epsilon: float = EPSILON,
    burn_in: int = 200,
    min_relative_amplitude: float = 0.01,
    damping_threshold: float = 0.05,
) -> Dict[str, Any]:
    """Detect oscillation via extrema-based analysis (run_002 algorithm).

    Identifies local peaks/troughs per gene and classifies oscillations
    as 'sustained' or 'damped' based on the median amplitude trend
    across all oscillatory genes.

    When converged, analysis is restricted to the window
    [burn_in, conv_time) to exclude the converged tail.

    Returns dict with:
      oscillation_exists  bool
      oscillation_type    'sustained' | 'damped' | 'none'
      amplitude           float
      frequency           float
      damping_rate        float | None
      oscillatory_genes   List[int]
    """
    T_total = X_traj.shape[0]

    if converged and conv_time > burn_in + 100:
        X = X_traj[burn_in:conv_time].numpy()
    else:
        X = X_traj[burn_in:].numpy()

    T, G_dim = X.shape
    if T < 50:
        return {
            'oscillation_exists': False,
            'oscillation_type': 'none',
            'amplitude': 0.0,
            'frequency': 0.0,
            'damping_rate': None,
            'oscillatory_genes': [],
        }

    oscillatory_genes: List[int] = []
    amplitudes: List[float] = []
    frequencies: List[float] = []
    damping_rates: List[float] = []

    for g in range(G_dim):
        x = X[:, g]
        signal_range = float(x.max() - x.min())
        if signal_range < epsilon:
            continue

        extrema = np.zeros(T, dtype=np.int8)
        for t in range(1, T - 1):
            if x[t] > x[t - 1] and x[t] > x[t + 1]:
                extrema[t] = 1
            elif x[t] < x[t - 1] and x[t] < x[t + 1]:
                extrema[t] = -1

        ext_idx = np.where(extrema != 0)[0]
        if len(ext_idx) < 3:
            continue

        peak_pairs = []
        for i in range(len(ext_idx) - 1):
            s, e = ext_idx[i], ext_idx[i + 1]
            delta = abs(float(x[e] - x[s]))
            peak_pairs.append(delta)

        if len(peak_pairs) < 2:
            continue

        gene_amplitude = float(np.median(peak_pairs))
        relative_amplitude = (gene_amplitude / signal_range
                              if signal_range > epsilon else 0.0)

        if relative_amplitude < min_relative_amplitude:
            continue

        gene_freq = len(ext_idx) / (2.0 * T)

        if len(peak_pairs) >= 4:
            mid = len(peak_pairs) // 2
            early = np.mean(peak_pairs[:mid])
            late = np.mean(peak_pairs[mid:])
            damping = float((early - late) / early) if early > epsilon else 0.0
        else:
            damping = 0.0

        oscillatory_genes.append(g)
        amplitudes.append(gene_amplitude)
        frequencies.append(gene_freq)
        damping_rates.append(damping)

    if not oscillatory_genes:
        return {
            'oscillation_exists': False,
            'oscillation_type': 'none',
            'amplitude': 0.0,
            'frequency': 0.0,
            'damping_rate': None,
            'oscillatory_genes': [],
        }

    avg_damping = float(np.median(damping_rates))
    osc_type = 'damped' if avg_damping > damping_threshold else 'sustained'

    return {
        'oscillation_exists': True,
        'oscillation_type': osc_type,
        'amplitude': float(np.mean(amplitudes)),
        'frequency': float(np.mean(frequencies)),
        'damping_rate': avg_damping,
        'oscillatory_genes': oscillatory_genes,
    }


def classify_attractor(
    eq: Dict[str, Any],
    st: Dict[str, Any],
    osc: Dict[str, Any],
    slow_convergence_threshold: int = SLOW_CONVERGENCE_THRESHOLD,
) -> str:
    """Classify attractor type (Type A–G) per run_002 convention.

    Priority order:
      Type E — Runaway divergence (not bounded)
      Type F — Numerical collapse   (collapsed)
      Type D — Sustained oscillation
      Type C — Damped oscillation
      Type A — Fast convergence      (converged, conv_time ≤ threshold)
      Type B — Slow convergence      (converged, conv_time > threshold)
      Type G — Others / ambiguous
    """
    if st['divergence_existence']:
        return 'Type E'
    if st['numerical_collapse']:
        return 'Type F'
    if osc['oscillation_type'] == 'sustained':
        return 'Type D'
    if osc['oscillation_type'] == 'damped':
        return 'Type C'
    if eq['converged']:
        conv_time = eq['convergence_time']
        if 0 <= conv_time <= slow_convergence_threshold:
            return 'Type A'
        return 'Type B'
    return 'Type G'


def compute_spectral_info(
    world_dict: Dict[str, Any],
    G: int = G_DEFAULT,
) -> Tuple[float, Dict[str, List[float]]]:
    """Compute spectral radius and eigenvalue info from world_dict.

    Builds the full transition matrix A via build_A_from_world_dict,
    then computes eigenvalues and spectral radius ρ(A).

    Returns (spectral_radius, eigen_info_dict).
    """
    A, _ = build_A_from_world_dict(world_dict, G=G)

    eigvals = np.linalg.eigvals(A)
    abs_vals = np.abs(eigvals)
    spectral_radius = float(np.max(abs_vals))
    eig_dict = {
        'real_parts': [float(ev.real) for ev in eigvals],
        'imag_parts': [float(ev.imag) for ev in eigvals],
        'abs_values': [float(abs(ev)) for ev in eigvals],
        'spectral_radius': spectral_radius,
    }
    return spectral_radius, eig_dict


# ═══════════════════════════════════════════════════════════════════
# Tier 1.3 — Topology helpers
# ═══════════════════════════════════════════════════════════════════

def compute_scc_from_pgraph(
    pg: Dict[int, List[int]],
    G: int = G_DEFAULT,
) -> Dict[str, Any]:
    """Compute strongly connected components from P_graph adjacency dict.

    P_graph convention (from ddc_model_a):
        P_graph[i] = list of genes that regulate gene i
        key   (i) = target (regulated gene)
        value (j) = regulators (genes that regulate i), i.e. j → i

    Uses networkx.strongly_connected_components internally.

    Returns dict with keys:
        comp_labels : list[int]   — length G, component ID per gene
        scc_sizes   : list[int]   — sizes sorted descending
        largest_scc_id : int      — component ID of the largest SCC
        n_scc       : int         — total number of SCCs

    Source: adapted from Phase 0 run_004 (phase0_regime_analysis.py) and
    run_005 (run_005_analysis.py), unified to networkx.
    """
    dg = nx.DiGraph()
    for target, regulators in pg.items():
        for source in regulators:
            dg.add_edge(int(source), int(target))

    # Add isolated nodes (genes with zero in/out edges) explicitly
    dg.add_nodes_from(range(G))

    scc_list = list(nx.strongly_connected_components(dg))
    scc_list.sort(key=len, reverse=True)

    comp_labels = [-1] * G
    for cid, comp in enumerate(scc_list):
        for node in comp:
            comp_labels[node] = cid

    return {
        'comp_labels': comp_labels,
        'scc_sizes': [len(c) for c in scc_list],
        'largest_scc_id': 0,
        'n_scc': len(scc_list),
    }


def hierarchical_cluster_relative_l2(
    states: List[np.ndarray],
    threshold: float = 1e-3,
    eps: float = EPSILON,
) -> Tuple[List[int], Optional[List[List[float]]], Optional[List[List[float]]]]:
    """Agglomerative average-linkage clustering on pairwise relative L2 distance.

    Parameters
    ----------
    states    : list of (G,) arrays — state vectors to cluster
    threshold : float — distance threshold for fcluster (criterion='distance')
    eps       : float — numerical epsilon for relative-L2 denominator

    Returns
    -------
    labels : list[int]        — 0-indexed cluster labels, length n
    abs_dm : (n,n) list[list[float]] | None — absolute L2 distance matrix
    rel_dm : (n,n) list[list[float]] | None — relative L2 distance matrix

    Source: adapted from run_005 (run_005_analysis.py §3 _relative_l2 + _hierarchical_cluster).
    """
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform

    n = len(states)
    if n <= 1:
        return [0] * n, None, None

    abs_dm = np.zeros((n, n))
    rel_dm = np.zeros((n, n))
    for i in range(n):
        ai = np.array(states[i])
        na = float(np.linalg.norm(ai))
        for j in range(i + 1, n):
            aj = np.array(states[j])
            d = float(np.linalg.norm(ai - aj))
            abs_dm[i, j] = abs_dm[j, i] = d
            nb = float(np.linalg.norm(aj))
            r = d / max(na, nb, eps)
            rel_dm[i, j] = rel_dm[j, i] = r

    condensed = squareform(rel_dm, checks=False)
    Z = linkage(condensed, method='average')
    labels_raw = fcluster(Z, t=threshold, criterion='distance') - 1  # 0-indexed

    return (
        labels_raw.tolist(),
        abs_dm.tolist() if abs_dm is not None else None,
        rel_dm.tolist() if rel_dm is not None else None,
    )


# ═══════════════════════════════════════════════════════════════════
# Tier 1.5 — Edge-perturbation helpers
# ═══════════════════════════════════════════════════════════════════

def get_edge_list(A: np.ndarray) -> List[Tuple[int, int]]:
    """Return list of (src, tgt) tuples for all non-self-loop edges in A.

    A non-self-loop edge is defined as A[i, j] != 0 with i != j.
    """
    G = A.shape[0]
    edges: List[Tuple[int, int]] = []
    for i in range(G):
        for j in range(G):
            if i != j and A[i, j] != 0.0:
                edges.append((i, j))
    return edges


def delete_random_edges(
    A: np.ndarray,
    n_delete: int,
    rng: np.random.RandomState,
) -> Tuple[np.ndarray, List[Tuple[int, int]]]:
    """Randomly delete n_delete edges from A.

    Edge selection is uniform over the current set of non-self-loop
    edges (i.e., entries A[i, j] != 0, i != j).

    Parameters
    ----------
    A        : (G, G) transition matrix
    n_delete : number of edges to zero out
    rng      : np.random.RandomState for reproducibility

    Returns
    -------
    A_perturbed : copy of A with n_delete off-diagonal entries zeroed
    deleted_edges : list of (src, tgt) tuples that were removed
    """
    edges = get_edge_list(A)
    n_available = len(edges)

    if n_delete > n_available:
        print(f'  WARNING: only {n_available} edges available, deleting all')
        n_delete = n_available

    idx = rng.choice(n_available, size=n_delete, replace=False)
    deleted = [edges[i] for i in idx]

    A_pert = A.copy()
    for src, tgt in deleted:
        A_pert[src, tgt] = 0.0

    return A_pert, deleted


def analyze_trajectory(
    X_traj: torch.Tensor,
    clip_count: torch.Tensor,
    total_clips: int,
    world_dict: Dict[str, Any],
    *,
    epsilon: float = EPSILON,
    window: int = CONVERGENCE_WINDOW,
    divergence_threshold: float = DIVERGENCE_THRESHOLD,
    collapse_threshold: float = COLLAPSE_THRESHOLD,
    clipping_frac_threshold: float = CLIPPING_FRAC_THRESHOLD,
    slow_convergence_threshold: int = SLOW_CONVERGENCE_THRESHOLD,
    G: int = G_DEFAULT,
) -> Dict[str, Any]:
    """Convenience: run detect_equilibrium + analyze_stability +
    analyze_oscillation + classify_attractor + compute_spectral_info.

    Returns a dict with keys 'equilibrium', 'stability', 'oscillation',
    'attractor_type', 'spectral_radius', 'eigenvalues'.
    """
    eq = detect_equilibrium(X_traj, epsilon=epsilon, window=window)
    st = analyze_stability(X_traj, clip_count, total_clips,
                           divergence_threshold=divergence_threshold,
                           collapse_threshold=collapse_threshold,
                           clipping_frac_threshold=clipping_frac_threshold,
                           G=G)
    osc = analyze_oscillation(X_traj, eq['converged'],
                              eq['convergence_time'], epsilon=epsilon)
    attractor_type = classify_attractor(
        eq, st, osc,
        slow_convergence_threshold=slow_convergence_threshold,
    )
    spectral_r, eig_dict = compute_spectral_info(world_dict, G=G)

    return {
        'equilibrium': eq,
        'stability': st,
        'oscillation': osc,
        'attractor_type': attractor_type,
        'spectral_radius': spectral_r,
        'eigenvalues': eig_dict,
    }
