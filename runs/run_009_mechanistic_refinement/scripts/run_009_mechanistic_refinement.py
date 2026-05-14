#!/usr/bin/env python3
"""
Run 009 — Stage 0 Mechanistic Refinement

Three tasks from the feedback document:

  Task 1 — Hill Curve Positioning:
    Per-gene TFinput/K on the Hill curve, combined with network degree
    to distinguish "activated vs saturated" and "sensitive vs insensitive"

  Task 2 — Collapse Mode Classification:
    Gradual decay vs threshold crossing
    Find evidence of a "point of no return" in the TFinput timeseries

  Task 3 — Extended Simulation:
    Extend T=400 to T=800 for anomalous survival cases
    to determine slow-collapse vs quasi-steady

Author: zhanghl
Date: 2026-04-28
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')
from ddc import (
    sample_world, sample_initial_state, simulate_single_cell,
    World, TF_GENES, EPI_GENES, G, DTYPE, T as DEFAULT_T,
    compute_TFinput, normalize_protein,
    update_chromatin, update_mRNA, update_protein_raw,
    apply_resource_projection, update_fate,
)

BASE_DIR = Path('/home/zhanghl/projects/ddc_github/runs/run_009_mechanistic_refinement')
PLOTS_DIR = BASE_DIR / 'plots'
RESULTS_DIR = BASE_DIR / 'results'
DATA_DIR = BASE_DIR / 'data'

os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

STEADY_SEEDS = [2026, 4094, 4478, 7709]
LAMBDA_VALUES = [0.5, 0.1]
TERMINAL_WINDOW = 20
T_EXTENDED = 1000

WORLD_DIR = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_enable_resource_projection'
RUN008_TRAJ_DIR = '/home/zhanghl/projects/ddc_github/runs/run_008/data/trajectories'
RUN008_RESPONSE = '/home/zhanghl/projects/ddc_github/runs/run_008/data/response_metrics/attenuation_summary.tsv'


# ============================================================
# Common utilities
# ============================================================

def load_world(seed: int) -> World:
    traj = torch.load(f'{WORLD_DIR}/seed_{seed}_traj.pt', weights_only=False)
    w = World(seed)
    w.from_dict(traj['world'])
    return w


def load_traj(seed: int, prefix: str) -> np.ndarray:
    path = Path(RUN008_TRAJ_DIR) / f'seed_{seed}' / f'{prefix}_P_traj.tsv'
    if not path.exists():
        return None
    df = pd.read_csv(path, sep='\t')
    return df.drop(columns=['time']).values


def load_attenuation_summary() -> pd.DataFrame:
    return pd.read_csv(RUN008_RESPONSE, sep='\t')


def compute_out_degree(world: World) -> Dict[int, int]:
    deg = {g: 0 for g in TF_GENES}
    for i in range(G):
        for j in world.P_graph[i]:
            if j in TF_GENES:
                deg[j] += 1
    return deg


def hill_func(x, n=2.0):
    xn = x ** n
    return xn / (1.0 + xn)


def hill_slope(x, n=2.0):
    return n * x ** (n - 1) / (1.0 + x ** n) ** 2


# ============================================================
# Task 1 — Hill Curve Positioning
# ============================================================

def task1_hill_positioning():
    print("\n" + "=" * 60)
    print("Task 1 — Hill Curve Positioning")
    print("=" * 60)

    attn_df = load_attenuation_summary()
    tf_data = {}
    all_data = {}

    for seed in STEADY_SEEDS:
        world = load_world(seed)
        degree = compute_out_degree(world)
        P_baseline = load_traj(seed, 'baseline')

        tf_ratios = _full_trajectory_tfinput_k_per_gene(world, P_baseline, TF_GENES)
        all_ratios = _full_trajectory_tfinput_k_per_gene(world, P_baseline, list(range(G)))

        tf_gene_info = {}
        for g in TF_GENES:
            collapse_info = {}
            for lam in LAMBDA_VALUES:
                cond = attn_df[
                    (attn_df['seed'] == seed) &
                    (attn_df['gene_id'] == g) &
                    (attn_df['lambda'] == lam)
                ]
                collapse_info[lam] = (
                    bool(len(cond) > 0 and cond.iloc[0]['regime'] == 'collapse-like')
                )

            tf_gene_info[g] = {
                'degree': degree[g],
                'baseline_x': tf_ratios[g],
                'baseline_y': hill_func(tf_ratios[g]),
                'baseline_slope': hill_slope(tf_ratios[g]),
                'collapse': collapse_info,
            }

        tf_data[seed] = tf_gene_info
        all_data[seed] = all_ratios

    _plot_hill_positioning(tf_data)
    _plot_hill_all_genes(all_data)
    _plot_hill_mean_comparison(all_data)
    _save_hill_summary(tf_data)

    return tf_data


def _full_trajectory_tfinput_k_per_gene(world: World, P_traj: np.ndarray,
                                         gene_list: list) -> Dict[int, float]:
    T = P_traj.shape[0]
    ratios = {g: 0.0 for g in gene_list}
    for t in range(T):
        P = torch.tensor(P_traj[t], dtype=DTYPE)
        tilde_P = normalize_protein(P, world)
        tfinput = compute_TFinput(tilde_P, world)
        for g in gene_list:
            ratios[g] += float(tfinput[g] / world.K[g])
    return {g: v / T for g, v in ratios.items()}


def _best_offset(idx, _existing):
    offsets = [
        (15, 8), (-15, 12), (18, -10), (-18, 6),
        (10, -14), (-10, 14), (12, 12), (-12, -8),
    ]
    return offsets[idx % len(offsets)]


def _plot_hill_positioning(all_data: dict):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for idx, seed in enumerate(STEADY_SEEDS):
        ax = axes[idx // 2, idx % 2]
        gene_info = all_data[seed]

        all_xs = [gene_info[g]['baseline_x'] for g in TF_GENES]
        max_x = max(all_xs)
        x_vals = np.linspace(0, max_x * 1.15, 500)
        ax.plot(x_vals, hill_func(x_vals), 'k-', linewidth=1.5, alpha=0.4, zorder=0)

        high_slope_mask = (x_vals >= 0.3) & (x_vals <= 3.0)
        ax.fill_between(x_vals[high_slope_mask], 0, hill_func(x_vals)[high_slope_mask],
                        color='red', alpha=0.04)

        ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)
        ax.axvline(x=1.0, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)

        groups = {'both': [], 'only_01': [], 'none': []}
        for g in TF_GENES:
            if gene_info[g]['collapse'][0.5]:
                groups['both'].append(g)
            elif gene_info[g]['collapse'][0.1]:
                groups['only_01'].append(g)
            else:
                groups['none'].append(g)

        palette = {'both': '#E74C3C', 'only_01': '#F39C12', 'none': '#27AE60'}
        labels_map = {
            'both': 'collapse at \u03bb=0.5',
            'only_01': 'collapse at \u03bb=0.1 only',
            'none': 'survive \u03bb=0.1',
        }

        global_gi = 0
        for key, glist in groups.items():
            if not glist:
                continue
            xs = [gene_info[g]['baseline_x'] for g in glist]
            ys = [gene_info[g]['baseline_y'] for g in glist]
            ax.scatter(xs, ys, c=palette[key], s=120, zorder=5,
                       edgecolors='black', linewidth=0.8, label=labels_map[key])
            for g in glist:
                x0, y0 = gene_info[g]['baseline_x'], gene_info[g]['baseline_y']
                offset = _best_offset(global_gi, [])
                global_gi += 1
                ax.annotate(
                    f'TF{g}',
                    (x0, y0),
                    textcoords="offset points", xytext=offset,
                    fontsize=9, ha='center', va='center',
                    color=palette[key], fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                              edgecolor=palette[key], alpha=0.85),
                )

        mean_x = np.mean(all_xs)
        mean_y = hill_func(mean_x)
        ax.axvline(x=mean_x, color='#7F8C8D', linestyle='--', linewidth=1.2, alpha=0.7)
        ax.scatter([mean_x], [mean_y], marker='D', s=100, c='#7F8C8D',
                   edgecolors='black', linewidth=0.8, zorder=6, label='mean TFinput/K')
        ax.annotate(f'mean={mean_x:.1f}', xy=(mean_x, mean_y),
                    textcoords="offset points", xytext=(8, -8),
                    fontsize=9, color='#7F8C8D', fontweight='bold')

        ax.set_xlim(0, max_x * 1.15)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('TFinput / K', fontsize=12)
        ax.set_ylabel('Hill Activation', fontsize=12)
        ax.set_title(f'Seed {seed}', fontsize=14, fontweight='bold')
        ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
        ax.grid(True, alpha=0.2)

        ticks = list(ax.get_xticks())
        if 1.0 not in ticks:
            ticks.append(1.0)
            ticks.sort()
            ax.set_xticks(ticks)

    fig.suptitle('Task 1 — TF Positioning on the Hill Curve',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    path = PLOTS_DIR / 'task1_hill_curve_positioning.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def _plot_hill_all_genes(all_data: dict):
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    for idx, seed in enumerate(STEADY_SEEDS):
        ax = axes[idx // 2, idx % 2]
        ratios = all_data[seed]

        all_xs = list(ratios.values())
        max_x = max(all_xs)
        x_vals = np.linspace(0, max_x * 1.15, 500)
        ax.plot(x_vals, hill_func(x_vals), 'k-', linewidth=1.5, alpha=0.4, zorder=0)

        high_slope_mask = (x_vals >= 0.3) & (x_vals <= 3.0)
        ax.fill_between(x_vals[high_slope_mask], 0, hill_func(x_vals)[high_slope_mask],
                        color='red', alpha=0.04)

        ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)
        ax.axvline(x=1.0, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)

        tf_xs = [ratios[g] for g in TF_GENES]
        tf_ys = [hill_func(ratios[g]) for g in TF_GENES]
        non_tf_xs = [ratios[g] for g in range(G) if g not in TF_GENES]
        non_tf_ys = [hill_func(ratios[g]) for g in range(G) if g not in TF_GENES]

        ax.scatter(non_tf_xs, non_tf_ys, c='#BDC3C7', s=30, zorder=4,
                   edgecolors='none', label='non-TF (44 genes)')
        ax.scatter(tf_xs, tf_ys, c='#E74C3C', s=100, zorder=5,
                   edgecolors='black', linewidth=0.8, label='TF (6 genes)')

        for g in TF_GENES:
            x0, y0 = ratios[g], hill_func(ratios[g])
            ax.annotate(f'TF{g}', (x0, y0),
                        textcoords="offset points", xytext=(8, 6),
                        fontsize=9, ha='center', va='center',
                        color='#E74C3C', fontweight='bold')

        mean_all = np.mean(all_xs)
        mean_all_y = hill_func(mean_all)
        ax.axvline(x=mean_all, color='#7F8C8D', linestyle='--', linewidth=1.2, alpha=0.7)
        ax.scatter([mean_all], [mean_all_y], marker='D', s=100, c='#7F8C8D',
                   edgecolors='black', linewidth=0.8, zorder=6)
        ax.annotate(f'mean={mean_all:.1f}', xy=(mean_all, mean_all_y),
                    textcoords="offset points", xytext=(8, -8),
                    fontsize=9, color='#7F8C8D', fontweight='bold')

        ax.set_xlim(0, max_x * 1.15)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xlabel('TFinput / K', fontsize=12)
        ax.set_ylabel('Hill Activation', fontsize=12)
        ax.set_title(f'Seed {seed} (all 50 genes)', fontsize=14, fontweight='bold')
        ax.legend(loc='lower right', fontsize=9, framealpha=0.9)
        ax.grid(True, alpha=0.2)

        ticks = list(ax.get_xticks())
        if 1.0 not in ticks:
            ticks.append(1.0)
            ticks.sort()
            ax.set_xticks(ticks)

    fig.suptitle('Task 1 — All 50 Genes on the Hill Curve',
                 fontsize=14, fontweight='bold', y=1.01)
    plt.tight_layout()
    path = PLOTS_DIR / 'task1_hill_curve_all_genes.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def _plot_hill_mean_comparison(all_data: dict):
    fig, ax = plt.subplots(figsize=(10, 7))

    means = {}
    for seed in STEADY_SEEDS:
        ratios = all_data[seed]
        means[seed] = np.mean(list(ratios.values()))

    max_x = max(means.values()) * 1.3
    x_vals = np.linspace(0, max_x, 500)
    ax.plot(x_vals, hill_func(x_vals), 'k-', linewidth=1.5, alpha=0.4, zorder=0)

    high_slope_mask = (x_vals >= 0.3) & (x_vals <= 3.0)
    ax.fill_between(x_vals[high_slope_mask], 0, hill_func(x_vals)[high_slope_mask],
                    color='red', alpha=0.04)

    ax.axhline(y=0.5, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)
    ax.axvline(x=1.0, color='gray', linestyle=':', alpha=0.4, linewidth=0.8)

    colors = ['#E74C3C', '#3498DB', '#2ECC71', '#F39C12']
    for idx, seed in enumerate(STEADY_SEEDS):
        mean_val = means[seed]
        mean_y = hill_func(mean_val)
        ax.scatter([mean_val], [mean_y], marker='D', s=100, c=colors[idx],
                   edgecolors='black', linewidth=1.0, zorder=6,
                   label=f'Seed {seed}')
        ax.annotate(f'{mean_val:.1f}', (mean_val, mean_y),
                    textcoords="offset points", xytext=(-14, 10),
                    fontsize=11, ha='center', va='center',
                    color=colors[idx], fontweight='bold')

    ticks = list(ax.get_xticks())
    if 1.0 not in ticks:
        ticks.append(1.0)
        ticks.sort()
        ax.set_xticks(ticks)

    ax.set_xlim(0, max_x)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel('mean(TFinput / K)', fontsize=12)
    ax.set_ylabel('Hill Activation', fontsize=12)
    ax.set_title(
        'Mean TFinput/K Position Across 4 Steady Worlds\n'
        '(mean computed over all 50 genes, not only TF genes)',
        fontsize=12, fontweight='bold')
    ax.legend(loc='lower right', fontsize=10, framealpha=0.9)
    ax.grid(True, alpha=0.2)

    fig.tight_layout()
    path = PLOTS_DIR / 'task1_hill_mean_comparison.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def _save_hill_summary(all_data: dict):
    rows = []
    for seed in STEADY_SEEDS:
        for g in TF_GENES:
            info = all_data[seed][g]
            rows.append({
                'seed': seed, 'gene_id': g, 'degree': info['degree'],
                'TFinput_over_K': round(info['baseline_x'], 4),
                'hill_activation': round(info['baseline_y'], 4),
                'hill_slope': round(info['baseline_slope'], 4),
                'collapse_0.5': info['collapse'][0.5],
                'collapse_0.1': info['collapse'][0.1],
            })
    df = pd.DataFrame(rows)
    path = RESULTS_DIR / 'task1_hill_position_summary.tsv'
    df.to_csv(path, sep='\t', index=False)
    print(f"  Saved: {path}")


# ============================================================
# Task 2 — Collapse Mode Classification (gradual vs crossing)
# ============================================================

def task2_collapse_mode():
    print("\n" + "=" * 60)
    print("Task 2 — Collapse Mode Classification")
    print("=" * 60)

    attn_df = load_attenuation_summary()
    all_results = []

    for seed in STEADY_SEEDS:
        world = load_world(seed)
        for g in TF_GENES:
            for lam in LAMBDA_VALUES:
                cond = attn_df[
                    (attn_df['seed'] == seed) &
                    (attn_df['gene_id'] == g) &
                    (attn_df['lambda'] == lam)
                ]
                if len(cond) == 0:
                    continue
                did_collapse = cond.iloc[0]['regime'] == 'collapse-like'

                P_traj = load_traj(seed, f'gene_{g}_lambda_{lam}')
                if P_traj is None:
                    continue

                tfinput_series = _compute_tfcore_mean_tfinput_timeseries(world, P_traj)

                mode, crossing_t, final_ratio = _classify_collapse_mode(
                    tfinput_series, did_collapse
                )

                all_results.append({
                    'seed': seed, 'gene_id': g, 'lambda': lam,
                    'did_collapse': did_collapse,
                    'mode': mode,
                    'crossing_time': crossing_t,
                    'terminal_core_tfinput_k': round(final_ratio, 4),
                })

    df = pd.DataFrame(all_results)
    path = RESULTS_DIR / 'task2_collapse_mode.tsv'
    df.to_csv(path, sep='\t', index=False)
    print(f"  Saved: {path}")

    _plot_collapse_mode_summary(df)

    return df


def _compute_tfcore_mean_tfinput_timeseries(world: World, P_traj: np.ndarray) -> np.ndarray:
    """Compute mean TFinput/K across all TF genes for collapse classification.

    Using the average across all TF genes better reflects the overall health
    of the regulatory network, rather than the specific input to the perturbed
    gene which can drop immediately due to attenuation.
    """
    from ddc import TF_GENES
    T_steps = P_traj.shape[0]
    ratios = np.zeros(T_steps)
    for t in range(T_steps):
        P = torch.tensor(P_traj[t], dtype=DTYPE)
        tilde_P = normalize_protein(P, world)
        tfinput = compute_TFinput(tilde_P, world)
        # Average TFinput/K across all TF genes
        mean_ratio = float((tfinput[TF_GENES] / world.K[TF_GENES]).mean())
        ratios[t] = mean_ratio
    return ratios


def _classify_collapse_mode(ratios: np.ndarray, did_collapse: bool):
    final_ratio = ratios[-TERMINAL_WINDOW:].mean()
    CROSS_THRESHOLD = 0.3
    SUSTAINED_WINDOW = 10

    if not did_collapse:
        return 'no_collapse', -1, final_ratio

    smoothing = np.convolve(ratios, np.ones(5) / 5, mode='valid')
    crossing_t = None
    for t in range(len(smoothing) - SUSTAINED_WINDOW):
        if smoothing[t] < CROSS_THRESHOLD:
            if all(smoothing[t:t + SUSTAINED_WINDOW] < CROSS_THRESHOLD):
                # 为什么是 t + 2？ 
                # 因为前面的滑动窗口大小是 5，
                # 索引 t 的值实际上是由原数组 [t, t+1, t+2, t+3, t+4] 算出来的
                # 所以原数组真正的中心发生突变的时刻是 t + 2
                crossing_t = t + 2
                break

    if crossing_t is None:
        if final_ratio < CROSS_THRESHOLD:
            return 'gradual_decay', -1, final_ratio
        else:
            return 'unclassified', -1, final_ratio

    pre_decay_len = crossing_t
    post_decay_len = len(ratios) - crossing_t
    if pre_decay_len > 30 and post_decay_len > 20:
        return 'threshold_crossing', int(crossing_t), final_ratio
    else:
        return 'rapid_collapse', int(crossing_t), final_ratio


def _plot_collapse_mode_summary(df: pd.DataFrame):
    collapsed = df[df['did_collapse']]
    counts = collapsed['mode'].value_counts()
    modes = ['threshold_crossing', 'gradual_decay', 'rapid_collapse']
    colors = ['#E74C3C', '#F39C12', '#3498DB']
    labels = ['Threshold\nCrossing', 'Gradual\nDecay', 'Rapid\nCollapse']

    fig, ax = plt.subplots(figsize=(8, 5))
    vals = [counts.get(m, 0) for m in modes]
    bars = ax.bar(range(len(modes)), vals, color=colors, alpha=0.85)
    ax.bar_label(bars, fontsize=12, padding=3)
    ax.set_xticks(range(len(modes)))
    ax.set_xticklabels(labels, fontsize=11)
    ax.set_ylabel('Number of (seed, gene, λ) cases', fontsize=11)
    ax.set_title('Task 2 — Collapse Mode Distribution', fontsize=13, fontweight='bold')
    ax.set_ylim(0, max(vals) * 1.2)
    ax.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    path = PLOTS_DIR / 'task2_collapse_mode.png'
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


# ============================================================
# Task 3 — Extended KO Simulation (T=800)
#
# Anomalous survivors (steady after KO at T=400):
#   From run_007:
#     seed 4094, KO gene 4
#     seed 4094, KO gene 5
#   From run_007_sampling_extension:
#     seed 20,   KO gene 2
#     seed 20,   KO gene 3
#     seed 33,   KO gene 3
#
# Extend to T=800 or 1200 to determine slow-collapse vs quasi-steady.
# ============================================================

RUN007_RESPONSE = '/home/zhanghl/projects/ddc_github/runs/run_007/data/response_metrics/response_summary.tsv'
RUN007SE_RESPONSE = '/home/zhanghl/projects/ddc_github/runs/run_007_sampling_extension/data/response_metrics/response_summary.tsv'
WORLD_DIR_EXT = '/home/zhanghl/projects/ddc_github/test_convergence/v1_2_run_004_100seeds'


def task3_extended_simulation():
    """Extend KO to T=800 for anomalous survivors from run_007 + sampling_extension."""
    print("\n" + "=" * 60)
    print("Task 3 — Extended KO Simulation (T=800)")
    print("=" * 60)

    cases = [
        {'seed': 4094, 'gene_id': 4,  'source': 'run_007',       'T': 800},
        {'seed': 4094, 'gene_id': 5,  'source': 'run_007',       'T': 800},
        {'seed': 20,   'gene_id': 2,  'source': 'run_007_se',    'T': 800},
        {'seed': 20,   'gene_id': 3,  'source': 'run_007_se',    'T': 1200},
        {'seed': 33,   'gene_id': 3,  'source': 'run_007_se',    'T': 800},
    ]

    run007 = pd.read_csv(RUN007_RESPONSE, sep='\t')
    run007_tf = run007[run007['gene_type'] == 'TF']
    run007se = pd.read_csv(RUN007SE_RESPONSE, sep='\t')
    run007se_tf = run007se[run007se['gene_type'] == 'TF']

    results = []
    for case in cases:
        seed, gene_idx, src, t_ext = case['seed'], case['gene_id'], case['source'], case['T']

        if src == 'run_007':
            df_tf = run007_tf
        else:
            df_tf = run007se_tf

        r7 = df_tf[(df_tf['seed'] == seed) & (df_tf['gene_id'] == gene_idx)]
        was_steady_T400 = (len(r7) > 0 and
                           r7.iloc[0]['ko_regime'] == 'steady')
        baseline_regime_T400 = r7.iloc[0]['baseline_regime'] if len(r7) > 0 else '?'

        print(f"\n  Simulating: seed={seed}, KO gene={gene_idx}, "
              f"src={src}, T400={r7.iloc[0]['ko_regime'] if len(r7) > 0 else '?'}, "
              f"T={t_ext}")

        world = _load_world_for_ko(seed, src)
        X0, P0, Z0, N0 = sample_initial_state(seed + 1, world)

        world_perturbed = World(seed)
        world_perturbed.from_dict(_load_raw_world_dict(seed, src))
        world_perturbed.rho[gene_idx] = 0.0

        traj = simulate_single_cell(world_perturbed, X0, P0, Z0, N0,
                                     t_ext)

        X_terminal = traj['X_traj'][-1, :]
        n_active = (X_terminal > 0.1).sum().item()
        final_regime = 'steady' if n_active > 1 else 'collapse-like'

        if final_regime == 'steady':
            classification = 'quasi_steady'
        elif was_steady_T400:
            classification = 'slow_collapse'
        else:
            classification = 'confirmed_collapse'

        results.append({
            'seed': seed, 'gene_id': gene_idx, 'perturbation': 'KO',
            'source': src,
            'baseline_regime_T400': baseline_regime_T400,
            'T': t_ext,
            'was_steady_T400': was_steady_T400,
            'final_regime': final_regime,
            'n_active_genes': n_active,
            'classification': classification,
        })

        _save_extended_trajectory(traj, seed, gene_idx, 'KO', t_ext)

    df = pd.DataFrame(results)
    # Keep only genuine anomalous survivors (steady at T=400)
    df = df[df['was_steady_T400'] == True].copy()
    path = RESULTS_DIR / 'task3_extended_simulation.tsv'
    df.to_csv(path, sep='\t', index=False)
    print(f"\n  Saved: {path}")
    for _, r in df.iterrows():
        print(f"    Seed {int(r['seed'])}, KO Gene {int(r['gene_id'])} "
              f"({r['source']}): {r['final_regime']} → {r['classification']}")

    _plot_extended_trajectories(df)
    return df


def _load_raw_world_dict(seed: int, source: str) -> dict:
    if source == 'run_007':
        return torch.load(f'{WORLD_DIR}/seed_{seed}_traj.pt', weights_only=False)['world']
    else:
        return torch.load(f'{WORLD_DIR_EXT}/seed_{seed}_traj.pt', weights_only=False)['world']


def _load_world_for_ko(seed: int, source: str) -> World:
    w = World(seed)
    w.from_dict(_load_raw_world_dict(seed, source))
    return w


def _save_extended_trajectory(traj: dict, seed: int, gene_idx: int,
                                pert_label: str, t_val: int):
    out_dir = DATA_DIR / f'seed_{seed}'
    os.makedirs(out_dir, exist_ok=True)
    for key in ['X_traj', 'P_traj', 'Z_traj']:
        data = traj[key].numpy()
        gene_cols = [f'gene_{i}' for i in range(data.shape[1])]
        df = pd.DataFrame(data, columns=gene_cols)
        df.insert(0, 'time', np.arange(data.shape[0]))
        prefix = f'KO_gene_{gene_idx}_T{t_val}'
        df.to_csv(out_dir / f'{prefix}_{key}.tsv', sep='\t', index=False)


def _plot_extended_trajectories(results_df: pd.DataFrame):
    from matplotlib.lines import Line2D

    fig, axes = plt.subplots(3, 2, figsize=(18, 18))
    axes_flat = axes.flat

    gene_legend = [
        Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
        Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
        Line2D([0], [0], color='gray', label='Others'),
    ]

    plot_df = results_df[results_df['was_steady_T400'] == True]
    for idx, (_, row) in enumerate(plot_df.iterrows()):
        seed, gene_idx = int(row['seed']), int(row['gene_id'])
        t_row = int(row['T'])
        t_orig = 400
        ax = axes_flat[idx]

        tsv_path = (DATA_DIR / f'seed_{seed}' /
                    f'KO_gene_{gene_idx}_T{t_row}_P_traj.tsv')
        if not tsv_path.exists():
            ax.text(0.5, 0.5, 'Data not found', transform=ax.transAxes,
                    ha='center', va='center')
            continue

        traj_df = pd.read_csv(tsv_path, sep='\t')
        data = traj_df.drop(columns=['time']).values
        T_len = data.shape[0]
        t = np.arange(1, T_len)

        y_max = data[1:].max()

        for gg in range(data.shape[1]):
            if gg in TF_GENES:
                color, lw, alpha, z = plt.cm.Set1(0), 1.2, 0.85, 2
            elif gg in EPI_GENES:
                color, lw, alpha, z = plt.cm.Set1(3), 1.2, 0.85, 2
            else:
                color, lw, alpha, z = 'gray', 0.6, 0.35, 1
            ax.plot(t, data[1:, gg], color=color, linewidth=lw,
                    alpha=alpha, zorder=z)

        ax.axvline(x=t_orig, color='#E74C3C', linestyle='--',
                   alpha=0.6, linewidth=1.5,
                   label=f'T={t_orig} (original)')

        ax.set_ylim(0, y_max * 1.05)
        ax.set_xlim(0, T_len - 1)
        ax.set_xlabel('Time', fontsize=11)
        ax.set_ylabel('Protein', fontsize=11)
        ax.set_title(f'Seed {seed}, KO Gene {gene_idx}, T={t_row} '
                     f'({row["classification"]})',
                     fontsize=12, fontweight='bold')

        handles = gene_legend + [ax.get_legend_handles_labels()[0][-1]]
        labels_comb = ['TF', 'Epigenetics', 'Others',
                       f'T={t_orig}']
        ax.legend(handles=handles, labels=labels_comb, loc='upper right',
                  fontsize=8, framealpha=0.9)
        ax.grid(True, alpha=0.3)

    # Hide unused subplot (5 cases → 3×2 layout)
    for i in range(len(plot_df), 6):
        axes_flat[i].set_visible(False)

    fig.suptitle('Task 3 — Extended KO Simulation (T=800)\n'
                 'Anomalous KO survivors (steady at T=400, extended to T=800)',
                 fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    out_path = PLOTS_DIR / 'task3_extended_trajectories.png'
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


# ============================================================
# Main
# ============================================================

def main():
    print("\n" + "=" * 60)
    print("Run 009 — Stage 0 Mechanistic Refinement")
    print("=" * 60)

    task1_hill_positioning()
    task2_collapse_mode()
    task3_extended_simulation()

    print("\n" + "=" * 60)
    print("Run 009 completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()
