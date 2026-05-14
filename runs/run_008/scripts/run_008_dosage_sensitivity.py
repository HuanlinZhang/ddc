#!/usr/bin/env python3
"""
Run 008 - TF Dosage Sensitivity Analysis

Analyzes how the system responds to persistent partial attenuation of TF genes,
where protein synthesis rate is reduced by factor λ instead of complete knockout.

Intervention formula:
    P_i(t+1) = (1 - δ_p) * P_i(t) + λ * γ_i * X_i(t)

This maintains a new steady-state at ~λ times the original protein level,
rather than exponential decay to zero.

Author: zhanghl
Date: 2026-04-20
"""

import os
import sys
import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from decimal import Decimal, ROUND_HALF_UP
from typing import Dict, List, Tuple
from pathlib import Path

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')
from ddc import (
    sample_world, sample_initial_state, simulate_single_cell,
    World, Tensor, T as DEFAULT_T, G, DTYPE, TF_GENES, EPI_GENES,
    normalize_protein, compute_TFinput, update_chromatin,
    update_mRNA, update_protein_raw, apply_resource_projection, update_fate
)

BASE_DIR = Path('/home/zhanghl/projects/ddc_github/runs/run_008')
DATA_DIR = BASE_DIR / 'data'
RESULTS_DIR = BASE_DIR / 'results'
PLOTS_DIR = BASE_DIR / 'plots'

RUN_007_SUMMARY = Path('/home/zhanghl/projects/ddc_github/runs/run_007/data/response_metrics/response_summary.tsv')
df = pd.read_csv(RUN_007_SUMMARY, sep='\t')
steady_seeds_all = df.groupby('seed')['baseline_regime'].first().reset_index()
STEADY_SEEDS = steady_seeds_all[steady_seeds_all['baseline_regime'] == 'steady']['seed'].tolist()[:4]
LAMBDA_VALUES = [0.5, 0.1]
T = 400
TERMINAL_WINDOW = 20


def setup_directories():
    for dir_path in [DATA_DIR, DATA_DIR / 'trajectories', DATA_DIR / 'response_metrics',
                     RESULTS_DIR, RESULTS_DIR / 'aggregated', PLOTS_DIR]:
        dir_path.mkdir(parents=True, exist_ok=True)


def classify_regime(X_traj: Tensor) -> str:
    THRESHOLD = 0.1
    N_ACTIVE_THRESHOLD = 1
    X_terminal = X_traj[-1, :]
    n_active = (X_terminal > THRESHOLD).sum().item()
    return 'steady' if n_active > N_ACTIVE_THRESHOLD else 'collapse-like'


def simulate_with_attenuation(world: World, X0: Tensor, P0: Tensor, Z0: Tensor, N0: float,
                              gene_idx: int, lambda_val: float, t_steps: int = T) -> Dict[str, Tensor]:
    """
    Simulate with persistent TF attenuation.

    For target gene i:
        P_i(t+1) = (1 - δ_p) * P_i(t) + λ * γ_i * X_i(t)

    Other genes follow normal dynamics.
    """
    X, P, Z, N = X0.clone(), P0.clone(), Z0.clone(), N0

    X_traj = torch.zeros((t_steps + 1, G), dtype=DTYPE)
    P_traj = torch.zeros((t_steps + 1, G), dtype=DTYPE)
    Z_traj = torch.zeros((t_steps + 1, G), dtype=DTYPE)
    N_traj = torch.zeros((t_steps + 1,), dtype=DTYPE)

    X_traj[0], P_traj[0], Z_traj[0], N_traj[0] = X, P, Z, N

    for t in range(t_steps):
        tilde_P = normalize_protein(P, world)
        TFinput = compute_TFinput(tilde_P, world)
        Z_next = update_chromatin(tilde_P, world)
        X_next = update_mRNA(X, Z, TFinput, world)

        P_raw = update_protein_raw(P, X, world)
        P_raw[gene_idx] = (1 - world.delta_p[gene_idx]) * P[gene_idx] + lambda_val * world.gamma[gene_idx] * X[gene_idx]
        P_next = apply_resource_projection(P_raw, world)

        N_next = update_fate(N, world)

        X, P, Z, N = X_next, P_next, Z_next, N_next
        X_traj[t + 1], P_traj[t + 1], Z_traj[t + 1], N_traj[t + 1] = X, P, Z, N

    return {'X_traj': X_traj, 'P_traj': P_traj, 'Z_traj': Z_traj, 'N_traj': N_traj}


def save_trajectory_tsv(traj: Dict[str, Tensor], seed_dir: str, prefix: str):
    os.makedirs(seed_dir, exist_ok=True)
    T = traj['X_traj'].shape[0]
    time_col = np.arange(T).reshape(-1, 1)
    for key in ['X_traj', 'P_traj', 'Z_traj']:
        data = traj[key].numpy()
        gene_cols = np.array([f'gene_{i}' for i in range(data.shape[1])])
        df = pd.DataFrame(data, columns=gene_cols)
        df.insert(0, 'time', time_col)
        df.to_csv(os.path.join(seed_dir, f'{prefix}_{key}.tsv'), sep='\t', index=False)


def load_trajectory_tsv(seed_dir: str, prefix: str) -> Dict[str, np.ndarray]:
    traj = {}
    for key in ['X_traj', 'P_traj', 'Z_traj']:
        df = pd.read_csv(os.path.join(seed_dir, f'{prefix}_{key}.tsv'), sep='\t')
        traj[key] = df.drop(columns=['time']).values
    return traj


def load_run_007_ko_trajectory(seed_dir: str, gene_idx: int) -> Dict[str, np.ndarray]:
    data = {}
    for traj_type in ['P_traj', 'X_traj', 'Z_traj']:
        file_path = os.path.join(seed_dir, f'KO_gene_{gene_idx}_{traj_type}.tsv')
        if os.path.exists(file_path):
            traj_data = []
            with open(file_path, 'r') as f:
                header = f.readline()
                for line in f:
                    parts = line.strip().split('\t')
                    traj_data.append([float(x) for x in parts[1:]])
            data[traj_type] = np.array(traj_data)
    return data


def run_simulation():
    print("\n" + "=" * 60)
    print("Run 008 - TF Dosage Sensitivity Analysis")
    print("=" * 60)

    setup_directories()

    results = []

    for seed in STEADY_SEEDS:
        print(f"\nProcessing seed {seed}...")

        world = sample_world(seed)
        X0, P0, Z0, N0 = sample_initial_state(seed + 1, world)

        seed_dir = DATA_DIR / 'trajectories' / f'seed_{seed}'
        os.makedirs(seed_dir, exist_ok=True)

        baseline_traj = simulate_single_cell(world, X0.clone(), P0.clone(), Z0.clone(), N0, T)
        baseline_regime = classify_regime(baseline_traj['X_traj'])
        save_trajectory_tsv(baseline_traj, str(seed_dir), 'baseline')

        baseline_terminal_P = baseline_traj['P_traj'][-TERMINAL_WINDOW:].mean(dim=0)

        for gene_idx in TF_GENES:
            for lambda_val in LAMBDA_VALUES:
                traj = simulate_with_attenuation(world, X0.clone(), P0.clone(), Z0.clone(), N0,
                                                gene_idx, lambda_val)
                regime = classify_regime(traj['X_traj'])

                prefix = f'gene_{gene_idx}_lambda_{lambda_val}'
                save_trajectory_tsv(traj, str(seed_dir), prefix)

                terminal_P = traj['P_traj'][-TERMINAL_WINDOW:].mean(dim=0)
                delta_terminal_P = (terminal_P - baseline_terminal_P).mean().item()

                results.append({
                    'seed': seed,
                    'gene_id': gene_idx,
                    'lambda': lambda_val,
                    'regime': regime,
                    'mean_P_terminal': terminal_P.mean().item(),
                    'delta_mean_P_terminal': delta_terminal_P,
                    'baseline_regime': baseline_regime
                })

                print(f"  Gene {gene_idx}, λ={lambda_val}: {regime}")

    results_df = pd.DataFrame(results)
    output_file = DATA_DIR / 'response_metrics' / 'attenuation_summary.tsv'
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved: {output_file}")

    return results_df


def analyze_results(results_df: pd.DataFrame):
    print("\n" + "=" * 60)
    print("Analyzing Results")
    print("=" * 60)

    gene_results = []

    for gene_idx in TF_GENES:
        gene_data = results_df[results_df['gene_id'] == gene_idx]

        steady_baseline_data = gene_data[gene_data['baseline_regime'] == 'steady']

        for lambda_val in LAMBDA_VALUES:
            lambda_data = steady_baseline_data[steady_baseline_data['lambda'] == lambda_val]
            n_collapse = (lambda_data['regime'] == 'collapse-like').sum()
            n_steady = (lambda_data['regime'] == 'steady').sum()
            n_total = len(lambda_data)

            if n_total > 0:
                collapse_rate = n_collapse / n_total
            else:
                collapse_rate = 0.0

            gene_results.append({
                'gene_id': gene_idx,
                'lambda': lambda_val,
                'collapse_count': n_collapse,
                'steady_count': n_steady,
                'collapse_rate': collapse_rate
            })

    gene_results_df = pd.DataFrame(gene_results)
    pivot_0_5 = gene_results_df[gene_results_df['lambda'] == 0.5].set_index('gene_id')['collapse_rate']
    pivot_0_1 = gene_results_df[gene_results_df['lambda'] == 0.1].set_index('gene_id')['collapse_rate']

    collapse_rate_0_5 = pivot_0_5.reindex(TF_GENES).fillna(0).values
    collapse_rate_0_1 = pivot_0_1.reindex(TF_GENES).fillna(0).values

    run_007_response_file = '/home/zhanghl/projects/ddc_github/runs/run_007/data/response_metrics/response_summary.tsv'
    collapse_rate_KO = np.zeros(len(TF_GENES))
    if os.path.exists(run_007_response_file):
        run_007_df = pd.read_csv(run_007_response_file, sep='\t')
        filtered_df = run_007_df[run_007_df['seed'].isin(STEADY_SEEDS)]
        tf_filtered = filtered_df[filtered_df['gene_type'] == 'TF']
        for gene_idx in TF_GENES:
            gene_data = tf_filtered[tf_filtered['gene_id'] == gene_idx]
            n_collapse = (gene_data['ko_regime'] == 'collapse-like').sum()
            n_total = len(gene_data)
            collapse_rate_KO[gene_idx] = n_collapse / n_total if n_total > 0 else 0.0

    classifications = []
    for i, gene_idx in enumerate(TF_GENES):
        cr_0_5 = collapse_rate_0_5[i]
        cr_0_1 = collapse_rate_0_1[i]
        cr_KO = collapse_rate_KO[i]

        if cr_KO == 0:
            classification = 'non_essential'
        elif cr_0_1 == 0 and cr_0_5 == 0:
            classification = 'robust'
        elif cr_0_5 > 0:
            classification = 'haploinsufficient_like'
        elif cr_0_5 == 0 and cr_0_1 > 0:
            classification = 'intermediate_dosage_sensitive'
        elif cr_KO > 0 and cr_0_1 == 0:
            classification = 'KO_only_sensitive'
        else:
            classification = 'unknown'

        classifications.append({
            'gene_id': gene_idx,
            'collapse_rate_0.5': cr_0_5,
            'collapse_rate_0.1': cr_0_1,
            'collapse_rate_KO': cr_KO,
            'classification': classification
        })

    classification_df = pd.DataFrame(classifications)

    output_file = RESULTS_DIR / 'aggregated' / 'TF_dosage_sensitivity.tsv'
    classification_df.to_csv(output_file, sep='\t', index=False)
    print(f"\nSaved: {output_file}")

    print("\nClassification Results:")
    for _, row in classification_df.iterrows():
        print(f"  Gene {row['gene_id']}: λ=0.5 ({row['collapse_rate_0.5']:.2f}), "
              f"λ=0.1 ({row['collapse_rate_0.1']:.2f}), KO ({row['collapse_rate_KO']:.2f}) "
              f"→ {row['classification']}")

    return classification_df, gene_results_df


def generate_plots(classification_df: pd.DataFrame, gene_results_df: pd.DataFrame):
    print("\n" + "=" * 60)
    print("Generating Plots")
    print("=" * 60)

    generate_collapse_rate_plots(classification_df)
    generate_trajectory_comparison_plots()


def generate_collapse_rate_plots(classification_df: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax = axes[0]
    lambda_labels = ['λ=0.5', 'λ=0.1', 'KO']
    avg_collapse_rates = [
        classification_df['collapse_rate_0.5'].mean(),
        classification_df['collapse_rate_0.1'].mean(),
        classification_df['collapse_rate_KO'].mean()
    ]
    x = np.arange(len(lambda_labels))
    bars = ax.bar(x, avg_collapse_rates, color=['#3498DB', '#E74C3C', '#2C3E50'], alpha=0.8)
    labels = [f'{float(Decimal(str(v)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP)):.2f}' for v in avg_collapse_rates]
    ax.bar_label(bars, labels=labels, fontsize=10, padding=3)
    ax.set_xlabel('λ (Dosage Factor)', fontsize=12)
    ax.set_ylabel('Average Collapse Rate', fontsize=12)
    ax.set_title('Collapse Response to TF Attenuation', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(lambda_labels)
    ax.set_ylim(0, 1.1)
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    ax.grid(axis='y', alpha=0.3)
    
    ax = axes[1]
    # 1. 将数据转换为 seaborn 喜欢的“长格式 (Long-form)”
    plot_data = []
    for gene_idx in TF_GENES:
        row = classification_df[classification_df['gene_id'] == gene_idx].iloc[0]
        plot_data.extend([
            {'TF': f'TF {gene_idx}', 'λ': 'Baseline', 'Collapse Rate': 0.0},
            {'TF': f'TF {gene_idx}', 'λ': '0.5', 'Collapse Rate': row['collapse_rate_0.5']},
            {'TF': f'TF {gene_idx}', 'λ': '0.1', 'Collapse Rate': row['collapse_rate_0.1']},
            {'TF': f'TF {gene_idx}', 'λ': 'KO', 'Collapse Rate': row['collapse_rate_KO']}
        ])
    df_plot = pd.DataFrame(plot_data)

    # 2. 使用 pointplot 配合 dodge 参数 (解决离散重叠问题)
    # dodge=0.15 会让同X轴的点稍微并排错开，互不遮挡
    sns.pointplot(data=df_plot, x='λ', y='Collapse Rate', hue='TF',
                  dodge=0.15, marker='o', 
                  linewidth=1.5,      # 把线调细一点 (默认较粗)
                  markersize=4,     # 把点调小一点
                  ax=ax)

    ax.set_title('per-TF Sensitivity Trajectories', fontsize=14)
    ax.set_xlabel('λ (Dosage Factor)', fontsize=12)
    ax.set_ylabel('Collapse Rate', fontsize=12)
    ax.set_ylim(-0.05, 1.05)
    ax.set_yticks(np.arange(0, 1.1, 0.25))
    ax.grid(True, alpha=0.3)
    
    # 优化图例位置，防止挡住线条
    ax.legend(title='TF ID', loc='upper left', fontsize=9)

    plt.tight_layout()
    plot_file = PLOTS_DIR / 'collapse_rate_vs_lambda.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_file}")


def generate_trajectory_comparison_plots():
    from matplotlib.lines import Line2D

    print("\n  Generating trajectory comparison plots...")

    run_007_traj_dir = Path('/home/zhanghl/projects/ddc_github/runs/run_007/data/trajectories')
    baseline_traj_dir = DATA_DIR / 'trajectories'

    TF_GENES_LIST = list(TF_GENES)

    for seed in STEADY_SEEDS:
        seed_traj_dir = baseline_traj_dir / f'seed_{seed}'
        run_007_seed_traj_dir = run_007_traj_dir / f'seed_{seed}'

        if not seed_traj_dir.exists():
            print(f"    Warning: Seed trajectory directory not found: {seed_traj_dir}")
            continue

        baseline_data = load_trajectory_tsv(str(seed_traj_dir), 'baseline')
        P_baseline = baseline_data['P_traj']
        T, G = P_baseline.shape

        for lambda_val in LAMBDA_VALUES:
            fig, axes = plt.subplots(3, 3, figsize=(15, 12))
            fig.suptitle(f'Seed {seed}, λ={lambda_val}: Trajectory Comparison', fontsize=14, fontweight='bold', y=0.98)

            time_axis_plot = np.arange(1, T)

            def plot_trajectory(P_data, time_axis, gene_idx_to_highlight=None, ax=None):
                if ax is None:
                    ax = axes[0, 0]

                y_max = P_data.max()
                epi_gene_set = set(EPI_GENES)
                tf_gene_set = set(TF_GENES)

                def get_gene_style(g, highlight_idx):
                    if highlight_idx is not None and g == highlight_idx:
                        return 'black', 3.0, 1.0
                    elif g in tf_gene_set:
                        return plt.cm.Set1(0), 1.2, 0.7
                    elif g in epi_gene_set:
                        return plt.cm.Set1(3), 1.2, 0.7
                    else:
                        return 'gray', 0.6, 0.4

                for g in range(P_data.shape[1]):
                    color, lw, alpha = get_gene_style(g, gene_idx_to_highlight)
                    ax.plot(time_axis, P_data[:, g], '-', color=color, alpha=alpha, linewidth=lw)

                ax.set_ylim(0, y_max * 1.05)
                ax.set_xlim(0, T - 1)
                return ax

            gene_legend_baseline = [
                Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                Line2D([0], [0], color='gray', label='Others')
            ]

            axes[0, 0] = plot_trajectory(P_baseline[1:], time_axis_plot, ax=axes[0, 0])
            axes[0, 0].legend(handles=gene_legend_baseline, loc='upper right', fontsize=8)
            axes[0, 0].set_title('Baseline', fontsize=12, fontweight='bold')
            axes[0, 0].set_xlabel('Time')
            axes[0, 0].set_ylabel('Protein')
            axes[0, 0].grid(True, alpha=0.3)

            axes[0, 1].axis('off')
            axes[0, 2].axis('off')

            plot_idx = 0
            for gene_idx in TF_GENES_LIST:
                row = 1 + plot_idx // 3
                col = plot_idx % 3

                lambda_data = load_trajectory_tsv(str(seed_traj_dir), f'gene_{gene_idx}_lambda_{lambda_val}')
                P_lambda = lambda_data['P_traj']

                gene_legend_per_plot = [
                    Line2D([0], [0], color='black', linewidth=3, label=f'TF {gene_idx}'),
                    Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                    Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                    Line2D([0], [0], color='gray', label='Others')
                ]

                ax = axes[row, col]
                ax = plot_trajectory(P_lambda[1:], time_axis_plot, gene_idx_to_highlight=gene_idx, ax=ax)

                ax.legend(handles=gene_legend_per_plot, loc='upper right', fontsize=8, ncol=2)
                ax.set_title(f'Gene {gene_idx} (TF): λ={lambda_val}', fontsize=12)
                ax.set_xlabel('Time')
                ax.set_ylabel('Protein')
                ax.grid(True, alpha=0.3)

                plot_idx += 1

            for remaining in range(plot_idx, 6):
                row = 1 + remaining // 3
                col = remaining % 3
                axes[row, col].axis('off')

            plt.tight_layout()
            plot_file = PLOTS_DIR / 'trajectories' / f'P_traj_seed_{seed}_lambda_{lambda_val}.png'
            plot_file.parent.mkdir(parents=True, exist_ok=True)
            plt.savefig(plot_file, dpi=150, bbox_inches='tight')
            plt.close()
            print(f"    Saved: {plot_file}")


def main():
    results_df = run_simulation()
    classification_df, gene_results_df = analyze_results(results_df)
    generate_plots(classification_df, gene_results_df)

    print("\n" + "=" * 60)
    print("Run 008 completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()
