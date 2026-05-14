"""
Run 007 - KO-induced Collapse Analysis
=======================================

Objective: Analyze the long-term dynamics of steady systems under single-gene KO (rho_i = 0)
T = 100

Core Questions:
    - Which genes' KO leads to steady -> collapse-like transition?
    - Are there systematic differences between TF and epigenetic genes?

Perturbation Design:
    - KO (knockout): rho_i = 0 via apply_perturbation
    - Target genes: all TF genes + all epigenetic modifier genes
    - Baseline: steady worlds from run_005

Key Metrics:
    - collapse_rate per gene
    - delta_mean_P_terminal
    - deviation_score
    - transition_type (steady->steady / steady->collapse)

Author: zhanghl
Date: 2026-04-10
"""

import os
import sys
import copy
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
from typing import Dict, List, Tuple, Any

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')

import ddc
from ddc import TF_GENES, EPI_GENES, World, sample_world, sample_initial_state, simulate_single_cell

warnings.filterwarnings('ignore')

BASE_DIR = '/home/zhanghl/projects/ddc_github/runs/run_007_100steps'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
PLOTS_DIR = os.path.join(BASE_DIR, 'plots')

RUN_005_DIR = '/home/zhanghl/projects/ddc_github/runs/run_005'
RUN_005_PROXY_METRICS = os.path.join(RUN_005_DIR, 'results/gain_proxy_analysis/proxy_metrics.tsv')

T = 100
TERMINAL_WINDOW = 20
N_ACTIVE_THRESHOLD = 1
THRESHOLD = 0.1

TARGET_GENES = TF_GENES + EPI_GENES

COLORS = {
    'TF': '#E74C3C',
    'EPI': '#3498DB',
    'collapse': '#E74C3C',
    'steady': '#27AE60'
}


def setup_directories():
    for d in [DATA_DIR, RESULTS_DIR, PLOTS_DIR]:
        os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join(DATA_DIR, 'trajectories'), exist_ok=True)
    os.makedirs(os.path.join(DATA_DIR, 'response_metrics'), exist_ok=True)
    os.makedirs(os.path.join(DATA_DIR, 'regime_results'), exist_ok=True)
    os.makedirs(os.path.join(RESULTS_DIR, 'per_seed'), exist_ok=True)
    os.makedirs(os.path.join(RESULTS_DIR, 'aggregated'), exist_ok=True)
    os.makedirs(os.path.join(PLOTS_DIR, 'trajectories'), exist_ok=True)
    os.makedirs(os.path.join(PLOTS_DIR, 'regime_transition'), exist_ok=True)
    os.makedirs(os.path.join(PLOTS_DIR, 'class_comparison'), exist_ok=True)


def save_trajectory_tsv(traj: Dict[str, torch.Tensor], seed_dir: str, prefix: str):
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


def load_steady_seeds(proxy_metrics_path: str, n_worlds: int = 4) -> List[int]:
    df = pd.read_csv(proxy_metrics_path, sep='\t')
    steady_df = df[df['regime'] == 'steady'].head(n_worlds)
    return steady_df['seed'].tolist()


def classify_regime(X_traj: torch.Tensor) -> str:
    X_terminal = X_traj[-1, :]
    n_active = (X_terminal > THRESHOLD).sum().item()
    return 'steady' if n_active > 1 else 'collapse-like'


def simulate_with_ko(world_seed: int, X0: torch.Tensor, P0: torch.Tensor, Z0: torch.Tensor, N0: float,
                     gene_idx: int) -> Dict[str, torch.Tensor]:
    world = sample_world(world_seed)
    state = (X0.clone(), P0.clone(), Z0.clone(), N0)
    world_ko, _ = ddc.apply_perturbation(world, state, {'knockout': [gene_idx]})
    
    traj = simulate_single_cell(world_ko, X0.clone(), P0.clone(), Z0.clone(), N0, T)
    return traj


def compute_response_metrics(P_baseline: torch.Tensor, P_ko: torch.Tensor) -> Dict[str, float]:
    terminal_mean_P_baseline = P_baseline[-TERMINAL_WINDOW:].mean().item()
    terminal_mean_P_ko = P_ko[-TERMINAL_WINDOW:].mean().item()
    
    delta_mean_P_terminal = terminal_mean_P_ko - terminal_mean_P_baseline
    
    deviation = torch.abs(P_ko - P_baseline).sum(dim=1)
    mean_deviation = deviation.mean().item()
    
    return {
        'terminal_mean_P_baseline': terminal_mean_P_baseline,
        'terminal_mean_P_ko': terminal_mean_P_ko,
        'delta_mean_P_terminal': delta_mean_P_terminal,
        'mean_deviation': mean_deviation
    }


def run_analysis():
    print("=" * 60)
    print("Run 007 - KO-induced Collapse Analysis")
    print("=" * 60)
    
    setup_directories()
    
    print("\n[1/5] Loading steady seeds from run_005...")
    steady_seeds = load_steady_seeds(RUN_005_PROXY_METRICS, n_worlds=4)
    print(f"    Loaded {len(steady_seeds)} steady seeds: {steady_seeds}")
    
    results = []
    
    print("\n[2/5] Running baseline and KO simulations...")
    for seed in steady_seeds:
        print(f"\n    Processing seed {seed}...")
        
        world = sample_world(seed)
        cell_seed = seed + 1
        X0, P0, Z0, N0 = sample_initial_state(cell_seed, world)
        
        seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
        os.makedirs(seed_dir, exist_ok=True)
        
        print(f"        Running baseline trajectory...")
        traj_baseline = simulate_single_cell(world, X0.clone(), P0.clone(), Z0.clone(), N0, T)
        
        baseline_regime = classify_regime(traj_baseline['X_traj'])
        print(f"        Baseline regime: {baseline_regime}")

        save_trajectory_tsv(traj_baseline, seed_dir, 'baseline')

        for gene_idx in TARGET_GENES:
            gene_type = 'TF' if gene_idx in TF_GENES else 'EPI'
            print(f"        KO gene {gene_idx} ({gene_type})...")
            
            traj_ko = simulate_with_ko(seed, X0, P0, Z0, N0, gene_idx)
            
            ko_regime = classify_regime(traj_ko['X_traj'])
            transition_type = f"{baseline_regime}->{ko_regime}"
            
            metrics = compute_response_metrics(traj_baseline['P_traj'], traj_ko['P_traj'])
            
            result = {
                'seed': seed,
                'gene_id': gene_idx,
                'gene_type': gene_type,
                'baseline_regime': baseline_regime,
                'ko_regime': ko_regime,
                'transition_type': transition_type,
                **metrics
            }
            results.append(result)

            save_trajectory_tsv(traj_ko, seed_dir, f'KO_gene_{gene_idx}')

    print("\n[3/5] Saving response metrics...")
    results_df = pd.DataFrame(results)

    per_seed_dir = os.path.join(RESULTS_DIR, 'per_seed')
    for seed in steady_seeds:
        seed_df = results_df[results_df['seed'] == seed]
        seed_df = seed_df[['seed', 'gene_id', 'gene_type', 'baseline_regime', 'ko_regime',
                           'transition_type', 'delta_mean_P_terminal', 'mean_deviation']]
        seed_df = seed_df.rename(columns={'seed': 'seed_id', 'mean_deviation': 'deviation_score', 'delta_mean_P_terminal': 'Δmean_P_terminal'})
        seed_df.to_csv(os.path.join(per_seed_dir, f'seed_{seed}_results.tsv'), sep='\t', index=False)
    print(f"    Saved per-seed results to {per_seed_dir}/")

    response_metrics_path = os.path.join(DATA_DIR, 'response_metrics', 'response_summary.tsv')
    response_df = results_df.rename(columns={'mean_deviation': 'deviation_score', 'delta_mean_P_terminal': 'Δmean_P_terminal'})
    response_df = response_df[['seed', 'gene_id', 'gene_type', 'terminal_mean_P_baseline', 'terminal_mean_P_ko', 'Δmean_P_terminal', 'deviation_score', 'baseline_regime', 'ko_regime', 'transition_type']]
    response_df.to_csv(response_metrics_path, sep='\t', index=False)
    print(f"    Saved to {response_metrics_path}")
    
    print("\n[4/5] Computing aggregated statistics...")
    gene_summary = results_df.groupby(['gene_id', 'gene_type']).agg({
        'transition_type': lambda x: (x == 'steady->collapse-like').sum(),
        'delta_mean_P_terminal': 'mean',
        'mean_deviation': 'mean'
    }).reset_index()
    gene_summary.columns = ['gene_id', 'gene_type', 'collapse_count', 'mean_Δmean_P_terminal', 'mean_deviation_score']
    gene_summary['collapse_rate'] = gene_summary['collapse_count'] / len(steady_seeds)
    gene_summary = gene_summary[['gene_id', 'gene_type', 'collapse_count', 'collapse_rate', 'mean_Δmean_P_terminal', 'mean_deviation_score']]
    
    collapse_rate_path = os.path.join(RESULTS_DIR, 'aggregated', 'gene_level_summary.tsv')
    gene_summary.to_csv(collapse_rate_path, sep='\t', index=False)
    
    class_summary = results_df.groupby('gene_type').agg({
        'transition_type': lambda x: (x == 'steady->collapse-like').sum(),
        'delta_mean_P_terminal': 'mean',
        'mean_deviation': 'mean'
    }).reset_index()
    class_summary.columns = ['gene_type', 'total_collapse_count', 'average_Δmean_P_terminal', 'average_deviation_score']
    class_summary['n_genes'] = class_summary['gene_type'].map(
        {'TF': len(TF_GENES), 'EPI': len(EPI_GENES)}
    )
    class_summary['average_collapse_rate'] = class_summary['total_collapse_count'] / (len(steady_seeds) * class_summary['n_genes'])
    class_summary = class_summary[['gene_type', 'average_collapse_rate', 'average_Δmean_P_terminal', 'average_deviation_score']]
    
    class_summary_path = os.path.join(RESULTS_DIR, 'aggregated', 'class_level_summary.tsv')
    class_summary.to_csv(class_summary_path, sep='\t', index=False)
    
    print(f"\n    Gene-level summary:")
    print(gene_summary.to_string(index=False))
    print(f"\n    Class-level summary:")
    print(class_summary.to_string(index=False))
    
    print("\n[5/5] Generating visualizations...")
    generate_plots(results_df, gene_summary, class_summary, steady_seeds)
    
    print("\n" + "=" * 60)
    print("Run 007 completed successfully!")
    print("=" * 60)
    
    return results_df, gene_summary, class_summary


def generate_plots(results_df: pd.DataFrame, gene_summary: pd.DataFrame, class_summary: pd.DataFrame, steady_seeds: List[int]):
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    ax1 = axes[0, 0]
    for gene_type in ['TF', 'EPI']:
        subset = gene_summary[gene_summary['gene_type'] == gene_type]
        color = COLORS.get(gene_type, 'gray')
        ax1.bar(subset['gene_id'].astype(str), subset['collapse_rate'], 
                 color=color, alpha=0.7, label=gene_type)
    ax1.set_xlabel('Gene ID')
    ax1.set_ylabel('Collapse Rate')
    ax1.set_title('Collapse Rate per Gene')
    ax1.legend()
    ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    
    ax2 = axes[0, 1]
    x = np.arange(len(gene_summary))
    width = 0.35
    ax2.bar(x - width/2, gene_summary['mean_Δmean_P_terminal'], width, label='Δmean_P', color='steelblue')
    ax2.set_xlabel('Gene ID')
    ax2.set_ylabel('Δmean_P Terminal')
    ax2.set_title('Terminal Expression Shift (KO vs Baseline)')
    ax2.set_xticks(x)
    ax2.set_xticklabels(gene_summary['gene_id'].astype(str), rotation=45)
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
    
    ax3 = axes[1, 0]
    transition_counts = results_df.groupby(['gene_type', 'transition_type']).size().unstack(fill_value=0)
    transition_counts.plot(kind='bar', ax=ax3, color=[COLORS['collapse'], COLORS['steady']])
    ax3.set_xlabel('Gene Type')
    ax3.set_ylabel('Count')
    ax3.set_title('Regime Transition Distribution by Gene Type')
    ax3.legend(title='Transition')
    ax3.set_xticklabels(ax3.get_xticklabels(), rotation=0)
    
    ax4 = axes[1, 1]
    class_summary_plot = class_summary.set_index('gene_type')
    class_summary_plot[['average_collapse_rate']].plot(kind='bar', ax=ax4, color=[COLORS['TF'], COLORS['EPI']])
    ax4.set_xlabel('Gene Class')
    ax4.set_ylabel('Collapse Rate')
    ax4.set_title('Collapse Rate: TF vs Epigenetic Genes')
    ax4.legend().remove()
    ax4.set_xticklabels(ax4.get_xticklabels(), rotation=0)
    ax4.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    comparison_path = os.path.join(PLOTS_DIR, 'class_comparison', 'TF_vs_EPI_comparison.png')
    plt.savefig(comparison_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {comparison_path}")
    
    fig2, axes2 = plt.subplots(1, 2, figsize=(12, 5))
    
    ax5 = axes2[0]
    transition_types = results_df['transition_type'].value_counts()
    colors_transition = [COLORS['collapse'] if 'collapse' in t else COLORS['steady'] for t in transition_types.index]
    ax5.pie(transition_types.values, labels=transition_types.index, autopct='%1.1f%%', colors=colors_transition)
    ax5.set_title('Overall Transition Type Distribution')
    
    ax6 = axes2[1]
    gene_summary_sorted = gene_summary.sort_values('collapse_rate', ascending=True)
    colors_genes = [COLORS.get(g, 'gray') for g in gene_summary_sorted['gene_type']]
    ax6.barh(gene_summary_sorted['gene_id'].astype(str), gene_summary_sorted['collapse_rate'], color=colors_genes)
    ax6.set_xlabel('Collapse Rate')
    ax6.set_ylabel('Gene ID')
    ax6.set_title('Gene Collapse Rate Ranking')
    ax6.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    summary_path = os.path.join(PLOTS_DIR, 'regime_transition', 'transition_summary.png')
    plt.savefig(summary_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {summary_path}")

    print("\n    Generating trajectory comparison plots...")
    for seed in steady_seeds[:2]:
        generate_trajectory_comparison_plot(seed, results_df)


def generate_trajectory_comparison_plot(seed: int, results_df: pd.DataFrame):
    """
    为单个seed生成轨迹对比图（4x3布局）

    子图布局：
        (0,0): Baseline轨迹（run_007基准系统）
        (0,1): Run_006扰动轨迹（含扰动类型、目标基因、regime转换信息）
        (0,2): 空置
        (1-3行): 9个KO轨迹，分别展示各基因被KO后的系统动态

    Args:
        seed: 当前seed值
        results_df: 包含所有seed结果的DataFrame，用于获取transition_type信息
    """
    seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
    baseline_traj = load_trajectory_tsv(seed_dir, 'baseline')
    P_baseline = baseline_traj['P_traj']
    T, G = P_baseline.shape

    seed_results = results_df[results_df['seed'] == seed]
    if len(seed_results) > 0:
        baseline_regime = seed_results['baseline_regime'].values[0]
    else:
        baseline_regime = 'unknown'

    # 从run_006读取同一seed的扰动信息（扰动类型、目标基因、recovery状态）
    run_006_response_path = '/home/zhanghl/projects/ddc_github/runs/run_006/data/response_metrics/response_summary.tsv'
    run_006_seed_info = None
    if os.path.exists(run_006_response_path):
        df_run006_resp = pd.read_csv(run_006_response_path, sep='\t')
        run_006_seed_info = df_run006_resp[df_run006_resp['seed'] == seed]
        if len(run_006_seed_info) > 0:
            run_006_seed_info = run_006_seed_info.iloc[0]

    # 加载run_006的轨迹数据（用于子图0,1对比）
    run_006_traj_path = '/home/zhanghl/projects/ddc_github/runs/run_006/data/trajectories/perturbed.tsv'
    P_run006 = None
    if os.path.exists(run_006_traj_path):
        df_run006 = pd.read_csv(run_006_traj_path, sep='\t')
        df_run006_seed = df_run006[df_run006['seed'] == seed]
        if len(df_run006_seed) > 0:
            df_wide = df_run006_seed.pivot(index='time', columns='gene_id', values='P')
            P_run006 = df_wide.values
            if P_run006.shape[0] > T:
                P_run006 = P_run006[:T]

    # 绘制起始时间从t=1开始，但x轴范围保持从0开始以显示完整时间尺度
    time_axis_plot = np.arange(1, T)
    tf_gene_set = set(TF_GENES)
    epi_gene_set = set(EPI_GENES)

    from matplotlib.lines import Line2D

    fig = plt.figure(figsize=(18, 16))
    gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.3, wspace=0.2)

    Y_MAX_NORMAL = 0.35

    def plot_trajectory_with_broken_axis(P_data, time_axis, tf_gene_set, epi_gene_set, gene_idx_to_highlight=None, row=None, col=None):
        y_max = P_data.max()

        def get_gene_style(g, highlight_idx):
            if highlight_idx is not None and g == highlight_idx:
                if g in epi_gene_set:
                    return plt.cm.Set1(3), 3.0, 1.0
                else:
                    return '#FF0000', 3.0, 1.0
            elif g in tf_gene_set:
                return plt.cm.Set1(0), 1.2, 0.7
            elif g in epi_gene_set:
                return plt.cm.Set1(3), 1.2, 0.7
            else:
                return 'gray', 0.6, 0.4

        def draw_curves(ax, highlight_idx):
            for g in range(P_data.shape[1]):
                color, lw, alpha = get_gene_style(g, highlight_idx)
                ax.plot(time_axis, P_data[:, g], '-', color=color, alpha=alpha, linewidth=lw)

        if row is not None:
            ax = fig.add_subplot(gs[row, col])
        else:
            ax = fig.add_subplot(111)

        draw_curves(ax, gene_idx_to_highlight)
        ax.set_ylim(0, y_max * 1.05)
        return ax

    # ============== 子图 (0,0): Baseline轨迹 ==============
    ax = plot_trajectory_with_broken_axis(P_baseline[1:], time_axis_plot, tf_gene_set, epi_gene_set, row=0, col=0)
    gene_legend = [Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                   Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                   Line2D([0], [0], color='gray', label='Others')]
    ax.legend(handles=gene_legend, loc='upper right', fontsize=8)
    ax.set_title(f'Baseline: {baseline_regime}', fontsize=12, fontweight='bold')
    ax.set_xlabel('Time')
    ax.set_ylabel('Protein')
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0, T - 1)

    # ============== 子图 (0,1): Run_006扰动轨迹 ==============
    if P_run006 is not None:
        ax = plot_trajectory_with_broken_axis(P_run006[1:], time_axis_plot, tf_gene_set, epi_gene_set, row=0, col=1)
        gene_legend_run006 = [Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                              Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                              Line2D([0], [0], color='gray', label='Others')]
        ax.legend(handles=gene_legend_run006, loc='upper right', fontsize=8)
        if run_006_seed_info is not None:
            perturb_type = run_006_seed_info['perturbation_type']
            target_gene = int(run_006_seed_info['target_gene'])
            target_gene_type = 'TF' if target_gene in TF_GENES else 'EPI'
            perturb_type_full = 'OE' if perturb_type == 'OE' else 'KD'
            recovery = run_006_seed_info.get('recovery', 1)
            transition_run006 = 'steady->steady' if recovery >= 1 else 'steady->collapse-like'
            title_color = '#FF0000' if transition_run006 == 'steady->collapse-like' else 'black'
            ax.set_title(f'{perturb_type_full} Gene {target_gene} ({target_gene_type}): {transition_run006}', fontsize=12, fontweight='bold', color=title_color)
        else:
            ax.set_title('Perturbed (run_006)', fontsize=12, fontweight='bold')
        ax.set_xlabel('Time')
        ax.set_ylabel('Protein')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, T - 1)

    # ============== 子图 (0,2): 空置 ==============
    ax_empty = fig.add_subplot(gs[0, 2])
    ax_empty.axis('off')

    # ============== 子图 (1-3行): KO轨迹（共9个基因）==============
    for row in range(1, 4):
        for col in range(3):
            idx = (row - 1) * 3 + col
            if idx < len(TARGET_GENES):
                gene_idx = TARGET_GENES[idx]
                gene_type = 'TF' if gene_idx in TF_GENES else 'EPI'

                gene_result = seed_results[seed_results['gene_id'] == gene_idx]
                if len(gene_result) > 0:
                    transition_type = gene_result['transition_type'].values[0]
                else:
                    transition_type = 'unknown'

                ko_traj = load_trajectory_tsv(seed_dir, f'KO_gene_{gene_idx}')
                P_ko = ko_traj['P_traj']

                ax = plot_trajectory_with_broken_axis(P_ko[1:], time_axis_plot, tf_gene_set, epi_gene_set, gene_idx_to_highlight=gene_idx, row=row, col=col)

                if gene_idx in epi_gene_set:
                    gene_legend = [Line2D([0], [0], color=plt.cm.Set1(3), linewidth=3, label=f'KO Gene {gene_idx}'),
                                   Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                                   Line2D([0], [0], color='gray', label='Others'),
                                   Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics')]
                else:
                    gene_legend = [Line2D([0], [0], color='#FF0000', linewidth=3, label=f'KO Gene {gene_idx}'),
                                   Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                                   Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                                   Line2D([0], [0], color='gray', label='Others')]
                ax.legend(handles=gene_legend, loc='upper right', fontsize=8, ncol=2)
                title_color = '#FF0000' if transition_type == 'steady->collapse-like' else 'black'
                ax.set_title(f'KO Gene {gene_idx} ({gene_type}): {transition_type}', fontsize=12, fontweight='bold', color=title_color)
                ax.set_xlabel('Time')
                ax.set_ylabel('Protein')
                ax.grid(True, alpha=0.3)
                ax.set_xlim(0, T - 1)
            else:
                ax_empty = fig.add_subplot(gs[row, col])
                ax_empty.axis('off')

    plt.suptitle(f'Seed {seed}: Trajectory Comparison', fontsize=16, fontweight='bold', y=0.92)

    traj_path = os.path.join(PLOTS_DIR, 'trajectories', f'traj_comparison_seed_{seed}.png')
    plt.savefig(traj_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {traj_path}")


if __name__ == '__main__':
    results_df, gene_summary, class_summary = run_analysis()