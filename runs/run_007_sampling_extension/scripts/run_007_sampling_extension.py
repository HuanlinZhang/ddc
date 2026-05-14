"""
Run 007 - Sampling Extension: KO-induced Collapse Analysis
=========================================================

Based on run_004/03_sampling_extension and run_005 steady worlds.

Objective: Analyze the long-term dynamics of steady systems under single-gene KO (rho_i = 0)
T = 400

Core Questions:
    - Which genes' KO leads to steady -> collapse-like transition?
    - Are there systematic differences between TF and epigenetic genes?

Perturbation Design:
    - KO (knockout): rho_i = 0 via apply_perturbation
    - Target genes: all TF genes + all epigenetic modifier genes (9 genes)
    - Baseline: steady worlds from run_004/03_sampling_extension and run_005

Key Metrics:
    - collapse_rate per gene
    - delta_mean_P_terminal
    - deviation_score
    - transition_type (steady->steady / steady->collapse-like)

Author: zhanghl
Date: 2026-04-16
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
import networkx as nx
from typing import Dict, List, Tuple, Any

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')

import ddc
from ddc import TF_GENES, EPI_GENES, World, sample_world, sample_initial_state, simulate_single_cell

warnings.filterwarnings('ignore')

BASE_DIR = '/home/zhanghl/projects/ddc_github/runs/run_007_sampling_extension'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
PLOTS_DIR = os.path.join(BASE_DIR, 'plots')

RUN_004_DATA = '/home/zhanghl/projects/ddc_github/runs/run_004/03_sampling_extension/data'
RUN_005_PROXY_METRICS = '/home/zhanghl/projects/ddc_github/runs/run_005/results/resource_ablation/regime_comparison.tsv'

T = 400
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


def load_steady_seeds_stratified(n_per_group: int = 4) -> pd.DataFrame:
    run005_df = pd.read_csv(RUN_005_PROXY_METRICS, sep='\t')
    run005_steady = run005_df[run005_df['original_regime'] == 'steady'][['seed', 'N_active_TF_original']].drop_duplicates()
    run005_steady = run005_steady.rename(columns={'N_active_TF_original': 'N_active_TF'})
    run005_steady['source'] = 'run_005'

    run004_df = pd.read_csv(os.path.join(RUN_004_DATA, 'world_regime_distribution.tsv'), sep='\t')
    run004_steady = run004_df[run004_df['regime'] == 'steady'][['seed', 'N_active_TF']].drop_duplicates()
    run004_steady['source'] = 'run_004_03'

    all_steady = pd.concat([run005_steady, run004_steady], ignore_index=True)

    def get_tf_group(n_active: int) -> str:
        if n_active <= 2:
            return 'low_tf'
        elif n_active <= 4:
            return 'mid_tf'
        else:
            return 'high_tf'

    all_steady['tf_group'] = all_steady['N_active_TF'].apply(get_tf_group)

    np.random.seed(42)
    low_tf = all_steady[all_steady['tf_group'] == 'low_tf']['seed'].tolist()
    mid_tf = all_steady[all_steady['tf_group'] == 'mid_tf']['seed'].tolist()
    high_tf = all_steady[all_steady['tf_group'] == 'high_tf']['seed'].tolist()

    np.random.shuffle(low_tf)
    np.random.shuffle(mid_tf)
    np.random.shuffle(high_tf)

    selected = []
    selected.extend(low_tf[:n_per_group] if len(low_tf) >= n_per_group else low_tf)
    selected.extend(mid_tf[:n_per_group] if len(mid_tf) >= n_per_group else mid_tf)
    selected.extend(high_tf[:n_per_group] if len(high_tf) >= n_per_group else high_tf)

    selected_info = all_steady[all_steady['seed'].isin(selected)][['seed', 'N_active_TF', 'tf_group', 'source']].copy()
    selected_info = selected_info.sort_values('seed').reset_index(drop=True)

    return selected_info


def classify_regime(X_traj: torch.Tensor) -> str:
    X_terminal = X_traj[-1, :]
    n_active = (X_terminal > THRESHOLD).sum().item()
    return 'steady' if n_active > N_ACTIVE_THRESHOLD else 'collapse-like'


def simulate_with_ko(world: World, X0: torch.Tensor, P0: torch.Tensor, Z0: torch.Tensor, N0: float,
                     gene_idx: int) -> Dict[str, torch.Tensor]:
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


def compute_tf_network_metrics_for_world(world: World, seed: int) -> pd.DataFrame:
    """Compute TF network metrics (degree, SCC) for a given world"""
    a_ij = world.a_ij

    G = nx.DiGraph()
    G.add_nodes_from(TF_GENES)

    for source in TF_GENES:
        if source not in a_ij:
            continue
        for target in TF_GENES:
            if target in a_ij[source] and a_ij[source][target] > 0:
                G.add_edge(source, target, weight=a_ij[source][target])

    in_degrees = dict(G.in_degree())
    out_degrees = dict(G.out_degree())
    total_degrees = dict(G.degree())

    scc = list(nx.strongly_connected_components(G))
    largest_scc = max(scc, key=len) if scc else set()
    scc_id_map = {}
    for idx, component in enumerate(scc):
        for node in component:
            scc_id_map[node] = idx

    metrics = []
    for tf_node in TF_GENES:
        metrics.append({
            'seed': seed,
            'tf_node': tf_node,
            'in_degree': in_degrees.get(tf_node, 0),
            'out_degree': out_degrees.get(tf_node, 0),
            'total_degree': total_degrees.get(tf_node, 0),
            'SCC_id': scc_id_map.get(tf_node, -1),
            'in_largest_SCC': tf_node in largest_scc,
            'largest_SCC_size': len(largest_scc)
        })

    return pd.DataFrame(metrics)


def run_analysis():
    print("=" * 60)
    print("Run 007 - Sampling Extension: KO-induced Collapse Analysis")
    print("=" * 60)

    setup_directories()

    print("\n[1/5] Loading steady seeds (combined run_005 + run_004/03_sampling_extension)...")
    steady_seeds_info = load_steady_seeds_stratified(n_per_group=4)
    steady_seeds = steady_seeds_info['seed'].tolist()
    print(f"    Selected {len(steady_seeds)} steady seeds (stratified by N_active_TF):")
    print(steady_seeds_info.to_string(index=False))

    results = []

    print("\n[2/5] Running baseline and KO simulations...")
    for seed in steady_seeds:
        print(f"\n    Processing seed {seed}...")

        world = sample_world(seed)
        cell_seed = seed + 1
        X0, P0, Z0, N0 = sample_initial_state(cell_seed, world)

        seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
        os.makedirs(seed_dir, exist_ok=True)

        baseline_path = os.path.join(seed_dir, 'baseline_P_traj.tsv')
        if os.path.exists(baseline_path):
            print(f"        Loading cached baseline trajectory...")
            traj_baseline = load_trajectory_tsv(seed_dir, 'baseline')
            baseline_regime = classify_regime(torch.tensor(traj_baseline['X_traj']))
        else:
            print(f"        Running baseline trajectory...")
            traj_baseline = simulate_single_cell(world, X0.clone(), P0.clone(), Z0.clone(), N0, T)
            save_trajectory_tsv(traj_baseline, seed_dir, 'baseline')
            baseline_regime = classify_regime(traj_baseline['X_traj'])
        print(f"        Baseline regime: {baseline_regime}")

        seed_results = {
            'seed': seed,
            'baseline_regime': baseline_regime,
            'baseline_P': torch.tensor(traj_baseline['P_traj'])
        }

        for gene_idx in TARGET_GENES:
            gene_type = 'TF' if gene_idx in TF_GENES else 'EPI'
            ko_path = os.path.join(seed_dir, f'KO_gene_{gene_idx}_P_traj.tsv')

            if os.path.exists(ko_path):
                print(f"        Loading cached KO gene {gene_idx} ({gene_type})...")
                traj_ko = load_trajectory_tsv(seed_dir, f'KO_gene_{gene_idx}')
            else:
                print(f"        KO gene {gene_idx} ({gene_type})...")
                traj_ko = simulate_with_ko(world, X0, P0, Z0, N0, gene_idx)
                save_trajectory_tsv(traj_ko, seed_dir, f'KO_gene_{gene_idx}')

            traj_ko_t = torch.tensor(traj_ko['P_traj']) if isinstance(traj_ko, dict) else traj_ko
            ko_regime = classify_regime(torch.tensor(traj_ko['X_traj']) if isinstance(traj_ko, dict) else traj_ko['X_traj'])
            transition_type = f"{baseline_regime}->{ko_regime}"

            metrics = compute_response_metrics(seed_results['baseline_P'], traj_ko_t)

            results.append({
                'seed': seed,
                'gene_id': gene_idx,
                'gene_type': gene_type,
                'baseline_regime': baseline_regime,
                'ko_regime': ko_regime,
                'transition_type': transition_type,
                **metrics
            })

    results_df = pd.DataFrame(results)

    print("\n[3/5] Saving response metrics...")
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
        'transition_type': [
            lambda x: (x == 'steady->collapse-like').sum(),
            lambda x: (x == 'steady->steady').sum()
        ],
        'delta_mean_P_terminal': 'mean',
        'mean_deviation': 'mean'
    }).reset_index()
    gene_summary.columns = ['gene_id', 'gene_type', 'collapse_count', 'steady_count', 'mean_Δmean_P_terminal', 'mean_deviation_score']
    gene_summary['collapse_rate'] = gene_summary['collapse_count'] / (gene_summary['collapse_count'] + gene_summary['steady_count'])
    gene_summary = gene_summary[['gene_id', 'gene_type', 'collapse_count', 'collapse_rate', 'mean_Δmean_P_terminal', 'mean_deviation_score']]

    class_summary = results_df.groupby('gene_type').agg({
        'transition_type': [
            lambda x: (x == 'steady->collapse-like').sum(),
            lambda x: (x == 'steady->steady').sum()
        ],
        'delta_mean_P_terminal': 'mean',
        'mean_deviation': 'mean'
    }).reset_index()
    class_summary.columns = ['gene_type', 'total_collapse_count', 'total_steady_count', 'average_Δmean_P_terminal', 'average_deviation_score']
    class_summary['average_collapse_rate'] = class_summary['total_collapse_count'] / (class_summary['total_collapse_count'] + class_summary['total_steady_count'])
    class_summary = class_summary[['gene_type', 'average_collapse_rate', 'average_Δmean_P_terminal', 'average_deviation_score']]

    gene_summary_path = os.path.join(RESULTS_DIR, 'aggregated', 'gene_level_summary.tsv')
    gene_summary.to_csv(gene_summary_path, sep='\t', index=False)
    print(f"    Saved: {gene_summary_path}")

    class_summary_path = os.path.join(RESULTS_DIR, 'aggregated', 'class_level_summary.tsv')
    class_summary.to_csv(class_summary_path, sep='\t', index=False)
    print(f"    Saved: {class_summary_path}")

    transition_matrix = results_df.groupby(['gene_id', 'transition_type']).size().unstack(fill_value=0)
    transition_matrix['total'] = transition_matrix.sum(axis=1)
    transition_matrix_path = os.path.join(RESULTS_DIR, 'aggregated', 'regime_transition_matrix.tsv')
    transition_matrix.to_csv(transition_matrix_path, sep='\t')
    print(f"    Saved: {transition_matrix_path}")

    print(f"\n    Gene-level summary:")
    print(gene_summary.to_string(index=False))
    print(f"\n    Class-level summary:")
    print(class_summary.to_string(index=False))

    print("\n[5/5] Generating visualizations...")
    generate_visualizations(results_df, steady_seeds, steady_seeds_info, gene_summary, class_summary)

    print("\n[6/6] Running surviving TF analysis...")
    run_surviving_tf_analysis(gene_summary)

    print("\n[Extra] Running cross-run comparison...")
    run_cross_run_comparison(gene_summary)

    print("\n" + "=" * 60)
    print("Run 007 Sampling Extension completed successfully!")
    print("=" * 60)

    return results_df, gene_summary, class_summary


def run_cross_run_comparison(gene_summary: pd.DataFrame):
    print("\n" + "=" * 60)
    print("Cross-run Comparison (run_006 vs run_007)")
    print("=" * 60)

    run_006_path = '/home/zhanghl/projects/ddc_github/runs/run_006/data/response_metrics/response_summary.tsv'
    if not os.path.exists(run_006_path):
        print(f"    run_006 data not found at {run_006_path}, skipping cross-run comparison")
        return

    run_006_df = pd.read_csv(run_006_path, sep='\t')
    run_006_agg = run_006_df.groupby('target_gene').agg({
        'recovery': 'mean'
    }).reset_index()
    run_006_agg.columns = ['gene_id', 'recovery_score']

    cross_run = gene_summary[gene_summary['gene_type'] == 'TF'][['gene_id', 'collapse_rate', 'mean_Δmean_P_terminal']].merge(
        run_006_agg, on='gene_id', how='outer'
    )
    cross_run = cross_run[['gene_id', 'recovery_score', 'collapse_rate', 'mean_Δmean_P_terminal']]

    cross_run_dir = os.path.join(BASE_DIR, 'results', 'cross_run')
    os.makedirs(cross_run_dir, exist_ok=True)
    cross_run_path = os.path.join(cross_run_dir, 'run_006_vs_run_007_comparison.tsv')
    cross_run.to_csv(cross_run_path, sep='\t', index=False)
    print(f"    Saved: {cross_run_path}")


def run_surviving_tf_analysis(gene_summary: pd.DataFrame):
    print("\n" + "=" * 60)
    print("Surviving TF Analysis")
    print("=" * 60)

    surviving_tf_dir = os.path.join(BASE_DIR, 'data', 'surviving_TF')
    results_surviving_dir = os.path.join(BASE_DIR, 'results', 'surviving_TF')
    plots_surviving_dir = os.path.join(PLOTS_DIR, 'surviving_TF')
    os.makedirs(surviving_tf_dir, exist_ok=True)
    os.makedirs(results_surviving_dir, exist_ok=True)
    os.makedirs(plots_surviving_dir, exist_ok=True)

    tf_genes = gene_summary[gene_summary['gene_type'] == 'TF'].copy()
    tf_genes['group'] = tf_genes['collapse_rate'].apply(lambda x: 'surviving' if x < 1 else 'collapsing')

    surviving_list = tf_genes[tf_genes['group'] == 'surviving'][['gene_id', 'collapse_rate']].copy()
    surviving_list['gene_type'] = 'TF'
    surviving_list_path = os.path.join(surviving_tf_dir, 'surviving_TF_list.tsv')
    surviving_list.to_csv(surviving_list_path, sep='\t', index=False)
    print(f"    Saved: {surviving_list_path}")

    steady_seeds_info = load_steady_seeds_stratified(n_per_group=4)
    all_seed_metrics = []

    for _, row in steady_seeds_info.iterrows():
        seed = row['seed']
        world = sample_world(seed)
        seed_metrics = compute_tf_network_metrics_for_world(world, seed)
        all_seed_metrics.append(seed_metrics)

    combined_metrics = pd.concat(all_seed_metrics, ignore_index=True)

    aggregated_metrics = combined_metrics.groupby('tf_node').agg({
        'in_degree': 'mean',
        'out_degree': 'mean',
        'total_degree': 'mean',
        'SCC_id': lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else x.iloc[0],
        'largest_SCC_size': 'first'
    }).reset_index()
    aggregated_metrics.columns = ['gene_id', 'in_degree', 'out_degree', 'total_degree', 'SCC_id', 'largest_SCC_size']

    network_metrics_path = os.path.join(surviving_tf_dir, 'surviving_TF_network_metrics.tsv')
    aggregated_metrics.to_csv(network_metrics_path, sep='\t', index=False)
    print(f"    Saved: {network_metrics_path}")

    surviving_vs_collapsing = tf_genes.merge(aggregated_metrics, on='gene_id')
    surviving_vs_collapsing = surviving_vs_collapsing[['gene_id', 'collapse_rate', 'out_degree', 'in_degree', 'SCC_id', 'group']]

    surviving_vs_collapsing_path = os.path.join(results_surviving_dir, 'surviving_vs_collapsing_TF.tsv')
    surviving_vs_collapsing.to_csv(surviving_vs_collapsing_path, sep='\t', index=False)
    print(f"    Saved: {surviving_vs_collapsing_path}")

    plot_surviving_tf_analysis(surviving_vs_collapsing, plots_surviving_dir)


def generate_visualizations(results_df: pd.DataFrame, steady_seeds: List[int],
                           steady_seeds_info: pd.DataFrame,
                           gene_summary: pd.DataFrame, class_summary: pd.DataFrame):
    print("\n    Generating trajectory plots...")
    for seed in steady_seeds:
        generate_trajectory_plot(seed, results_df)
    print(f"    Saved trajectory plots")

    print("\n    Generating global summary plot...")
    plot_global_summary(steady_seeds)

    print("\n    Generating regime transition summary...")
    plot_regime_transition_summary(results_df, gene_summary)

    print("\n    Generating gene-level effects plot...")
    plot_gene_level_effects(gene_summary)

    print("\n    Generating class comparison plot...")
    plot_class_comparison(class_summary)

    print("\n    Generating delta P over time plot...")
    plot_delta_mean_P_over_time(steady_seeds)

    print("\n    Generating TF network topology plot...")
    plot_tf_network_topology(steady_seeds_info, results_df)


def generate_trajectory_plot(seed: int, results_df: pd.DataFrame):
    seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
    baseline_traj = load_trajectory_tsv(seed_dir, 'baseline')
    P_baseline = baseline_traj['P_traj']
    T, G = P_baseline.shape

    seed_results = results_df[results_df['seed'] == seed]
    baseline_regime = seed_results['baseline_regime'].values[0] if len(seed_results) > 0 else 'unknown'

    tf_gene_set = set(TF_GENES)
    epi_gene_set = set(EPI_GENES)
    time_axis = np.arange(1, T)

    from matplotlib.lines import Line2D

    gene_legend_elements = [Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                            Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                            Line2D([0], [0], color='gray', label='Others')]

    fig = plt.figure(figsize=(18, 16))
    gs = gridspec.GridSpec(4, 3, figure=fig, hspace=0.3, wspace=0.2)

    def get_gene_style(g, highlight_idx):
        if highlight_idx is not None and g == highlight_idx:
            return 'black', 3.0, 1.0
        elif g in tf_gene_set:
            return plt.cm.Set1(0), 1.2, 0.7
        elif g in epi_gene_set:
            return plt.cm.Set1(3), 1.2, 0.7
        else:
            return 'gray', 0.6, 0.4

    def draw_trajectory(ax, P_data, time_axis, highlight_idx=None, show_legend=True):
        y_max = P_data.max()
        for g in range(P_data.shape[1]):
            color, lw, alpha = get_gene_style(g, highlight_idx)
            ax.plot(time_axis, P_data[:, g], '-', color=color, alpha=alpha, linewidth=lw)
        ax.set_ylim(0, y_max * 1.05)
        ax.set_xlim(0, len(time_axis))
        if show_legend:
            ax.legend(handles=gene_legend_elements, loc='upper right', fontsize=8)

    ax = fig.add_subplot(gs[0, 0])
    draw_trajectory(ax, P_baseline[1:], time_axis)
    baseline_title_color = '#E74C3C' if 'collapse' in baseline_regime else 'black'
    ax.set_title(f'Baseline (Seed {seed}): {baseline_regime}', fontsize=12, fontweight='bold', color=baseline_title_color)
    ax.set_xlabel('Time')
    ax.set_ylabel('Protein')
    ax.grid(True, alpha=0.3)

    ax_empty = fig.add_subplot(gs[0, 1])
    ax_empty.axis('off')

    ax_empty2 = fig.add_subplot(gs[0, 2])
    ax_empty2.axis('off')

    for row in range(1, 4):
        for col in range(3):
            idx = (row - 1) * 3 + col
            if idx < len(TARGET_GENES):
                gene_idx = TARGET_GENES[idx]
                gene_type = 'TF' if gene_idx in TF_GENES else 'EPI'

                gene_result = seed_results[seed_results['gene_id'] == gene_idx]
                transition_type = gene_result['transition_type'].values[0] if len(gene_result) > 0 else 'unknown'

                ko_traj = load_trajectory_tsv(seed_dir, f'KO_gene_{gene_idx}')
                P_ko = ko_traj['P_traj']

                ax = fig.add_subplot(gs[row, col])
                draw_trajectory(ax, P_ko[1:], time_axis, highlight_idx=gene_idx, show_legend=True)

                title_color = '#E74C3C' if 'collapse' in transition_type else 'black'
                ax.set_title(f'KO Gene {gene_idx} ({gene_type}): {transition_type}', fontsize=10, fontweight='bold', color=title_color)
                ax.set_xlabel('Time')
                ax.set_ylabel('Protein')
                ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, 'trajectories', f'seed_{seed}_trajectory.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()


def plot_global_summary(steady_seeds: List[int]):
    n_seeds = len(steady_seeds)
    n_rows = 2
    n_cols = (n_seeds + n_rows - 1) // n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 12))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    tf_gene_set = set(TF_GENES)
    epi_gene_set = set(EPI_GENES)

    for idx, seed in enumerate(steady_seeds):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
        baseline_traj = load_trajectory_tsv(seed_dir, 'baseline')
        P_baseline = baseline_traj['P_traj']

        time_axis = np.arange(1, P_baseline.shape[0])
        P_baseline_mean = P_baseline[1:].mean(axis=1)
        ax.plot(time_axis, P_baseline_mean, '-', color='blue', label='Baseline', linewidth=2)

        for gene_idx in TARGET_GENES:
            ko_traj = load_trajectory_tsv(seed_dir, f'KO_gene_{gene_idx}')
            P_ko = ko_traj['P_traj']
            P_ko_mean = P_ko[1:].mean(axis=1)

            if gene_idx in tf_gene_set:
                color = plt.cm.Set1(0)
            elif gene_idx in epi_gene_set:
                color = plt.cm.Set1(3)
            else:
                color = 'gray'

            ax.plot(time_axis, P_ko_mean, '-', color=color, linewidth=1.5, alpha=0.8)

        ax.set_xlabel('Time')
        ax.set_ylabel('Protein')
        ax.set_title(f'Seed {seed}: Mean Protein Trajectory')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, P_baseline.shape[0] - 1)

        from matplotlib.lines import Line2D
        gene_legend = [Line2D([0], [0], color='blue', label='Baseline'),
                       Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                       Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics')]
        ax.legend(handles=gene_legend, loc='upper right', fontsize=10, bbox_to_anchor=(1, 0.9))

    for idx in range(n_seeds, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout()
    global_path = os.path.join(PLOTS_DIR, 'global_summary.png')
    plt.savefig(global_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {global_path}")


def plot_regime_transition_summary(results_df: pd.DataFrame, gene_summary: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    transition_types = results_df['transition_type'].value_counts()
    colors_transition = [COLORS['collapse'] if 'collapse' in t else COLORS['steady'] for t in transition_types.index]
    ax1.pie(transition_types.values, labels=transition_types.index, autopct='%1.1f%%', colors=colors_transition)
    ax1.set_title('Overall Transition Type Distribution')

    ax2 = axes[1]
    gene_summary_sorted = gene_summary.sort_values('collapse_rate', ascending=True)
    colors_genes = [COLORS.get(g, 'gray') for g in gene_summary_sorted['gene_type']]
    ax2.barh(gene_summary_sorted['gene_id'].astype(str), gene_summary_sorted['collapse_rate'], color=colors_genes)
    ax2.set_xlabel('Collapse Rate')
    ax2.set_ylabel('Gene ID')
    ax2.set_title('Gene Collapse Rate Ranking')
    ax2.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)

    plt.tight_layout()
    summary_path = os.path.join(PLOTS_DIR, 'regime_transition', 'transition_summary.png')
    plt.savefig(summary_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {summary_path}")


def plot_gene_level_effects(gene_summary: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    ax1 = axes[0]
    gene_ids = gene_summary['gene_id'].tolist()
    x_pos = np.arange(len(gene_ids))
    bars = ax1.bar(x_pos, gene_summary['collapse_rate'],
                   color=[COLORS['TF'] if g in TF_GENES else COLORS['EPI'] for g in gene_ids])
    ax1.set_xticks(x_pos)
    ax1.set_xticklabels([f'G{g}' for g in gene_ids], fontsize=9)
    ax1.set_xlabel('Gene', fontsize=11)
    ax1.set_ylabel('Collapse Rate', fontsize=11)
    ax1.set_title('Per-Gene Collapse Rate', fontsize=12, fontweight='bold')
    ax1.set_ylim(0, 1)
    ax1.grid(True, alpha=0.3, axis='y')
    for bar, rate in zip(bars, gene_summary['collapse_rate']):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
                f'{rate:.2f}', ha='center', va='bottom', fontsize=8)

    ax2 = axes[1]
    delta_P_terminal = gene_summary['mean_Δmean_P_terminal'].tolist()
    colors = [COLORS['TF'] if g in TF_GENES else COLORS['EPI'] for g in gene_ids]
    bars = ax2.bar(x_pos, delta_P_terminal, color=colors)
    ax2.set_xticks(x_pos)
    ax2.set_xticklabels([f'G{g}' for g in gene_ids], fontsize=9)
    ax2.set_xlabel('Gene', fontsize=11)
    ax2.set_ylabel('Δmean_P_terminal', fontsize=11)
    ax2.set_title('Per-Gene Δmean_P_terminal', fontsize=12, fontweight='bold')
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.5)
    ax2.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    plot_file = os.path.join(PLOTS_DIR, 'gene_level_effects.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {plot_file}")


def plot_class_comparison(class_summary: pd.DataFrame):
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax1 = axes[0]
    class_summary_plot = class_summary.set_index('gene_type')
    class_summary_plot[['average_collapse_rate']].plot(kind='bar', ax=ax1, color=[plt.cm.Set1(0), plt.cm.Set1(3)])
    ax1.set_xlabel('Gene Class')
    ax1.set_ylabel('Collapse Rate')
    ax1.set_title('Collapse Rate: TF vs Epigenetic Genes')
    ax1.legend().remove()
    ax1.set_xticklabels(ax1.get_xticklabels(), rotation=0)
    ax1.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

    ax2 = axes[1]
    class_summary_plot[['average_Δmean_P_terminal']].plot(kind='bar', ax=ax2, color=[plt.cm.Set1(0), plt.cm.Set1(3)])
    ax2.set_xlabel('Gene Class')
    ax2.set_ylabel('Δmean_P_terminal')
    ax2.set_title('Δmean_P_terminal: TF vs Epigenetic Genes')
    ax2.legend().remove()
    ax2.set_xticklabels(ax2.get_xticklabels(), rotation=0)
    ax2.axhline(y=0, color='gray', linestyle='-', alpha=0.5)

    plt.tight_layout()
    comparison_path = os.path.join(PLOTS_DIR, 'class_comparison', 'TF_vs_EPI_comparison.png')
    plt.savefig(comparison_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {comparison_path}")


def plot_delta_mean_P_over_time(steady_seeds: List[int]):
    n_seeds = len(steady_seeds)
    n_rows = 2
    n_cols = (n_seeds + n_rows - 1) // n_rows

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 12))
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    tf_gene_set = set(TF_GENES)
    epi_gene_set = set(EPI_GENES)

    for idx, seed in enumerate(steady_seeds):
        row = idx // n_cols
        col = idx % n_cols
        ax = axes[row, col]

        seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
        baseline_traj = load_trajectory_tsv(seed_dir, 'baseline')
        P_baseline = baseline_traj['P_traj']

        time_axis = np.arange(1, P_baseline.shape[0])

        for gene_idx in TARGET_GENES:
            ko_traj = load_trajectory_tsv(seed_dir, f'KO_gene_{gene_idx}')
            P_ko = ko_traj['P_traj']
            delta_P = P_ko[1:] - P_baseline[1:]
            delta_P_mean = delta_P.mean(axis=1)

            if gene_idx in tf_gene_set:
                color = plt.cm.Set1(0)
            elif gene_idx in epi_gene_set:
                color = plt.cm.Set1(3)
            else:
                color = 'gray'

            ax.plot(time_axis, delta_P_mean, '-', color=color, linewidth=1.5, alpha=0.8)

        ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)

        from matplotlib.lines import Line2D
        gene_legend = [Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                       Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics')]
        ax.legend(handles=gene_legend, loc='upper right', fontsize=12, bbox_to_anchor=(1, 0.9))

        ax.set_xlabel('Time')
        ax.set_ylabel('ΔProtein')
        ax.set_title(f'Seed {seed}: Δmean_P(t) per KO')
        ax.grid(True, alpha=0.3)
        ax.set_xlim(0, P_baseline.shape[0] - 1)

    for idx in range(n_seeds, n_rows * n_cols):
        row = idx // n_cols
        col = idx % n_cols
        axes[row, col].axis('off')

    plt.tight_layout()
    delta_path = os.path.join(PLOTS_DIR, 'delta_P_over_time.png')
    plt.savefig(delta_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {delta_path}")


def plot_surviving_tf_analysis(surviving_vs_collapsing: pd.DataFrame, plots_surviving_dir: str):
    surviving_degrees = surviving_vs_collapsing[surviving_vs_collapsing['group'] == 'surviving']['out_degree']
    collapsing_degrees = surviving_vs_collapsing[surviving_vs_collapsing['group'] == 'collapsing']['out_degree']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].bar(['surviving', 'collapsing'],
                [surviving_degrees.mean() if len(surviving_degrees) > 0 else 0,
                 collapsing_degrees.mean() if len(collapsing_degrees) > 0 else 0],
                yerr=[surviving_degrees.std() if len(surviving_degrees) > 0 else 0,
                      collapsing_degrees.std() if len(collapsing_degrees) > 0 else 0],
                color=['#27AE60', '#E74C3C'])
    axes[0].set_ylabel('Mean Out Degree')
    axes[0].set_title('Degree Comparison: Surviving vs Collapsing TF')

    scc_counts = surviving_vs_collapsing['group'].value_counts()
    axes[1].bar(scc_counts.index, scc_counts.values, color=['#27AE60', '#E74C3C'])
    axes[1].set_ylabel('Count')
    axes[1].set_title('SCC Membership Distribution')

    colors = surviving_vs_collapsing['group'].map({'surviving': '#27AE60', 'collapsing': '#E74C3C'})
    axes[2].scatter(surviving_vs_collapsing['out_degree'], surviving_vs_collapsing['collapse_rate'], c=colors)
    axes[2].set_xlabel('Out Degree')
    axes[2].set_ylabel('Collapse Rate')
    axes[2].set_title('Collapse Rate vs Out Degree')

    plt.tight_layout()
    plot_path = os.path.join(plots_surviving_dir, 'surviving_tf_analysis.png')
    plt.savefig(plot_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {plot_path}")


def plot_tf_network_topology(steady_seeds_info: pd.DataFrame, results_df: pd.DataFrame):
    seeds_collapse = []
    seeds_steady = []
    for _, row in steady_seeds_info.iterrows():
        seed = row['seed']
        seed_results = results_df[results_df['seed'] == seed]
        if len(seed_results) > 0:
            regime = seed_results['baseline_regime'].values[0]
        else:
            seed_dir = os.path.join(DATA_DIR, 'trajectories', f'seed_{seed}')
            baseline_traj = load_trajectory_tsv(seed_dir, 'baseline')
            regime = classify_regime(torch.tensor(baseline_traj['X_traj']))
        if 'collapse' in regime:
            seeds_collapse.append(seed)
        else:
            seeds_steady.append(seed)

    n_collapse = len(seeds_collapse)
    n_steady = len(seeds_steady)
    n_cols = 4

    if n_collapse > 0:
        n_collapse_rows = (n_collapse + n_cols - 1) // n_cols
    else:
        n_collapse_rows = 0
    if n_steady > 0:
        n_steady_rows = (n_steady + n_cols - 1) // n_cols
    else:
        n_steady_rows = 0
    n_rows = n_collapse_rows + n_steady_rows

    if n_rows == 0:
        print("    No networks to plot")
        return

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 5 * n_rows))
    axes = list(np.atleast_1d(axes).flatten())

    def plot_network(ax, seed, regime, world):
        a_ij = world.a_ij

        G = nx.DiGraph()
        G.add_nodes_from(TF_GENES)

        for source in TF_GENES:
            if source not in a_ij:
                continue
            for target in TF_GENES:
                if target in a_ij[source] and a_ij[source][target] > 0:
                    G.add_edge(source, target, weight=a_ij[source][target])

        scc = list(nx.strongly_connected_components(G))
        largest_scc = max(scc, key=len) if scc else set()
        scc_size = len(largest_scc)

        try:
            cycles = list(nx.simple_cycles(G))
            n_cycles = len(cycles)
        except:
            n_cycles = 0

        total_degrees = dict(G.degree())

        try:
            pos = nx.circular_layout(G, scale=0.8)
        except:
            pos = nx.random_layout(G, seed=42)

        node_colors = []
        node_sizes = []
        for node in G.nodes():
            if node in largest_scc:
                node_colors.append('#4C7CB8')
            else:
                node_colors.append('#CCCCCC')
            node_size = 700 + total_degrees.get(node, 0) * 200
            node_sizes.append(node_size)

        edge_colors = []
        for u, v in G.edges():
            weight = G[u][v].get('weight', 0)
            if weight > 0:
                edge_colors.append('#E04532')
            else:
                edge_colors.append('#4C7CB8')

        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                              width=3.0, alpha=0.7, arrows=True,
                              arrowsize=25, connectionstyle="arc3,rad=0.1")
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                              node_size=node_sizes, alpha=0.9)
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=16,
                               font_weight='bold')

        title_color = '#E74C3C' if 'collapse' in regime else '#27AE60'
        ax.set_title(f'{regime.upper()} - Seed {seed}\n'
                    f'SCC: {scc_size} | Cycles: {n_cycles}', fontsize=18, color=title_color)
        ax.axis('off')
        ax.grid(False)

    for idx, seed in enumerate(seeds_collapse):
        world = sample_world(seed)
        plot_network(axes[idx], seed, 'collapse', world)

    steady_start_idx = n_collapse_rows * n_cols
    for idx, seed in enumerate(seeds_steady):
        world = sample_world(seed)
        plot_network(axes[steady_start_idx + idx], seed, 'steady', world)

    for idx in range(n_collapse, steady_start_idx):
        axes[idx].axis('off')

    for idx in range(steady_start_idx + n_steady, len(axes)):
        axes[idx].axis('off')

    plt.suptitle('TF Regulatory Network Topology', fontsize=24, y=1)
    plt.tight_layout()

    plot_file = os.path.join(PLOTS_DIR, 'TF_network_topology.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"    Saved: {plot_file}")


if __name__ == "__main__":
    run_analysis()
