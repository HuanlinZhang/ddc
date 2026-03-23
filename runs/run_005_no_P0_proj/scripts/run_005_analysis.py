"""
Run 005 - Phase 0 Regime Analysis Pipeline
=========================================

Task 1: Resource Projection Ablation
Task 2: Effective Gain Proxy Analysis
Task 3: TF Support / Rescue
Task 4: Initial Condition Scaling

Author: zhanghl
Date: 2026-03-17
"""

import os
import sys
import json
import torch
from torch import Tensor
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')

import ddc
from ddc import compute_TFinput

BASE_DIR = '/home/zhanghl/projects/ddc_github/runs/run_005_no_P0_proj'
DATA_DIR = os.path.join(BASE_DIR, 'data')
PLOTS_DIR = os.path.join(BASE_DIR, 'plots')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')

DIR_ENABLE = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_enable_resource_projection'
DIR_DISABLE = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_disable_resource_projection'

THRESHOLD = 0.1


def setup_directories():
    """
    Create necessary output directories according to Output Spec.
    
    Structure:
    run_005/
    ├── data/
    │   ├── world_parameters/
    │   ├── trajectories/
    │   └── derived_metrics/
    ├── results/
    │   ├── resource_ablation/
    │   ├── gain_proxy_analysis/
    │   ├── TF_rescue/
    │   └── initial_condition_scaling/
    ├── plots/
    ├── scripts/
    └── execution_record.md
    """
    # Main directories
    for d in [DATA_DIR, PLOTS_DIR, RESULTS_DIR]:
        os.makedirs(d, exist_ok=True)
    
    # Data subdirectories (per Output Spec section 2)
    data_subdirs = ['world_parameters', 'trajectories', 'derived_metrics']
    for subdir in data_subdirs:
        os.makedirs(os.path.join(DATA_DIR, subdir), exist_ok=True)
    
    # Results subdirectories (per Output Spec section 2)
    results_subdirs = [
        'resource_ablation',
        'gain_proxy_analysis',
        'TF_rescue',
        'initial_condition_scaling'
    ]
    for subdir in results_subdirs:
        os.makedirs(os.path.join(RESULTS_DIR, subdir), exist_ok=True)

    # Plots subdirectories
    plots_subdirs = [
        'resource_ablation',
        'gain_proxy_analysis',
        'TF_rescue',
        'initial_condition_scaling'
    ]
    for subdir in plots_subdirs:
        os.makedirs(os.path.join(PLOTS_DIR, subdir), exist_ok=True)


def classify_regime(P_terminal: np.ndarray, X_terminal: np.ndarray, TF_genes: List[int], threshold: float = THRESHOLD) -> Tuple[str, int, int]:
    """
    Classify world as collapse or steady based on convergence data.
    Uses X_traj (mRNA) terminal expression to determine regime.
    
    Per run_004 standard:
    - collapse: N_active ≈ 0 (few or no genes with expression > threshold)
    - steady-state: N_active > 0 (multiple genes remain active)
    
    Args:
        P_terminal: Final protein levels, shape (G,) - not used directly
        X_terminal: Final mRNA levels, shape (G,)
        TF_genes: List of TF gene indices
        threshold: Expression threshold for "active" gene
    
    Returns:
        regime: 'collapse' or 'steady'
        n_active: Total number of active genes (all genes, not just TF)
        n_active_TF: Number of active TF genes
    """
    n_active = np.sum(X_terminal > threshold)
    
    tf_mask = np.array([i in TF_genes for i in range(len(X_terminal))])
    tf_expression = X_terminal[tf_mask]
    n_active_TF = np.sum(tf_expression > threshold)
    
    # regime = 'steady' if n_active > 0 else 'collapse'
    regime = 'steady' if n_active > 1 else 'collapse' # seed: 6798 n_active=1 but collapse
    
    return regime, n_active, n_active_TF


def load_seeds_from_directory(data_dir: str) -> List[int]:
    """Extract seed numbers from trajectory files in directory."""
    seeds = []
    for f in os.listdir(data_dir):
        if f.endswith('_traj.pt') and f.startswith('seed_'):
            seed_str = f.replace('seed_', '').replace('_traj.pt', '')
            seeds.append(int(seed_str))
    return sorted(seeds)


def load_results_from_directory(data_dir: str) -> pd.DataFrame:
    """Load simulation results from a directory."""
    seeds = load_seeds_from_directory(data_dir)
    
    results = []
    for seed in seeds:
        traj_path = os.path.join(data_dir, f'seed_{seed}_traj.pt')
        
        saved_data = torch.load(traj_path)
        P_terminal = saved_data['P_traj'][-1, :].numpy()
        X_terminal = saved_data['X_traj'][-1, :].numpy()
        
        regime, n_active, n_active_TF = classify_regime(P_terminal, X_terminal, ddc.TF_GENES)
        mean_X = float(np.mean(X_terminal))
        
        results.append({
            'seed': seed,
            'regime': regime,
            'n_active': n_active,
            'n_active_TF': n_active_TF,
            'mean_X_terminal': mean_X
        })
    
    return pd.DataFrame(results)


# ==========================================
# Task 1: Resource Projection Ablation - Helper Functions
# ==========================================

def generate_trajectory_data():
    """
    Generate trajectory data file for Task 1.
    Output: data/trajectories/resource_ablation_trajectories.tsv
    Format: seed, condition, time, gene_id, X, P, Z
    """
    seeds = load_seeds_from_directory(DIR_ENABLE)
    
    trajectory_rows = []
    
    for seed in seeds:
        # Load trajectories from both conditions
        traj_enable = torch.load(os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt'))
        traj_disable = torch.load(os.path.join(DIR_DISABLE, f'seed_{seed}_traj.pt'))
        
        # Extract data
        X_enable = traj_enable['X_traj'].numpy()  # (T, G)
        P_enable = traj_enable['P_traj'].numpy()  # (T, G)
        Z_enable = traj_enable['Z_traj'].numpy()  # (T, G)
        
        X_disable = traj_disable['X_traj'].numpy()
        P_disable = traj_disable['P_traj'].numpy()
        Z_disable = traj_disable['Z_traj'].numpy()
        
        T, G = X_enable.shape
        
        # Add rows for original (with projection)
        for t in range(T):
            for g in range(G):
                trajectory_rows.append({
                    'seed': seed,
                    'condition': 'original',
                    'time': t,
                    'gene_id': g,
                    'X': X_enable[t, g],
                    'P': P_enable[t, g],
                    'Z': Z_enable[t, g]
                })
        
        # Add rows for ablation (without projection)
        for t in range(T):
            for g in range(G):
                trajectory_rows.append({
                    'seed': seed,
                    'condition': 'projection_removed',
                    'time': t,
                    'gene_id': g,
                    'X': X_disable[t, g],
                    'P': P_disable[t, g],
                    'Z': Z_disable[t, g]
                })
    
    df_traj = pd.DataFrame(trajectory_rows)
    
    # Create trajectories directory if needed
    os.makedirs(os.path.join(DATA_DIR, 'trajectories'), exist_ok=True)
    
    output_path = os.path.join(DATA_DIR, 'trajectories', 'resource_ablation_trajectories.tsv')
    df_traj.to_csv(output_path, sep='\t', index=False)
    print(f"Trajectory data saved to: {output_path}")
    print(f"  Total rows: {len(df_traj)}")
    print(f"  Seeds: {len(seeds)}, Time points: {T}, Genes: {G}, Conditions: 2")


def plot_resource_ablation_trajectories(comparison: pd.DataFrame):
    """
    Generate visualization for Task 1.
    Output: plots/seedXXX_resource_ablation_trajectories.png (one per seed)
    
    Each file contains 2 panels:
    - Panel A: X_i(t) mRNA trajectories (all genes)
    - Panel B: P_TF(t) protein trajectories (TF genes only, each TF with distinct color)
    
    Color scheme:
    - Collapse: red
    - Steady: green
    - Line style: solid=original (enable), dashed=projection removed (disable)
    """
    all_seeds = comparison['seed'].tolist()
    
    # Colors for regimes
    COLORS = {
        'collapse': '#E74C3C',  # red
        'steady': '#27AE60',    # green
    }
    
    os.makedirs(PLOTS_DIR, exist_ok=True)
    
    for seed in all_seeds:
        regime_orig = comparison[comparison['seed'] == seed]['regime_original'].values[0]
        regime_abl = comparison[comparison['seed'] == seed]['regime_ablation'].values[0]
        
        fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
        
        # Load trajectory data
        traj_enable = torch.load(os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt'))
        traj_disable = torch.load(os.path.join(DIR_DISABLE, f'seed_{seed}_traj.pt'))
        
        X_enable = traj_enable['X_traj'].numpy()  # (T, G)
        X_disable = traj_disable['X_traj'].numpy()
        P_enable = traj_enable['P_traj'].numpy()
        P_disable = traj_disable['P_traj'].numpy()
        
        T, G = X_enable.shape
        time_axis = np.arange(T)
        
        TF_GENES = list(range(0, 6))
        EPI_GENES = list(range(17, 20))
        tf_gene_set = set(TF_GENES)
        epi_gene_set = set(EPI_GENES)
        
        added_labels_A = set()
        
        # ================== Panel A: X_i(t) mRNA ==================
        ax_a = axes[0]
        
        # Plot all genes: TF, Epigenetics, Others
        for g in range(G):
            if g in tf_gene_set:
                current_color = plt.cm.Set1(0)
                label = "TF"
                alpha, linewidth = 0.85, 1.8
            elif g in epi_gene_set:
                current_color = plt.cm.Set1(3)
                label = "Epigenetics"
                alpha, linewidth = 0.85, 1.8
            else:
                current_color = 'gray'
                label = "Others"
                alpha, linewidth = 0.5, 0.8
            
            plot_label = label if label not in added_labels_A else None
            if label:
                added_labels_A.add(label)
            
            ax_a.plot(time_axis, X_enable[:, g], '-', color=current_color, alpha=alpha, linewidth=linewidth, label=plot_label)
            ax_a.plot(time_axis, X_disable[:, g], '--', color=current_color, alpha=alpha*0.8, linewidth=linewidth*0.8)
        
        ax_a.set_xlabel('Time', fontsize=10)
        ax_a.set_ylabel('X_i(t)', fontsize=10)
        ax_a.set_title(f'Panel A: X_i(t) mRNA (all genes)\n— Original, -- Ablation', fontsize=11, fontweight='bold')
        ax_a.grid(True, alpha=0.3)
        
        from matplotlib.lines import Line2D
        gene_legend = [Line2D([0], [0], color=plt.cm.Set1(0), label='TF'),
                      Line2D([0], [0], color=plt.cm.Set1(3), label='Epigenetics'),
                      Line2D([0], [0], color='gray', label='Others')]
        ax_a.legend(handles=gene_legend, loc='upper right', fontsize=7, ncol=1)
        
        # ================== Panel B: P_TF(t) protein (TF genes only) ==================
        ax_b = axes[1]
        
        # Different colors for each TF gene
        tf_colors = plt.cm.tab10(np.linspace(0, 1, len(TF_GENES)))
        
        for idx, g in enumerate(TF_GENES):
            if g < G:
                color = tf_colors[idx]
                ax_b.plot(time_axis, P_enable[:, g], '-', color=color, alpha=0.9, linewidth=2,
                         label=f'TF{g}')
                ax_b.plot(time_axis, P_disable[:, g], '--', color=color, alpha=0.7, linewidth=1.5)
        
        ax_b.set_xlabel('Time', fontsize=10)
        ax_b.set_ylabel('P_TF(t)', fontsize=10)
        ax_b.set_title(f'Panel B: P_TF(t) protein (TF genes only)\n— Original, -- Ablation', fontsize=11, fontweight='bold')
        ax_b.grid(True, alpha=0.3)
        ax_b.legend(loc='upper right', fontsize=8, ncol=2, title='TF Genes')
        
        plt.suptitle(f'Seed {seed}: Resource Ablation Analysis', fontsize=14, fontweight='bold', y=1.02)
        plt.tight_layout()
        
        output_path = os.path.join(PLOTS_DIR, 'resource_ablation', f'seed{seed}_resource_ablation_trajectories.png')
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()
    
    print(f"Generated {len(all_seeds)} visualization files in: {PLOTS_DIR}")


def plot_resource_ablation_summary(comparison: pd.DataFrame):
    """
    Generate summary visualization for Task 1.
    Output: plots/resource_ablation_summary.png
    
    Panel: Regime transition heatmap
    
                    ablation
                 collapse   steady
    original
    collapse        n1        n2
    steady          n3        n4
    """
    all_seeds = comparison['seed'].tolist()
    
    n11 = len(comparison[(comparison['regime_original'] == 'collapse') & (comparison['regime_ablation'] == 'collapse')])
    n12 = len(comparison[(comparison['regime_original'] == 'collapse') & (comparison['regime_ablation'] == 'steady')])
    n21 = len(comparison[(comparison['regime_original'] == 'steady') & (comparison['regime_ablation'] == 'collapse')])
    n22 = len(comparison[(comparison['regime_original'] == 'steady') & (comparison['regime_ablation'] == 'steady')])
    
    heatmap_data = np.array([[n11, n12], [n21, n22]])
    
    fig, ax = plt.subplots(figsize=(6, 5))
    
    im = ax.imshow(heatmap_data, cmap='Blues', aspect='auto')
    
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(['collapse', 'steady'])
    ax.set_yticklabels(['collapse', 'steady'])
    
    ax.set_xlabel('Ablation', fontsize=12, fontweight='bold')
    ax.set_ylabel('Original', fontsize=12, fontweight='bold')
    ax.set_title('Regime Transition Heatmap', fontsize=14, fontweight='bold')
    
    for i in range(2):
        for j in range(2):
            text = ax.text(j, i, heatmap_data[i, j], ha='center', va='center', 
                          fontsize=16, fontweight='bold', color='white' if heatmap_data[i, j] > heatmap_data.max()/2 else 'black')
    
    plt.colorbar(im, ax=ax, label='Count')
    plt.tight_layout()
    
    output_path = os.path.join(PLOTS_DIR, 'resource_ablation', 'resource_ablation_summary.png')
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Summary visualization saved to: {output_path}")
    plt.close()


# ==========================================
# Task 1: Resource Projection Ablation
# ==========================================

def task1_compare_regimes():
    """
    Main Task 1: Compare regimes with and without resource projection.
    Read existing data from v1_1_1_enable_resource_projection and v1_1_1_disable_resource_projection.
    
    Outputs (per Output Spec):
    1. results/resource_ablation/regime_comparison.tsv
    2. data/trajectories/resource_ablation_trajectories.tsv
    3. plots/seedXXX_resource_ablation_trajectories.png
    """
    print("=" * 60)
    print("Task 1: Resource Projection Ablation")
    print("=" * 60)
    
    # Load summary data
    print("\n--- Loading data with resource projection ENABLED ---")
    df_enable = load_results_from_directory(DIR_ENABLE)
    print(f"Loaded {len(df_enable)} seeds")
    
    print("\n--- Loading data with resource projection DISABLED ---")
    df_disable = load_results_from_directory(DIR_DISABLE)
    print(f"Loaded {len(df_disable)} seeds")
    
    # Merge and create comparison table
    comparison = pd.merge(
        df_enable, df_disable,
        on='seed',
        suffixes=('_original', '_ablation')
    )
    
    # Step 1: Regime Comparison Table (per Output Spec 3.1)
    regime_table = comparison[[
        'seed',
        'regime_original', 'regime_ablation',
        'n_active_TF_original', 'n_active_TF_ablation',
        'mean_X_terminal_original', 'mean_X_terminal_ablation'
    ]].copy()
    
    # Rename columns to match Output Spec
    regime_table.columns = [
        'seed',
        'original_regime', 'ablation_regime',
        'N_active_TF_original', 'N_active_TF_ablation',
        'mean_X_terminal_original', 'mean_X_terminal_ablation'
    ]
    
    output_path = os.path.join(RESULTS_DIR, 'resource_ablation', 'regime_comparison.tsv')
    regime_table.to_csv(output_path, sep='\t', index=False)
    print(f"\nRegime comparison saved to: {output_path}")
    
    # Step 2: Trajectory Data (per Output Spec 3.2)
    print("\n--- Generating trajectory data ---")
    generate_trajectory_data()
    
    # Step 3: Visualization (per Output Spec 3.3)
    print("\n--- Generating visualizations ---")
    plot_resource_ablation_trajectories(comparison)
    plot_resource_ablation_summary(comparison)
    
    return comparison


# ==========================================
# Task 2: Effective Gain Proxy Analysis
# ==========================================

DTYPE = torch.float32

def compute_TFinput_at_t(P: Tensor, P_graph: dict, a_ij: dict, epsilon: float) -> Tensor:
    """Compute TFinput for a single time point."""
    G = len(P)
    tilde_P = P / (torch.sum(P) + epsilon)
    
    TFinput = torch.zeros(G, dtype=DTYPE)
    for i in range(G):
        d_i = len(P_graph.get(i, []))
        if d_i == 0:
            TFinput[i] = 0.0
        else:
            prod = 1.0
            for j in P_graph.get(i, []):
                prod *= (tilde_P[j] ** a_ij[i][j])
            TFinput[i] = prod ** (1.0 / d_i)
    return TFinput


def compute_mean_TFinput_over_K(P_traj: Tensor, K: Tensor, P_graph: dict, a_ij: dict, epsilon: float) -> Tensor:
    """Compute TFinput/K time series. Returns shape (T,)."""
    T = P_traj.shape[0]
    all_ratios = []
    for t in range(T):
        TFinput = compute_TFinput_at_t(P_traj[t], P_graph, a_ij, epsilon)
        ratio = (TFinput / (K + epsilon)).mean()
        all_ratios.append(ratio)
    return torch.stack(all_ratios)


def compute_TF_retention_ratio(P_traj: Tensor, epsilon: float = 1e-8) -> float:
    """Compute TF_retention_ratio = mean_TF(t=40~60) / max_TF(t=0~50)."""
    T = P_traj.shape[0]
    mean_TF_protein = torch.stack([P_traj[t, ddc.TF_GENES].mean() for t in range(T)])
    mean_terminal = mean_TF_protein[40:61].mean()
    max_TF = P_traj[:51, ddc.TF_GENES].max()
    return float((mean_terminal / (max_TF + epsilon)).item())


def compute_terminal_TFinput(P_traj: Tensor, P_graph: dict, a_ij: dict, epsilon: float) -> float:
    """Compute mean TFinput at final time point."""
    TFinput_final = compute_TFinput_at_t(P_traj[-1], P_graph, a_ij, epsilon)
    return float(TFinput_final.mean().item())


def compute_mean_Hill_activation(P_traj: Tensor, K: Tensor, n: Tensor, P_graph: dict, a_ij: dict, epsilon: float) -> float:
    """Compute mean Hill activation over time."""
    T = P_traj.shape[0]
    all_hill = []
    for t in range(T):
        TFinput = compute_TFinput_at_t(P_traj[t], P_graph, a_ij, epsilon)
        hill_term = (TFinput ** n) / ((K ** n) + (TFinput ** n) + epsilon)
        all_hill.append(hill_term)
    return float(torch.mean(torch.stack(all_hill)).item())


def compute_active_TF_count(P_traj: Tensor, threshold: float = 1e-4) -> Tensor:
    """Compute active TF count at each time point. Returns shape (T,), dtype=int64."""
    TF_expression = P_traj[:, ddc.TF_GENES]
    active_counts = (TF_expression > threshold).sum(dim=1)
    return active_counts


def compute_proxy_metrics(data_dir: str, seed: int) -> Dict[str, Any]:
    """
    Compute proxy metrics for a given seed.
    """
    traj_path = os.path.join(data_dir, f'seed_{seed}_traj.pt')
    saved_data = torch.load(traj_path)
    P_traj = saved_data['P_traj']
    X_traj = saved_data['X_traj']
    
    world_dict = saved_data.get('world', {})
    params = world_dict.get('parameters', {})
    P_graph = {int(k): v for k, v in world_dict.get('P_graph', {}).items()}
    a_ij = {int(k): {int(sub_k): float(sub_v) for sub_k, sub_v in v.items()}
            for k, v in params.get('a_ij', {}).items()}
    
    K = torch.tensor(params.get('K', []), dtype=DTYPE)
    n = torch.tensor(params.get('n', []), dtype=DTYPE)
    epsilon = params.get('epsilon', 1e-8)
    
    TFinput_over_K_timeseries = compute_mean_TFinput_over_K(P_traj, K, P_graph, a_ij, epsilon)
    mean_TFinput_over_K = TFinput_over_K_timeseries.mean().item()
    TF_retention_ratio = compute_TF_retention_ratio(P_traj, epsilon)
    terminal_TFinput = compute_terminal_TFinput(P_traj, P_graph, a_ij, epsilon)
    mean_Hill_activation = compute_mean_Hill_activation(P_traj, K, n, P_graph, a_ij, epsilon)
    active_TF_count = compute_active_TF_count(P_traj)
    
    X_terminal = X_traj[-1, :]
    n_active = (X_terminal > THRESHOLD).sum().item()
    regime = 'steady' if n_active > 1 else 'collapse'
    
    metrics = {
        'seed': seed,
        'regime': regime,
        'mean_TFinput_over_K': mean_TFinput_over_K,
        'TFinput_over_K_timeseries': TFinput_over_K_timeseries,
        'TF_retention_ratio': TF_retention_ratio,
        'terminal_TFinput': terminal_TFinput,
        'mean_Hill_activation': mean_Hill_activation,
        'active_TF_count': active_TF_count
    }
    
    return metrics


def generate_task2_visualizations(df_metrics: pd.DataFrame):
    """Generate visualizations for Task 2 proxy metrics."""

    metric_names = ['mean_TFinput_over_K', 'TF_retention_ratio', 'terminal_TFinput',
                    'mean_Hill_activation']
    metric_labels = ['mean(TFinput/K)', 'TF Retention Ratio', 'Terminal TFinput',
                     'mean Hill Activation']

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    axes_flat = axes.flatten()

    colors = {
        ('original', 'collapse'): '#E74C3C',
        ('original', 'steady'): '#27AE60',
        ('projection_removed', 'collapse'): '#E74C3C',
        ('projection_removed', 'steady'): '#27AE60'
    }

    patterns = {
        'original': '',
        'projection_removed': '///'
    }

    for idx, (metric, label) in enumerate(zip(metric_names, metric_labels)):
        ax = axes_flat[idx]

        data_to_plot = []
        tick_labels = []
        box_colors = []

        for condition in ['original', 'projection_removed']:
            for regime_type in ['collapse', 'steady']:
                subset = df_metrics[(df_metrics['condition'] == condition) & (df_metrics['regime'] == regime_type)]
                vals = subset[metric].values
                data_to_plot.append(vals)
                tick_labels.append(f'{condition[:4]}-{regime_type[:4]}')
                box_colors.append(colors[(condition, regime_type)])

        bp = ax.boxplot(data_to_plot, tick_labels=tick_labels, patch_artist=True)

        for i, box in enumerate(bp['boxes']):
            box.set_facecolor(box_colors[i])
            condition_key = tick_labels[i].split('-')[0]
            if condition_key == 'orig':
                box.set_hatch('')
            else:
                box.set_hatch('///')

        ax.set_ylabel(label, fontsize=10)
        ax.set_title(label, fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.tick_params(axis='x', rotation=45)

    plt.suptitle('Task 2: Proxy Metrics - Collapse vs Steady (Original vs Projection Removed)', fontsize=14, fontweight='bold', y=1)
    plt.tight_layout()
    boxplot_path = os.path.join(PLOTS_DIR, 'gain_proxy_analysis', 'gain_proxy_boxplots.png')
    plt.savefig(boxplot_path, dpi=150, bbox_inches='tight')
    print(f"Boxplot saved to: {boxplot_path}")
    plt.close()

    print("\n--- Generating active TF count time series ---")

    fig, ax = plt.subplots(figsize=(10, 6))

    time_axis = np.arange(df_metrics['active_TF_count'].iloc[0].shape[0])

    has_legend = {('original', 'collapse'): False, ('original', 'steady'): False,
                  ('projection_removed', 'collapse'): False, ('projection_removed', 'steady'): False}
    for idx, row in df_metrics.iterrows():
        active_TF_series = row['active_TF_count'].numpy()
        condition = row['condition']
        regime = row['regime']
        color = '#E74C3C' if regime == 'collapse' else '#27AE60'
        linestyle = '--' if condition == 'projection_removed' else '-'
        key = (condition, regime)
        label = f'{condition}-{regime}' if not has_legend[key] else None
        if not has_legend[key]:
            has_legend[key] = True
        ax.plot(time_axis, active_TF_series, color=color, alpha=0.5, linewidth=1.5,
                linestyle=linestyle, label=label)

    ax.set_xlabel('Time (t)', fontsize=12)
    ax.set_ylabel('Active TF Count', fontsize=12)
    ax.set_title('Active TF Count Over Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='right', fontsize=10)
    
    scatter_path = os.path.join(PLOTS_DIR, 'gain_proxy_analysis', 'active_TF_count_timeseries.png')
    plt.savefig(scatter_path, dpi=150, bbox_inches='tight')
    print(f"Active TF count time series saved to: {scatter_path}")
    plt.close()
    
    print("\n--- Generating TFinput/K time series scatter plot ---")

    fig, ax = plt.subplots(figsize=(10, 6))

    colors = {'collapse': '#E74C3C', 'steady': '#27AE60'}
    has_legend = {('original', 'collapse'): False, ('original', 'steady'): False,
                  ('projection_removed', 'collapse'): False, ('projection_removed', 'steady'): False}

    for idx, row in df_metrics.iterrows():
        condition = row['condition']
        regime = row['regime']
        key = (condition, regime)

        tfinput_series = row['TFinput_over_K_timeseries'].numpy()
        time_axis = np.arange(len(tfinput_series))

        color = colors[regime]
        linestyle = '--' if condition == 'projection_removed' else '-'
        label = f'{condition}-{regime}' if not has_legend[key] else None
        if not has_legend[key]:
            has_legend[key] = True

        ax.plot(time_axis, tfinput_series, color=color, alpha=0.5, linewidth=1.5,
                linestyle=linestyle, label=label)
    
    ax.set_xlabel('Time (t)', fontsize=12)
    ax.set_ylabel('TFinput / K', fontsize=12)
    ax.set_title('TFinput/K Over Time', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(loc='right', fontsize=10)
    
    tfinput_scatter_path = os.path.join(PLOTS_DIR, 'gain_proxy_analysis', 'tfinput_over_k_timeseries.png')
    plt.savefig(tfinput_scatter_path, dpi=150, bbox_inches='tight')
    print(f"TFinput/K time series saved to: {tfinput_scatter_path}")
    plt.close()


def task2_effective_gain_proxy():
    """
    Task 2: Identify dynamical metrics that capture effective regulatory gain.

    Goal: Find dynamical metrics that can distinguish collapse vs steady-state worlds.

    Computes proxy metrics:
    - mean_TFinput_over_K
    - TF_retention_ratio
    - terminal_TFinput
    - mean_Hill_activation
    - active_TF_count

    Data: Original system data (with resource projection) and ablation data (without projection)
    Comparison: collapse vs steady for each condition

    Outputs:
    - results/gain_proxy_analysis/proxy_metrics.tsv
    - results/gain_proxy_analysis/viability_proxy_summary.tsv
    - plots/gain_proxy_analysis/gain_proxy_boxplots.png
    """
    print("=" * 60)
    print("Task 2: Effective Gain Proxy Analysis")
    print("=" * 60)

    os.makedirs(os.path.join(RESULTS_DIR, 'gain_proxy_analysis'), exist_ok=True)

    seeds_enable = load_seeds_from_directory(DIR_ENABLE)
    seeds_disable = load_seeds_from_directory(DIR_DISABLE)

    print(f"\n--- Computing proxy metrics for {len(seeds_enable)} seeds (original) ---")

    all_metrics = []

    for seed in seeds_enable:
        metrics = compute_proxy_metrics(DIR_ENABLE, seed)
        metrics['condition'] = 'original'
        all_metrics.append(metrics)
        print(f"  Seed {seed}: regime={metrics['regime']}, TF_retention={metrics['TF_retention_ratio']:.4f}")

    print(f"\n--- Computing proxy metrics for {len(seeds_disable)} seeds (projection_removed) ---")

    for seed in seeds_disable:
        metrics = compute_proxy_metrics(DIR_DISABLE, seed)
        metrics['condition'] = 'projection_removed'
        all_metrics.append(metrics)
        print(f"  Seed {seed}: regime={metrics['regime']}, TF_retention={metrics['TF_retention_ratio']:.4f}")

    df_metrics = pd.DataFrame(all_metrics)

    df_for_save = df_metrics.copy()
    df_for_save['TFinput_over_K_timeseries'] = df_for_save['TFinput_over_K_timeseries'].apply(
        lambda x: x.numpy().tolist() if hasattr(x, 'numpy') else x)
    df_for_save['active_TF_count'] = df_for_save['active_TF_count'].apply(
        lambda x: x.numpy().tolist() if hasattr(x, 'numpy') else x)
    df_for_save['N_active_TF'] = df_for_save['active_TF_count'].apply(
        lambda x: int(x[-1].item()) if hasattr(x, '__iter__') and hasattr(x[-1], 'item') else int(x[-1]))

    output_path = os.path.join(RESULTS_DIR, 'gain_proxy_analysis', 'proxy_metrics.tsv')
    df_for_save.to_csv(output_path, sep='\t', index=False)
    print(f"Proxy metrics saved to: {output_path}")

    viability_proxy_df = df_for_save[['seed', 'condition', 'mean_TFinput_over_K', 'TF_retention_ratio',
                                      'terminal_TFinput', 'mean_Hill_activation', 'N_active_TF', 'regime']].copy()
    viability_proxy_path = os.path.join(RESULTS_DIR, 'gain_proxy_analysis', 'viability_proxy_summary.tsv')
    viability_proxy_df.to_csv(viability_proxy_path, sep='\t', index=False)
    print(f"Viability proxy summary saved to: {viability_proxy_path}")

    print("\n--- Regime distribution by condition ---")
    print(df_metrics.groupby(['condition', 'regime']).size().unstack(fill_value=0))

    print("\n--- Generating visualizations ---")

    generate_task2_visualizations(df_metrics)

    print("\n--- Summary Statistics ---")
    print(df_metrics[['seed', 'condition', 'regime', 'TF_retention_ratio']].to_string(index=False))

    print("\nTask 2 Complete!")


# ==========================================
# Task 3: TF Support / Rescue Test
# ==========================================

def apply_tf_overexpression(world: Any, X0: Tensor, P0: Tensor, Z0: Tensor, N0: float,
                           zeta: float) -> Dict[str, Tensor]:
    """Apply TF overexpression intervention via mRNA."""
    X0_overexpressed = X0.clone()
    X0_overexpressed[ddc.TF_GENES] *= zeta
    return ddc.simulate_single_cell(world, X0_overexpressed, P0, Z0, N0)


def apply_tf_clamp(world: Any, X0: Tensor, P0: Tensor, Z0: Tensor, N0: float,
                  clamp_window: int, threshold: float) -> Dict[str, Tensor]:
    """Apply TF clamp intervention in early time window."""
    T = ddc.T
    X, P, Z, N = X0.clone(), P0.clone(), Z0.clone(), N0
    P[ddc.TF_GENES] = torch.clamp(P[ddc.TF_GENES], min=threshold)
    ### 0323 Do not perform projection on P0
    # P = ddc.apply_resource_projection(P, world)
    ###
    X_traj = torch.zeros((T+1, ddc.G), dtype=ddc.DTYPE)
    P_traj = torch.zeros((T+1, ddc.G), dtype=ddc.DTYPE)
    Z_traj = torch.zeros((T+1, ddc.G), dtype=ddc.DTYPE)
    N_traj = torch.zeros((T+1,), dtype=ddc.DTYPE)
    
    X_traj[0], P_traj[0], Z_traj[0], N_traj[0] = X, P, Z, N
    
    for t in range(T):
        tilde_P = ddc.normalize_protein(P, world)
        TFinput = ddc.compute_TFinput(tilde_P, world)
        Z_next = ddc.update_chromatin(tilde_P, world)
        X_next = ddc.update_mRNA(X, Z, TFinput, world)
        P_raw = ddc.update_protein_raw(P, X, world)
        P_next = ddc.apply_resource_projection(P_raw, world)

        if t < clamp_window:
            P_next[ddc.TF_GENES] = torch.clamp(P_next[ddc.TF_GENES], min=threshold)
            P_next = ddc.apply_resource_projection(P_next, world)

        N_next = ddc.update_fate(N, world)
        
        X, P, Z, N = X_next, P_next, Z_next, N_next
        X_traj[t+1], P_traj[t+1], Z_traj[t+1], N_traj[t+1] = X, P, Z, N
    
    return {'X_traj': X_traj, 'P_traj': P_traj, 'Z_traj': Z_traj, 'N_traj': N_traj}


def task3_tf_rescue():
    """
    Task 3: Test whether increased TF support can rescue collapsed worlds.
    
    Methods:
    - Method A: TF Overexpression (P_TF_initial = ζ × original)
    - Method B: TF Clamp (enforce P_TF >= threshold in early window)
    """
    
    print("=" * 60)
    print("Task 3: TF Support / Rescue Test")
    print("=" * 60)
    
    proxy_metrics_path = os.path.join(RESULTS_DIR, 'gain_proxy_analysis', 'proxy_metrics.tsv')
    df_proxy = pd.read_csv(proxy_metrics_path, sep='\t')
    collapse_df = df_proxy[df_proxy['regime'] == 'collapse'].sort_values('TF_retention_ratio')
    all_collapse_seeds = collapse_df['seed'].tolist()
    collapse_seeds = all_collapse_seeds
    
    print(f"All {len(collapse_seeds)} collapse seeds:")
    for s in collapse_seeds:
        ratio = collapse_df[collapse_df['seed'] == s]['TF_retention_ratio'].values[0]
        print(f"  seed {s}: TF_retention_ratio = {ratio:.6f}")
    
    os.makedirs(os.path.join(RESULTS_DIR, 'TF_rescue'), exist_ok=True)
    
    zeta_values = [2, 3, 4, 5]
    print("\n=== Method A: TF Overexpression ===")
    for seed in collapse_seeds:
        print(f"\n--- Processing seed {seed} ---")
        
        traj_path = os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt')
        saved_data = torch.load(traj_path)
        
        world_dict = saved_data['world']
        X0 = saved_data['X_traj'][0]
        P0 = saved_data['P_traj'][0]
        Z0 = saved_data['Z_traj'][0]
        N0 = saved_data['N_traj'][0].item()
        
        world = ddc.World(seed)
        world.from_dict(world_dict)
         
        for zeta in zeta_values:
            print(f"  ζ={zeta}: ", end="")
            
            result = apply_tf_overexpression(world, X0, P0, Z0, N0, zeta)
            
            X_terminal = result['X_traj'][-1]
            n_active = (X_terminal > THRESHOLD).sum().item()
            regime = 'steady' if n_active > 1 else 'collapse'
            
            print(f"regime={regime}")
            
            save_data = {
                'seed': seed,
                'method': 'overexpression',
                'zeta': zeta,
                'regime': regime,
                'X_traj': result['X_traj'],
                'P_traj': result['P_traj'],
                'Z_traj': result['Z_traj'],
                'N_traj': result['N_traj'],
            }
            save_path = os.path.join(RESULTS_DIR, 'TF_rescue', f'seed{seed}_overexp_zeta{zeta}.pt')
            torch.save(save_data, save_path)
    
    clamp_window = 50
    clamp_thresholds = [0.02, 0.04, 0.06, 0.08]
    print("\n=== Method B: TF Clamp ===")
    for seed in collapse_seeds:
        print(f"\n--- Processing seed {seed} ---")
        
        traj_path = os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt')
        saved_data = torch.load(traj_path)
        
        world_dict = saved_data['world']
        X0 = saved_data['X_traj'][0]
        P0 = saved_data['P_traj'][0]
        Z0 = saved_data['Z_traj'][0]
        N0 = saved_data['N_traj'][0].item()
        
        world = ddc.World(seed)
        world.from_dict(world_dict)
        
        for threshold in clamp_thresholds:
            print(f"  window={clamp_window}, threshold={threshold}: ", end="")
            
            result_clamp = apply_tf_clamp(world, X0, P0, Z0, N0, clamp_window, threshold)
            
            X_terminal = result_clamp['X_traj'][-1]
            n_active = (X_terminal > THRESHOLD).sum().item()
            regime = 'steady' if n_active > 1 else 'collapse'
            
            print(f"regime={regime}")
            
            save_data = {
                'seed': seed,
                'method': 'clamp',
                'clamp_window': clamp_window,
                'clamp_threshold': threshold,
                'regime': regime,
                'X_traj': result_clamp['X_traj'],
                'P_traj': result_clamp['P_traj'],
                'Z_traj': result_clamp['Z_traj'],
                'N_traj': result_clamp['N_traj'],
            }
            save_path = os.path.join(RESULTS_DIR, 'TF_rescue', 
                                    f'seed{seed}_clamp_win{clamp_window}_thresh{threshold}.pt')
            torch.save(save_data, save_path)
    
    generate_task3_output_files(collapse_seeds, zeta_values, clamp_window, clamp_thresholds)
    
    print("\nTask 3 Complete!")


def generate_task3_output_files(collapse_seeds: list, zeta_values: list, clamp_window: int, clamp_thresholds: list):
    """Generate output files for Task 3: rescue_summary.tsv and TF_rescue_trajectories.tsv"""
    print("\n--- Generating rescue_summary.tsv ---")
    rows = []
    for seed in collapse_seeds:
        original_regime = 'collapse'
        for zeta in zeta_values:
            path = os.path.join(RESULTS_DIR, 'TF_rescue', f'seed{seed}_overexp_zeta{zeta}.pt')
            data = torch.load(path)
            rows.append({
                'seed': seed,
                'world_type': 'collapse',
                'support_type': 'overexpression',
                'support_target': 'X of TF',
                'support_strength': zeta,
                'original_regime': original_regime,
                'rescued_regime': data['regime']
            })
        for thresh in clamp_thresholds:
            path = os.path.join(RESULTS_DIR, 'TF_rescue', f'seed{seed}_clamp_win{clamp_window}_thresh{thresh}.pt')
            data = torch.load(path)
            rows.append({
                'seed': seed,
                'world_type': 'collapse',
                'support_type': 'clamp',
                'support_target': 'P_next of TF',
                'support_strength': f'window:{clamp_window}, threshold:{thresh}',
                'original_regime': original_regime,
                'rescued_regime': data['regime']
            })
    df_rescue = pd.DataFrame(rows)
    rescue_summary_path = os.path.join(RESULTS_DIR, 'TF_rescue', 'rescue_summary.tsv')
    df_rescue.to_csv(rescue_summary_path, sep='\t', index=False)
    print(f"  Saved: {rescue_summary_path}")
    
    print("\n--- Generating TF_rescue_trajectories.tsv ---")
    traj_rows = []
    for seed in collapse_seeds:
        for zeta in zeta_values:
            path = os.path.join(RESULTS_DIR, 'TF_rescue', f'seed{seed}_overexp_zeta{zeta}.pt')
            data = torch.load(path)
            for t in range(data['P_traj'].shape[0]):
                for gene in range(data['P_traj'].shape[1]):
                    traj_rows.append({
                        'seed': seed,
                        'method': 'overexpression',
                        'support_strength': zeta,
                        'time': t,
                        'gene_id': gene,
                        'X': data['X_traj'][t, gene].item(),
                        'P': data['P_traj'][t, gene].item(),
                        'Z': data['Z_traj'][t, gene].item()
                    })
        for thresh in clamp_thresholds:
            path = os.path.join(RESULTS_DIR, 'TF_rescue', f'seed{seed}_clamp_win{clamp_window}_thresh{thresh}.pt')
            data = torch.load(path)
            for t in range(data['P_traj'].shape[0]):
                for gene in range(data['P_traj'].shape[1]):
                    traj_rows.append({
                        'seed': seed,
                        'method': 'clamp',
                        'support_strength': thresh,
                        'time': t,
                        'gene_id': gene,
                        'X': data['X_traj'][t, gene].item(),
                        'P': data['P_traj'][t, gene].item(),
                        'Z': data['Z_traj'][t, gene].item()
                    })
    df_traj = pd.DataFrame(traj_rows)
    traj_path = os.path.join(DATA_DIR, 'trajectories', 'TF_rescue_trajectories.tsv')
    df_traj.to_csv(traj_path, sep='\t', index=False)
    print(f"  Saved: {traj_path}")
    
    generate_task3_visualizations(collapse_seeds, zeta_values, clamp_window, clamp_thresholds)


def generate_task3_visualizations(collapse_seeds: list, zeta_values: list, clamp_window: int, clamp_thresholds: list):
    """Generate visualizations for Task 3 TF Rescue results."""
    TF_RESCUE_DIR = os.path.join(RESULTS_DIR, 'TF_rescue')
    PLOT_SUBDIR = os.path.join(PLOTS_DIR, 'TF_rescue')
    os.makedirs(PLOT_SUBDIR, exist_ok=True)
    
    TF_COLOR = plt.cm.Set1(0)
    EPI_COLOR = plt.cm.Set1(3)
    OTHER_COLOR = 'gray'
    
    def plot_single_panel(ax, P_traj):
        """Plot P_traj on a single axis with gene categories."""
        T_steps = P_traj.shape[0]
        t = np.arange(T_steps)
        
        added_labels = set()
        for gene in range(P_traj.shape[1]):
            if gene in ddc.TF_GENES:
                color = TF_COLOR
                label = "TF" if "TF" not in added_labels else None
                added_labels.add("TF")
                alpha, lw, zorder = 0.85, 1.8, 2
            elif gene in range(17, 20):
                color = EPI_COLOR
                label = "Epigenetics" if "Epigenetics" not in added_labels else None
                added_labels.add("Epigenetics")
                alpha, lw, zorder = 0.85, 1.8, 2
            else:
                color = OTHER_COLOR
                label = "Others" if "Others" not in added_labels else None
                added_labels.add("Others")
                alpha, lw, zorder = 0.5, 0.8, 1
            
            ax.plot(t, P_traj[:, gene], color=color, alpha=alpha, lw=lw, zorder=zorder, label=label)
        
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        
        handles, labels = ax.get_legend_handles_labels()
        order_map = {"TF": 0, "Epigenetics": 1, "Others": 2}
        pairs = sorted(zip(handles, labels), key=lambda x: order_map.get(x[1], 99))
        if pairs:
            h_sorted, l_sorted = zip(*pairs)
            ax.legend(h_sorted, l_sorted, fontsize=14, loc='upper right', frameon=True)
    
    print("\n--- Generating per-seed visualizations (3x4 grid) ---")
    for seed in collapse_seeds:
        fig, axes = plt.subplots(3, 4, figsize=(24, 16))
        
        traj_path = os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt')
        saved_data = torch.load(traj_path)
        P_traj_original = saved_data['P_traj'].numpy()
        
        ax_orig = axes[0, 0]
        plot_single_panel(ax_orig, P_traj_original)
        ax_orig.set_title(f'Original (collapse)', fontsize=14, fontweight='bold')
        ax_orig.set_ylabel('Protein Expression', fontsize=14)
        ax_orig.set_xlabel('Time Step', fontsize=14)    
        
        for col_idx, zeta in enumerate(zeta_values):
            ax = axes[1, col_idx]
            path = os.path.join(TF_RESCUE_DIR, f'seed{seed}_overexp_zeta{zeta}.pt')
            data = torch.load(path)
            P_traj = data['P_traj'].numpy()
            plot_single_panel(ax, P_traj)
            ax.set_title(f'ζ={zeta} (collapse)', fontsize=14, fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel('Protein Expression', fontsize=14)
            ax.set_xlabel('Time Step', fontsize=14)
        
        for col_idx, thresh in enumerate(clamp_thresholds):
            ax = axes[2, col_idx]
            path = os.path.join(TF_RESCUE_DIR, f'seed{seed}_clamp_win{clamp_window}_thresh{thresh}.pt')
            data = torch.load(path)
            P_traj = data['P_traj'].numpy()
            plot_single_panel(ax, P_traj)
            ax.set_title(f'window={clamp_window} clamp={thresh} (collapse)', fontsize=14, fontweight='bold')
            if col_idx == 0:
                ax.set_ylabel('Protein Expression', fontsize=14)
            ax.set_xlabel('Time Step', fontsize=14)
        
        axes[0, 1].axis('off')
        axes[0, 2].axis('off')
        axes[0, 3].axis('off')
        
        plt.suptitle(f'Task 3: Seed {seed} — Original & TF Rescue Interventions', 
                     fontsize=18, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.98])
        
        save_path = os.path.join(PLOT_SUBDIR, f'seed{seed}_TF_rescue_trajectories.png')
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {save_path}")
    
    print("\nTask 3 visualizations complete!")


# ==========================================
# Task 4: Initial Condition Scaling
# ==========================================

def apply_initial_condition_scaling(world: Any, X0: Tensor, P0: Tensor, Z0: Tensor, N0: float,
                                    gamma: float) -> Dict[str, Tensor]:
    """Apply initial condition scaling: X0 → γ×X0, P0 → γ×P0."""
    X0_scaled = gamma * X0
    P0_scaled = gamma * P0
    ### 0323 Do not perform projection on P0_scaled
    # P0_scaled = ddc.apply_resource_projection(P0_scaled, world)
    ###
    return ddc.simulate_single_cell(world, X0_scaled, P0_scaled, Z0, N0)


def task4_initial_condition_scaling():
    """
    Task 4: Test whether collapse depends on initial condition scaling.
    
    If scaling initial conditions (X0, P0) changes the regime, 
    the system may have multiple attractors.
    """
    
    print("=" * 60)
    print("Task 4: Initial Condition Scaling Test")
    print("=" * 60)
    
    proxy_metrics_path = os.path.join(RESULTS_DIR, 'gain_proxy_analysis', 'proxy_metrics.tsv')
    df_proxy = pd.read_csv(proxy_metrics_path, sep='\t')
    collapse_df = df_proxy[df_proxy['regime'] == 'collapse'].sort_values('TF_retention_ratio')
    collapse_seeds = collapse_df['seed'].tolist()
    
    os.makedirs(os.path.join(RESULTS_DIR, 'initial_condition_scaling'), exist_ok=True)
    
    gamma_values = [0.5, 1.0, 2.0, 5.0]
    
    print("\n=== Initial Condition Scaling ===")
    
    summary_rows = []
    traj_rows = []
    
    for seed in collapse_seeds:
        print(f"\n--- Processing seed {seed} ---")
        
        traj_path = os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt')
        saved_data = torch.load(traj_path)
        
        world_dict = saved_data['world']
        X0_original = saved_data['X_traj'][0]
        P0_original = saved_data['P_traj'][0]
        Z0 = saved_data['Z_traj'][0]
        N0 = saved_data['N_traj'][0].item()
        original_regime = saved_data.get('regime', 'collapse')
        
        world = ddc.World(seed)
        world.from_dict(world_dict)
        
        for gamma in gamma_values:
            print(f"  γ={gamma}: ", end="")
            
            result = apply_initial_condition_scaling(world, X0_original, P0_original, Z0, N0, gamma)
            
            X_terminal = result['X_traj'][-1]
            P_terminal = result['P_traj'][-1]
            mean_X_terminal = X_terminal.mean().item()
            mean_P_terminal = P_terminal.mean().item()
            
            n_active = (X_terminal > THRESHOLD).sum().item()
            scaled_regime = 'steady' if n_active > 1 else 'collapse'
            
            print(f"regime={scaled_regime}")
            
            summary_rows.append({
                'seed': seed,
                'gamma': gamma,
                'original_regime': original_regime,
                'scaled_regime': scaled_regime,
                'mean_X_terminal': mean_X_terminal,
                'mean_P_terminal': mean_P_terminal
            })
            
            for t in range(result['X_traj'].shape[0]):
                for gene in range(result['X_traj'].shape[1]):
                    traj_rows.append({
                        'seed': seed,
                        'gamma': gamma,
                        'time': t,
                        'gene_id': gene,
                        'X': result['X_traj'][t, gene].item(),
                        'P': result['P_traj'][t, gene].item(),
                        'Z': result['Z_traj'][t, gene].item()
                    })
            
            save_data = {
                'seed': seed,
                'gamma': gamma,
                'regime': scaled_regime,
                'X_traj': result['X_traj'],
                'P_traj': result['P_traj'],
                'Z_traj': result['Z_traj'],
                'N_traj': result['N_traj']
            }
            save_path = os.path.join(RESULTS_DIR, 'initial_condition_scaling', f'seed{seed}_gamma{gamma}.pt')
            torch.save(save_data, save_path)
    
    print("\n--- Saving scaling_summary.tsv ---")
    df_summary = pd.DataFrame(summary_rows)
    summary_path = os.path.join(RESULTS_DIR, 'initial_condition_scaling', 'scaling_summary.tsv')
    df_summary.to_csv(summary_path, sep='\t', index=False)
    print(f"  Saved: {summary_path}")
    
    print("\n--- Saving scaling_trajectories.tsv ---")
    df_traj = pd.DataFrame(traj_rows)
    traj_path = os.path.join(DATA_DIR, 'trajectories', 'scaling_trajectories.tsv')
    os.makedirs(os.path.dirname(traj_path), exist_ok=True)
    df_traj.to_csv(traj_path, sep='\t', index=False)
    print(f"  Saved: {traj_path}")
    
    generate_task4_visualizations(collapse_seeds, gamma_values)
    
    print("\nTask 4 Complete!")


def generate_task4_visualizations(collapse_seeds: list, gamma_values: list):
    """Generate visualizations for Task 4 scaling results."""
    SCALING_DIR = os.path.join(RESULTS_DIR, 'initial_condition_scaling')
    PLOT_SUBDIR = os.path.join(PLOTS_DIR, 'initial_condition_scaling')
    os.makedirs(PLOT_SUBDIR, exist_ok=True)
    
    TF_COLOR = plt.cm.Set1(0)
    EPI_COLOR = plt.cm.Set1(3)
    OTHER_COLOR = 'gray'
    
    def plot_single_panel(ax, P_traj, title=None):
        """Plot P_traj on a single axis with gene categories."""
        T_steps = P_traj.shape[0]
        t = np.arange(T_steps)
        
        added_labels = set()
        for gene in range(P_traj.shape[1]):
            if gene in ddc.TF_GENES:
                color = TF_COLOR
                label = "TF" if "TF" not in added_labels else None
                added_labels.add("TF")
                alpha, lw, zorder = 0.85, 1.8, 2
            elif gene in range(17, 20):
                color = EPI_COLOR
                label = "Epigenetics" if "Epigenetics" not in added_labels else None
                added_labels.add("Epigenetics")
                alpha, lw, zorder = 0.85, 1.8, 2
            else:
                color = OTHER_COLOR
                label = "Others" if "Others" not in added_labels else None
                added_labels.add("Others")
                alpha, lw, zorder = 0.5, 0.8, 1
            
            ax.plot(t, P_traj[:, gene], color=color, alpha=alpha, lw=lw, zorder=zorder, label=label)
        
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        
        if title:
            ax.set_title(title, fontsize=14, fontweight='bold')
        
        handles, labels = ax.get_legend_handles_labels()
        order_map = {"TF": 0, "Epigenetics": 1, "Others": 2}
        pairs = sorted(zip(handles, labels), key=lambda x: order_map.get(x[1], 99))
        if pairs:
            h_sorted, l_sorted = zip(*pairs)
            ax.legend(h_sorted, l_sorted, fontsize=14, loc='upper right', frameon=True)
    
    print("\n--- Generating per-seed visualizations (1x5 grid) ---")
    
    for seed in collapse_seeds:
        fig, axes = plt.subplots(1, 5, figsize=(25, 5))
        
        traj_path = os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt')
        saved_data = torch.load(traj_path)
        P_traj_original = saved_data['P_traj'].numpy()
        regime_original = 'unknown'
        
        plot_single_panel(axes[0], P_traj_original, f'Original (regime={regime_original})')
        axes[0].set_ylabel('Protein Expression', fontsize=14)
        
        for col_idx, gamma in enumerate(gamma_values):
            path = os.path.join(SCALING_DIR, f'seed{seed}_gamma{gamma}.pt')
            data = torch.load(path)
            P_traj = data['P_traj'].numpy()
            regime = data.get('regime', 'unknown')
            plot_single_panel(axes[col_idx + 1], P_traj, f'γ={gamma} (regime={regime})')
        
        plt.suptitle(f'Task 4: Seed {seed} — Initial Condition Scaling', fontsize=18, fontweight='bold')
        plt.tight_layout(rect=[0, 0, 1, 0.98])
        
        save_path = os.path.join(PLOT_SUBDIR, f'seed{seed}_scaling_trajectories.png')
        fig.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  Saved: {save_path}")
    
    print("\nTask 4 visualizations complete!")


def main():
    """Main entry point for run_005 analysis."""
    setup_directories()
    
    task1_compare_regimes()
    task2_effective_gain_proxy()
    task3_tf_rescue()
    task4_initial_condition_scaling()
    
    print("\n" + "=" * 60)
    print("Run 005 Analysis Complete")
    print("=" * 60)


if __name__ == '__main__':
    main()
