"""
Run 006 - External Perturbation Response Analysis
=================================================

Objective: Verify (state, perturbation) -> response mapping

Perturbation Design:
    - steady systems: TF knockdown (KD) at t=50
    - collapse systems: TF overexpression (OE) at t=5

Response Metrics:
    - delta_TF: immediate response
    - recovery: most critical metric
    - amplification: auxiliary metric

Author: zhanghl
Date: 2026-03-31
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
from typing import Dict, List, Tuple, Any

# Ensure ddc.py can be found in the system path
sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')

import ddc
from ddc import TF_GENES, EPI_GENES

# ==========================================
# Directory and Global Configuration
# ==========================================
BASE_DIR = '/home/zhanghl/projects/ddc_github/runs/run_006'
DATA_DIR = os.path.join(BASE_DIR, 'data')
RESULTS_DIR = os.path.join(BASE_DIR, 'results')
PLOTS_DIR = os.path.join(BASE_DIR, 'plots')

DIR_ENABLE = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_enable_resource_projection'

T_PERTURB_STEADY = 50
T_PERTURB_COLLAPSE = 5
T_PRE = 50
T_POST = 50
EPSILON_RECOVERY = 0.1
EPSILON_AMP = 1e-4

# Output file paths
BASELINE_TSV_PATH = os.path.join(DATA_DIR, 'trajectories', 'baseline.tsv')
PERTURBED_TSV_PATH = os.path.join(DATA_DIR, 'trajectories', 'perturbed.tsv')


def setup_directories():
    for d in [DATA_DIR, RESULTS_DIR, PLOTS_DIR]:
        os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join(DATA_DIR, 'trajectories'), exist_ok=True)
    os.makedirs(os.path.join(DATA_DIR, 'response_metrics'), exist_ok=True)
    os.makedirs(os.path.join(RESULTS_DIR, 'aggregated'), exist_ok=True)
    os.makedirs(os.path.join(RESULTS_DIR, 'per_seed'), exist_ok=True)

    if os.path.exists(BASELINE_TSV_PATH): os.remove(BASELINE_TSV_PATH)
    if os.path.exists(PERTURBED_TSV_PATH): os.remove(PERTURBED_TSV_PATH)


def load_selected_seeds(proxy_metrics_path: str, n_per_regime: int = 4) -> Tuple[List[int], List[int]]:
    df = pd.read_csv(proxy_metrics_path, sep='\t')
    collapse_seeds = sorted(df[df['regime'] == 'collapse']['seed'].tolist()[:n_per_regime])
    steady_seeds = sorted(df[df['regime'] == 'steady']['seed'].tolist()[:n_per_regime])
    return collapse_seeds, steady_seeds


def select_target_gene(P_traj: torch.Tensor, T_perturb: int) -> int:
    P_TF_at_T = P_traj[T_perturb, TF_GENES]
    target_idx = int(torch.argmax(P_TF_at_T).item())
    return TF_GENES[target_idx]


def compute_response_metrics(
    P_baseline: torch.Tensor,
    P_perturbed: torch.Tensor,
    T_perturb: int
) -> Dict[str, float]:
    """Compute response metrics strictly according to Output Spec (start aligned to T_perturb + 1)"""
    mean_TF = lambda P: P[:, TF_GENES].mean(dim=1)

    delta_TF = (mean_TF(P_perturbed)[T_perturb + 5] - mean_TF(P_perturbed)[T_perturb + 1]).item()

    P_baseline_end = P_baseline[-1].mean().item()
    P_perturbed_end = P_perturbed[-1].mean().item()
    recovery = 1 if abs(P_perturbed_end - P_baseline_end) < EPSILON_RECOVERY else 0

    deviation = torch.abs(P_perturbed - P_baseline).sum(dim=1)
    initial_deviation = deviation[T_perturb + 1].item()
    max_deviation = deviation[T_perturb + 1:].max().item()

    if initial_deviation < EPSILON_AMP:
        amplification = float('nan')
    else:
        amplification = max_deviation / initial_deviation

    return {
        'delta_TF': delta_TF,
        'recovery': recovery,
        'amplification': amplification,
        'initial_deviation': initial_deviation,
        'max_deviation': max_deviation
    }


def save_trajectory_tsv(traj_dict: Dict[str, torch.Tensor], seed: int, save_path: str):
    """Convert 3D trajectory into a standard long-format TSV and append to file"""
    T_steps, G = traj_dict['X_traj'].shape
    
    X_np = traj_dict['X_traj'].numpy()
    P_np = traj_dict['P_traj'].numpy()
    Z_np = traj_dict['Z_traj'].numpy()

    records = []
    for t in range(T_steps):
        for g in range(G):
            records.append({
                'seed': seed,
                'time': t,
                'gene_id': g,
                'X': X_np[t, g],
                'P': P_np[t, g],
                'Z': Z_np[t, g]
            })
            
    df = pd.DataFrame(records)
    file_exists = os.path.isfile(save_path)
    df.to_csv(save_path, sep='\t', mode='a', header=not file_exists, index=False)


def run_single_sample(seed: int, regime: str, traj_path: str = None) -> Dict[str, Any]:
    if traj_path is None:
        traj_path = os.path.join(DIR_ENABLE, f'seed_{seed}_traj.pt')

    saved_data = torch.load(traj_path)
    world_dict = saved_data['world']

    T_perturb = T_PERTURB_STEADY if regime == 'steady' else T_PERTURB_COLLAPSE

    # Physically isolated initial states
    X0 = saved_data['X_traj'][0].clone()
    P0 = saved_data['P_traj'][0].clone()
    Z0 = saved_data['Z_traj'][0].clone()
    N0 = saved_data['N_traj'][0].item()

    world = ddc.World(seed)
    world.from_dict(world_dict)

    # 1. Baseline simulation
    result_baseline = ddc.simulate_single_cell(
        world, X0.clone(), P0.clone(), Z0.clone(), N0,
        t_steps=T_PRE + T_POST
    )
    P_baseline = result_baseline['P_traj']
    
    # 2. Target gene identification
    target_gene = select_target_gene(P_baseline, T_perturb)

    # 3. Perturbation simulation
    scale_factor = 0.5 if regime == 'steady' else 2.0
    intervention_config = {'scale_P': [(target_gene, scale_factor)]}

    result_perturbed = ddc.simulate_single_cell(
        world, X0.clone(), P0.clone(), Z0.clone(), N0,
        t_steps=T_PRE + T_POST,
        intervention_time=T_perturb,
        intervention_config=intervention_config
    )
    P_perturbed = result_perturbed['P_traj']

    # Note: Since ddc.py v1.1, the post-intervention state is automatically recorded
    # in the trajectory at time T_perturb. No manual modification needed here.

    # 4. Compute metrics (Safe, as all metrics calculation starts from T_perturb + 1)
    response = compute_response_metrics(P_baseline, P_perturbed, T_perturb)

    # 5. Save trajectories to TSV
    save_trajectory_tsv(result_baseline, seed, BASELINE_TSV_PATH)
    save_trajectory_tsv(result_perturbed, seed, PERTURBED_TSV_PATH)

    # Return record for aggregation
    sample_record = {
        'seed': seed,
        'regime': regime,
        'target_gene': target_gene,
        'perturbation_type': 'KD' if regime == 'steady' else 'OE',
        'perturbation_strength': scale_factor,
        'T_perturb': T_perturb,
        'delta_TF': response['delta_TF'],
        'recovery': response['recovery'],
        'amplification': response['amplification'],
        'initial_deviation': response['initial_deviation'],
        'max_deviation': response['max_deviation']
    }

    return sample_record


# ==========================================
# Visualization Functions
# ==========================================

TF_COLOR = plt.cm.Set1(0)
EPI_COLOR = plt.cm.Set1(3)
OTHER_COLOR = 'gray'


def plot_single_panel(ax, P_traj, title=None, fontsize=12, legend_fontsize=10, t_start: int = 1,
                      t_end: int = None, T_perturb: int = None):
    """Plot P_traj on a single axis with gene categories.
    
    Args:
        ax: Matplotlib axis to plot on
        P_traj: Protein trajectory array of shape (T+1, G)
        title: Optional title for the panel
        fontsize: Base font size
        legend_fontsize: Legend font size
        t_start: Starting time step index (default 1, meaning t=1,2,3,...)
        t_end: Ending time step index (inclusive). If None, plot to the end.
        T_perturb: Perturbation time point to mark with vertical line (optional)
    """
    # Slice P_traj based on t_start and t_end
    if t_end is None:
        t_end = P_traj.shape[0] - 1
    P_traj_plot = P_traj[t_start:t_end + 1]
    T_steps = P_traj_plot.shape[0]
    t = np.arange(t_start, t_start + T_steps)
    
    added_labels = set()
    for gene in range(P_traj_plot.shape[1]):
        if gene in TF_GENES:
            color = TF_COLOR
            label = "TF" if "TF" not in added_labels else None
            alpha, lw, zorder = 0.85, 1.8, 2
        elif gene in EPI_GENES:
            color = EPI_COLOR
            label = "Epigenetics" if "Epigenetics" not in added_labels else None
            alpha, lw, zorder = 0.85, 1.8, 2
        else:
            color = OTHER_COLOR
            label = "Others" if "Others" not in added_labels else None
            alpha, lw, zorder = 0.5, 0.8, 1
        
        if label:
            added_labels.add(label.split()[0] if 'Epigenetics' in label else label)
        ax.plot(t, P_traj_plot[:, gene], color=color, alpha=alpha, lw=lw, zorder=zorder, label=label)
    
    # Add perturbation time marker on x-axis (vertical line only, no label)
    if T_perturb is not None and T_perturb >= t_start:
        ax.axvline(x=T_perturb, color='#3574B1', linestyle='--', alpha=0.8, linewidth=1.5)
    
    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Time Step', fontsize=fontsize)
    ax.set_ylabel('Protein', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize-2)
    
    # Set x-axis ticks: adaptive based on range
    x_max = t_start + T_steps - 1
    x_range = x_max - t_start + 1
    if x_range <= 12:
        # Small range (zoomed view): ticks from t_start, every 1 step
        x_ticks = np.arange(t_start, x_max + 1, 1)
        ax.set_xticks(x_ticks)
        ax.set_xlim(t_start - 0.5, x_max + 0.5)
    else:
        # Large range (full view): ticks from 0, every 10 steps
        x_ticks = np.arange(0, x_max + 1, 10)
        ax.set_xticks(x_ticks)
    
    if title:
        ax.set_title(title, fontsize=fontsize+2, fontweight='bold')
    
    handles, labels = ax.get_legend_handles_labels()
    order_map = {"TF": 0, "Epigenetics": 1, "Others": 2}
    pairs = sorted(zip(handles, labels), key=lambda x: order_map.get(x[1], 99))
    if pairs:
        h_sorted, l_sorted = zip(*pairs)
        ax.legend(h_sorted, l_sorted, fontsize=legend_fontsize, loc='upper right', frameon=True)


def plot_trajectory_comparison(baseline_df: pd.DataFrame, perturbed_df: pd.DataFrame, 
                                response_df: pd.DataFrame, t_start: int = 1):
    """
    Generate trajectory comparison plots for each seed.
    
    Per Output Spec section 9.1:
    - baseline vs perturbed
    - collapse vs steady
    
    Args:
        baseline_df: DataFrame with baseline trajectory data
        perturbed_df: DataFrame with perturbed trajectory data
        response_df: DataFrame with response metrics
        t_start: Starting time step index (default 1, meaning t=1,2,3,...)
    
    Output: plots/seed{N}_trajectory.png (per-seed plots)
    """
    print("\n--- Generating trajectory comparison plots ---")
    
    seeds = response_df['seed'].unique()
    
    for seed in seeds:
        seed_baseline = baseline_df[baseline_df['seed'] == seed]
        seed_perturbed = perturbed_df[perturbed_df['seed'] == seed]
        seed_info = response_df[response_df['seed'] == seed].iloc[0]
        
        regime = seed_info['regime']
        perturb_type = seed_info['perturbation_type']
        target_gene = int(seed_info['target_gene'])
        T_perturb = int(seed_info['T_perturb'])
        
        T_max = int(seed_baseline['time'].max())
        G = int(seed_baseline['gene_id'].max()) + 1
        
        # Reshape to (T+1, G) matrices
        P_baseline = np.zeros((T_max + 1, G))
        P_perturbed = np.zeros((T_max + 1, G))
        
        for _, row in seed_baseline.iterrows():
            t, g = int(row['time']), int(row['gene_id'])
            P_baseline[t, g] = row['P']
        
        for _, row in seed_perturbed.iterrows():
            t, g = int(row['time']), int(row['gene_id'])
            P_perturbed[t, g] = row['P']
        
        # Create figure with 2x2 panels
        fig, axes = plt.subplots(2, 2, figsize=(14, 10))
        
        # Row 1: Full trajectory (keep original)
        # Panel A: Baseline trajectory (full)
        plot_single_panel(axes[0, 0], P_baseline, 'Baseline Trajectory (Full)', fontsize=12, 
                         legend_fontsize=9, t_start=t_start, T_perturb=T_perturb)
        
        # Panel B: Perturbed trajectory (full)
        plot_single_panel(axes[0, 1], P_perturbed, 'Perturbed Trajectory (Full)', fontsize=12,
                         legend_fontsize=9, t_start=t_start, T_perturb=T_perturb)
        
        # Row 2: Zoomed view around T_perturb (3 steps before, 6 steps after = 10 total)
        zoom_pre = 3   # 3 steps before T_perturb
        zoom_post = 6  # 6 steps after T_perturb
        t_zoom_start = max(t_start, T_perturb - zoom_pre)
        t_zoom_end = min(T_max, T_perturb + zoom_post)
        
        # Panel C: Baseline trajectory (zoomed)
        plot_single_panel(axes[1, 0], P_baseline, f'Baseline (T={t_zoom_start}~{t_zoom_end})', fontsize=12, 
                         legend_fontsize=9, t_start=t_zoom_start, t_end=t_zoom_end, T_perturb=T_perturb)
        
        # Panel D: Perturbed trajectory (zoomed)
        plot_single_panel(axes[1, 1], P_perturbed, f'Perturbed (T={t_zoom_start}~{t_zoom_end})', fontsize=12,
                         legend_fontsize=9, t_start=t_zoom_start, t_end=t_zoom_end, T_perturb=T_perturb)
        
        regime_color = '#E04532' if regime == 'collapse' else '#4C7CB8'
        plt.suptitle(f'Seed {seed} ({regime.upper()}) — {perturb_type} on TF{target_gene}\nT_perturb={T_perturb}', 
                     fontsize=14, fontweight='bold', color=regime_color, y=1)
        plt.tight_layout()
        
        save_path = os.path.join(PLOTS_DIR, f'seed{seed}_trajectory.png')
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {save_path}")


def plot_metrics_comparison(response_df: pd.DataFrame):
    """
    Generate metrics distribution plots.
    
    Per Output Spec section 9.2:
    - plots/metrics.png
    
    Show delta_TF (boxplot), recovery (bar plot), amplification (boxplot)
    """
    print("\n--- Generating metrics distribution plots ---")
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    colors = {'collapse': '#E04532', 'steady': '#4C7CB8'}
    
    # Panel 1: delta_TF (boxplot)
    ax = axes[0]
    collapse_data = response_df[response_df['regime'] == 'collapse']['delta_TF'].values
    steady_data = response_df[response_df['regime'] == 'steady']['delta_TF'].values
    
    bp = ax.boxplot([collapse_data, steady_data], 
                    tick_labels=['Collapse', 'Steady'],
                    patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['collapse'])
    bp['boxes'][1].set_facecolor(colors['steady'])
    
    ax.set_ylabel('delta_TF', fontsize=12)
    ax.set_title('delta_TF Distribution', fontsize=13, fontweight='bold')
    ax.tick_params(axis='x', labelsize=12)
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Panel 2: recovery rate (bar plot)
    ax = axes[1]
    collapse_recovery = response_df[response_df['regime'] == 'collapse']['recovery'].mean()
    steady_recovery = response_df[response_df['regime'] == 'steady']['recovery'].mean()
    
    x = np.arange(2)
    recovery_rates = [collapse_recovery, steady_recovery]
    bar_colors = [colors['collapse'], colors['steady']]
    
    bars = ax.bar(x, recovery_rates, width=0.4, color=bar_colors, alpha=0.8, edgecolor='black', linewidth=1.5)
    
    ax.set_xticks(x)
    ax.set_xticklabels(['Collapse', 'Steady'], fontsize=12)
    ax.set_ylabel('Recovery Rate', fontsize=12)
    ax.set_title('Recovery Rate', fontsize=13, fontweight='bold')
    ax.set_ylim(0, 1.1)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar, rate in zip(bars, recovery_rates):
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height,
                f'{rate:.2%}',
                ha='center', va='bottom', fontsize=11, fontweight='bold')
    
    # Panel 3: amplification (boxplot)
    ax = axes[2]
    collapse_amp = response_df[response_df['regime'] == 'collapse']['amplification'].dropna().values
    steady_amp = response_df[response_df['regime'] == 'steady']['amplification'].dropna().values
    
    bp = ax.boxplot([collapse_amp, steady_amp],
                    tick_labels=['Collapse', 'Steady'],
                    patch_artist=True)
    bp['boxes'][0].set_facecolor(colors['collapse'])
    bp['boxes'][1].set_facecolor(colors['steady'])
    
    ax.set_ylabel('Amplification', fontsize=12)
    ax.set_title('Amplification Distribution', fontsize=13, fontweight='bold')
    ax.tick_params(axis='x', labelsize=12)
    ax.grid(True, alpha=0.3)
    
    plt.suptitle('Run 006: Response Metrics Comparison', fontsize=15, fontweight='bold', y=1)
    plt.tight_layout()
    
    save_path = os.path.join(PLOTS_DIR, 'metrics.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def plot_merged_trajectory(baseline_df: pd.DataFrame, perturbed_df: pd.DataFrame,
                           response_df: pd.DataFrame, t_start: int = 1):
    """
    Generate a merged trajectory plot with all seeds in a 4x4 grid.
    
    Layout: 4 rows x 4 columns
    - Columns 0-1 (Left): Collapse seeds (baseline + perturbed)
    - Columns 2-3 (Right): Steady seeds (baseline + perturbed)
    """
    print("\n--- Generating merged trajectory plot ---")
    
    # Separate collapse and steady seeds
    collapse_seeds = response_df[response_df['regime'] == 'collapse']['seed'].tolist()
    steady_seeds = response_df[response_df['regime'] == 'steady']['seed'].tolist()
    
    # Sort seeds
    collapse_seeds = sorted(collapse_seeds)
    steady_seeds = sorted(steady_seeds)
    
    # 增加整体维度，保持图表清晰
    fig = plt.figure(figsize=(24, 13)) 
    
    # 💡 恢复为 1行2列 的外层宏观网格
    gs_outer = gridspec.GridSpec(1, 2, figure=fig, wspace=0.18)
    
    # 左右两区的内层网格 (4行2列)
    gs_left = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=gs_outer[0], hspace=0.5, wspace=0.18)
    gs_right = gridspec.GridSpecFromSubplotSpec(4, 2, subplot_spec=gs_outer[1], hspace=0.5, wspace=0.18)
    
    # ==========================================
    # 💡 终极定位方案：利用子图自身的相对坐标 (transAxes)
    # ==========================================
    # 左半区透明画布
    ax_left_title = fig.add_subplot(gs_outer[0])
    ax_left_title.axis('off') 
    # x=0.5 保证水平绝对居中。y=1.05 代表放置在包围盒上方 3% 的位置，确保完美悬浮在子图之上且不留大片空白
    ax_left_title.text(0.5, 1.05, 'COLLAPSE', transform=ax_left_title.transAxes, 
                       fontsize=18, fontweight='bold', color='#E04532', ha='center', va='bottom')
    
    # 右半区透明画布
    ax_right_title = fig.add_subplot(gs_outer[1])
    ax_right_title.axis('off')
    ax_right_title.text(0.5, 1.05, 'STEADY', transform=ax_right_title.transAxes, 
                        fontsize=18, fontweight='bold', color='#4C7CB8', ha='center', va='bottom')
    # ==========================================
    
    # Pre-create a 4x4 axes matrix to maintain compatibility with existing loop code
    axes = np.empty((4, 4), dtype=object)
    for r in range(4):
        axes[r, 0] = fig.add_subplot(gs_left[r, 0])
        axes[r, 1] = fig.add_subplot(gs_left[r, 1])
        axes[r, 2] = fig.add_subplot(gs_right[r, 0])
        axes[r, 3] = fig.add_subplot(gs_right[r, 1])

    # Plot collapse seeds on the left (columns 0-1)
    for row_idx, seed in enumerate(collapse_seeds):
        seed_baseline = baseline_df[baseline_df['seed'] == seed]
        seed_perturbed = perturbed_df[perturbed_df['seed'] == seed]
        seed_info = response_df[response_df['seed'] == seed].iloc[0]
        
        perturb_type = seed_info['perturbation_type']
        target_gene = int(seed_info['target_gene'])
        T_perturb = int(seed_info['T_perturb'])
        
        T_max = int(seed_baseline['time'].max())
        G = int(seed_baseline['gene_id'].max()) + 1
        
        # Reshape to (T+1, G) matrices
        P_baseline = np.zeros((T_max + 1, G))
        P_perturbed = np.zeros((T_max + 1, G))
        
        for _, row in seed_baseline.iterrows():
            t, g = int(row['time']), int(row['gene_id'])
            P_baseline[t, g] = row['P']
            
        for _, row in seed_perturbed.iterrows():
            t, g = int(row['time']), int(row['gene_id'])
            P_perturbed[t, g] = row['P']
            
        # Plot baseline (column 0)
        ax = axes[row_idx, 0]
        title = f'Seed {seed} - Baseline'
        plot_single_panel(ax, P_baseline, title, fontsize=11, legend_fontsize=8, 
                          t_start=t_start, T_perturb=T_perturb)
        ax.title.set_color('#E04532')
        
        # Plot perturbed (column 1)
        ax = axes[row_idx, 1]
        title_pert = f'Seed {seed} - Perturbed\n{perturb_type} on TF{target_gene}, T={T_perturb}'
        plot_single_panel(ax, P_perturbed, title_pert, fontsize=11, legend_fontsize=8,
                          t_start=t_start, T_perturb=T_perturb)
        ax.title.set_color('#E04532')
    
    # Plot steady seeds on the right (columns 2-3)
    for row_idx, seed in enumerate(steady_seeds):
        seed_baseline = baseline_df[baseline_df['seed'] == seed]
        seed_perturbed = perturbed_df[perturbed_df['seed'] == seed]
        seed_info = response_df[response_df['seed'] == seed].iloc[0]
        
        perturb_type = seed_info['perturbation_type']
        target_gene = int(seed_info['target_gene'])
        T_perturb = int(seed_info['T_perturb'])
        
        T_max = int(seed_baseline['time'].max())
        G = int(seed_baseline['gene_id'].max()) + 1
        
        # Reshape to (T+1, G) matrices
        P_baseline = np.zeros((T_max + 1, G))
        P_perturbed = np.zeros((T_max + 1, G))
        
        for _, row in seed_baseline.iterrows():
            t, g = int(row['time']), int(row['gene_id'])
            P_baseline[t, g] = row['P']
            
        for _, row in seed_perturbed.iterrows():
            t, g = int(row['time']), int(row['gene_id'])
            P_perturbed[t, g] = row['P']
            
        # Plot baseline (column 2)
        ax = axes[row_idx, 2]
        title = f'Seed {seed} - Baseline'
        plot_single_panel(ax, P_baseline, title, fontsize=11, legend_fontsize=8, 
                          t_start=t_start, T_perturb=T_perturb)
        ax.title.set_color('#4C7CB8')
        
        # Plot perturbed (column 3)
        ax = axes[row_idx, 3]
        title_pert = f'Seed {seed} - Perturbed\n{perturb_type} on TF{target_gene}, T={T_perturb}'
        plot_single_panel(ax, P_perturbed, title_pert, fontsize=11, legend_fontsize=8,
                          t_start=t_start, T_perturb=T_perturb)
        ax.title.set_color('#4C7CB8')

    # 💡 留出顶部的空间 (rect=[..., 0.93]) 并放置主标题 (y=0.98)
    plt.suptitle('Trajectory Comparison', fontsize=24, fontweight='bold', y=0.98)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning) 
        plt.tight_layout(rect=[0, 0, 1, 0.93]) 
    
    save_path = os.path.join(PLOTS_DIR, 'trajectory.png')
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {save_path}")


def generate_all_visualizations(t_start: int = 1):
    """
    Generate all visualizations for run_006.
    
    Args:
        t_start: Starting time step index (default 1, meaning t=1,2,3,...)
                  Set to 0 to include t=0 (initial state)
    
    Outputs:
    1. plots/seed{N}_trajectory.png (per-seed)
    2. plots/merged_trajectory.png (4x4 grid for all seeds)
    3. plots/metrics.png (delta_TF, recovery, amplification)
    """
    print("\n" + "=" * 60)
    print(f"Generating Visualizations (t_start={t_start})")
    print("=" * 60)
    
    # Load data
    baseline_df = pd.read_csv(BASELINE_TSV_PATH, sep='\t')
    perturbed_df = pd.read_csv(PERTURBED_TSV_PATH, sep='\t')
    response_df = pd.read_csv(os.path.join(DATA_DIR, 'response_metrics', 'response_summary.tsv'), sep='\t')
    
    # Generate visualizations
    # 1. Per-seed trajectory plots
    plot_trajectory_comparison(baseline_df, perturbed_df, response_df, t_start=t_start)
    
    # 2. Merged trajectory plot (4x4 grid)
    plot_merged_trajectory(baseline_df, perturbed_df, response_df, t_start=t_start)
    
    # 3. Metrics comparison (delta_TF, recovery, amplification)
    plot_metrics_comparison(response_df)
    
    print("\nAll visualizations generated successfully!")


def main():
    print("=" * 60)
    print("Run 006: External Perturbation Response Analysis")
    print("=" * 60)

    setup_directories()

    proxy_metrics_path = '/home/zhanghl/projects/ddc_github/runs/run_005_no_P0_proj/results/gain_proxy_analysis/proxy_metrics.tsv'
    collapse_seeds, steady_seeds = load_selected_seeds(proxy_metrics_path, n_per_regime=4)

    print(f"\nTarget Collapse seeds: {collapse_seeds}")
    print(f"Target Steady seeds:   {steady_seeds}")

    collapse_records = []
    steady_records = []

    for seed in collapse_seeds:
        print(f"\n--- Processing collapse seed {seed} ---")
        record = run_single_sample(seed, 'collapse')
        collapse_records.append(record)
        print(f"  delta_TF={record['delta_TF']:.4f}, recovery={record['recovery']}, amp={record['amplification']:.4f}")

    for seed in steady_seeds:
        print(f"\n--- Processing steady seed {seed} ---")
        record = run_single_sample(seed, 'steady')
        steady_records.append(record)
        print(f"  delta_TF={record['delta_TF']:.4f}, recovery={record['recovery']}, amp={record['amplification']:.4f}")

    df_collapse = pd.DataFrame(collapse_records)
    df_steady = pd.DataFrame(steady_records)
    df_all = pd.concat([df_collapse, df_steady], ignore_index=True)

    response_summary_path = os.path.join(DATA_DIR, 'response_metrics', 'response_summary.tsv')
    df_all.to_csv(response_summary_path, sep='\t', index=False)
    print(f"\nResponse summary saved to: {response_summary_path}")

    df_agg = pd.DataFrame({
        'regime': ['collapse', 'steady'],
        'mean_delta_TF': [df_collapse['delta_TF'].mean(), df_steady['delta_TF'].mean()],
        'mean_amplification': [df_collapse['amplification'].mean(), df_steady['amplification'].mean()],
        'recovery_rate': [df_collapse['recovery'].mean(), df_steady['recovery'].mean()],
        # 'n_samples': [len(df_collapse), len(df_steady)]
    })
    agg_path = os.path.join(RESULTS_DIR, 'aggregated', 'regime_response.tsv')
    df_agg.to_csv(agg_path, sep='\t', index=False)
    print(f"Aggregated results saved to: {agg_path}")

    # Generate visualizations
    generate_all_visualizations()

    print("\n" + "=" * 60)
    print("Run 006 completed successfully!")
    print("=" * 60)


if __name__ == '__main__':
    main()