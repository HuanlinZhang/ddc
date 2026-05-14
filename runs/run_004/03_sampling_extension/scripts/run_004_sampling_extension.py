"""
Run 004 - Sampling Extension Script
=====================================

This script extends the run_004 analysis by testing whether the collapse
vs steady-state pattern observed in run_004 remains stable under larger
world sampling (~100 worlds).

Based on 03_Sampling_Extension.md specification.

Key features:
    - Generate 100 worlds using sample_world() with world_seed = 0..99
    - For each world, use cell_seed = world_seed + 1 per spec
    - Run Phase 0 simulation with T=200
    - Classify each world as collapse or steady based on final N_active
    - Save all world and trajectory data to test_convergence/v1_2_run_004_100seeds
    - Generate regime distribution statistics and visualizations

Output structure:
    - test_convergence/v1_2_run_004_100seeds/  (full data)
        - world_seed_{i}.json
        - seed_{i}_traj.tsv
    - run_004/03_sampling_extension/data/       (required outputs only)
        - world_regime_distribution.tsv
    - run_004/03_sampling_extension/results/
        - regime_summary.tsv
    - run_004/03_sampling_extension/plots/
        - regime_distribution.png
        - N_active_TF_vs_regime.png

Author: zhanghl
Date: 2026-04-13
"""

import os
import sys
import torch
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')
from ddc import TF_GENES, RBP_GENES, KINASE_GENES, PHOSPHATASE_GENES, EPI_GENES, CELL_CYCLE_GENES, APOPTOSIS_GENES, BACKGROUND_GENES, World, sample_world, sample_initial_state, simulate_single_cell

OUTPUT_BASE = '/home/zhanghl/projects/ddc_github/test_convergence/v1_2_run_004_100seeds'
SCRIPT_DIR = '/home/zhanghl/projects/ddc_github/runs/run_004/03_sampling_extension/scripts'
DATA_DIR = '/home/zhanghl/projects/ddc_github/runs/run_004/03_sampling_extension/data'
RESULTS_DIR = '/home/zhanghl/projects/ddc_github/runs/run_004/03_sampling_extension/results'
PLOTS_DIR = '/home/zhanghl/projects/ddc_github/runs/run_004/03_sampling_extension/plots'

os.makedirs(OUTPUT_BASE, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)

N_WORLDS = 100
T = 200
THRESHOLD = 0.1


def generate_worlds_and_simulate(n_worlds: int):
    """
    Generate n_worlds and run Phase 0 simulation for each.
    Save world and trajectory data to OUTPUT_BASE.
    """
    results = []

    for world_seed in range(n_worlds):
        print(f"  Processing world {world_seed + 1}/{n_worlds}...", end='\r')

        world = sample_world(world_seed)
        cell_seed = world_seed + 1
        X0, P0, Z0, N0 = sample_initial_state(cell_seed, world)

        traj = simulate_single_cell(world, X0.clone(), P0.clone(), Z0.clone(), N0, T)

        X_traj = traj['X_traj'].numpy()
        P_traj = traj['P_traj'].numpy()
        Z_traj = traj['Z_traj'].numpy()
        time = np.arange(T + 1)

        X_final = X_traj[-1, :]
        P_final = P_traj[-1, :]

        N_active = np.sum(X_final > THRESHOLD)

        tf_mask = np.array([i in TF_GENES for i in range(len(X_final))])
        tf_expression = X_final[tf_mask]
        N_active_TF = np.sum(tf_expression > THRESHOLD)

        mean_P_terminal = np.mean(P_final)

        regime = "steady" if N_active > 1 else "collapse"

        results.append({
            'seed': world_seed,
            'regime': regime,
            'N_active': int(N_active),
            'N_active_TF': int(N_active_TF),
            'mean_P_terminal': float(mean_P_terminal)
        })

        world_file = os.path.join(OUTPUT_BASE, f'world_seed_{world_seed}.pt')
        torch.save(world.to_dict(), world_file)

        traj_file = os.path.join(OUTPUT_BASE, f'seed_{world_seed}_traj.pt')
        torch.save({
            'X_traj': torch.tensor(X_traj),
            'P_traj': torch.tensor(P_traj),
            'Z_traj': torch.tensor(Z_traj),
            'N_traj': torch.tensor(N0) if isinstance(N0, (int, float)) else N0,
            'world': world.to_dict(),
        }, traj_file)

    print(f"\n  Completed {n_worlds} worlds processing")

    return results


def plot_P_traj_by_category(world_seed: int, P_traj: np.ndarray, time: np.ndarray, regime: str = None):
    """
    Plot P_traj for all genes, grouped by category.
    Style follows run_006_analysis.py plot_single_panel.
    """
    TF_COLOR = plt.cm.Set1(0)
    EPI_COLOR = plt.cm.Set1(3)
    OTHER_COLOR = 'gray'

    fig, ax = plt.subplots(figsize=(10, 4.5))

    added_labels = set()
    for gene in range(P_traj.shape[1]):
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
        ax.plot(time, P_traj[:, gene], color=color, alpha=alpha, lw=lw, zorder=zorder, label=label)

    ax.grid(True, alpha=0.3)
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Time Step', fontsize=10)
    ax.set_ylabel('Protein', fontsize=10)
    title = f'Seed {world_seed}: P_i(t) Protein Trajectory'
    if regime:
        title += f' ({regime.capitalize()})'
    ax.set_title(title, fontsize=11, fontweight='bold')

    from matplotlib.lines import Line2D
    gene_legend = [
        Line2D([0], [0], color=TF_COLOR, lw=1.8, label='TF'),
        Line2D([0], [0], color=EPI_COLOR, lw=1.8, label='Epigenetics'),
        Line2D([0], [0], color=OTHER_COLOR, lw=0.8, label='Others')
    ]
    ax.legend(handles=gene_legend, loc='upper right', fontsize=8, ncol=1)

    plt.tight_layout()

    plot_file = os.path.join(OUTPUT_BASE, f'seed_{world_seed}_P_traj.png')
    plt.savefig(plot_file, dpi=100, bbox_inches='tight')
    plt.close()


def compute_regime_statistics(results: list):
    """Compute regime distribution statistics."""
    df = pd.DataFrame(results)

    N_total = len(df)
    N_collapse = len(df[df['regime'] == 'collapse'])
    N_steady = len(df[df['regime'] == 'steady'])

    collapse_rate = N_collapse / N_total if N_total > 0 else 0
    steady_rate = N_steady / N_total if N_total > 0 else 0

    return {
        'N_total': N_total,
        'N_collapse': N_collapse,
        'N_steady': N_steady,
        'collapse_rate': collapse_rate,
        'steady_rate': steady_rate
    }


# def plot_combined(df: pd.DataFrame, stats: dict):
#     """Plot combined figure: regime distribution + N_active_TF vs regime."""
#     fig = plt.figure(figsize=(12, 5))

#     ax_bar = fig.add_subplot(1, 2, 1)
#     regimes = ['Collapse', 'Steady']
#     counts = [stats['N_collapse'], stats['N_steady']]
#     colors = ['#E04532', '#4C7CB8']

#     bars = ax_bar.bar(regimes, counts, color=colors, edgecolor='black', linewidth=1.2, width=0.6)

#     for bar, count in zip(bars, counts):
#         height = bar.get_height()
#         ax_bar.annotate(f'{count}\n({count/stats["N_total"]*100:.1f}%)',
#                         xy=(bar.get_x() + bar.get_width() / 2, height),
#                         xytext=(0, 3),
#                         textcoords="offset points",
#                         ha='center', va='bottom', fontsize=10, fontweight='bold')

#     ax_bar.set_xlabel('Regime', fontsize=11)
#     ax_bar.set_ylabel('Number of Worlds', fontsize=11)
#     ax_bar.set_title(f'Regime Distribution\n(N={stats["N_total"]})', fontsize=12, fontweight='bold')
#     ax_bar.grid(True, alpha=0.3, axis='y')
#     ax_bar.set_ylim(0, max(counts) * 1.25)

#     collapse_data = df[df['regime'] == 'collapse']['N_active_TF'].values
#     steady_data = df[df['regime'] == 'steady']['N_active_TF'].values

#     ax_box = fig.add_subplot(1, 2, 2)
#     bp = ax_box.boxplot([collapse_data, steady_data],
#                         tick_labels=['Collapse', 'Steady'],
#                         patch_artist=True,
#                         widths=0.5)
#     bp['boxes'][0].set_facecolor('#E04532')
#     bp['boxes'][1].set_facecolor('#4C7CB8')

#     ax_box.set_xlabel('Regime', fontsize=11)
#     ax_box.set_ylabel('N_active_TF', fontsize=11)
#     ax_box.set_title('N_active_TF by Regime', fontsize=12, fontweight='bold')
#     ax_box.grid(True, alpha=0.3, axis='y')

#     collapse_mean = np.mean(collapse_data) if len(collapse_data) > 0 else 0
#     steady_mean = np.mean(steady_data) if len(steady_data) > 0 else 0

#     y_max = max(collapse_data.max() if len(collapse_data) > 0 else 1,
#                 steady_data.max() if len(steady_data) > 0 else 1)
#     y_upper = y_max * 1.15 if y_max > 0 else 1.15
#     ax_box.set_ylim(-0.5, y_upper)

#     ax_box.annotate(f'mean: {collapse_mean:.2f}', xy=(1, collapse_mean),
#                     xytext=(1.15, min(collapse_mean * 1.1, y_upper * 0.95)), fontsize=9, va='top')
#     ax_box.annotate(f'mean: {steady_mean:.2f}', xy=(2, steady_mean),
#                     xytext=(2.15, min(steady_mean * 1.1, y_upper * 0.95)), fontsize=9, va='top')

#     plt.tight_layout()

#     plot_file = os.path.join(PLOTS_DIR, 'sampling_extension_combined.png')
#     plt.savefig(plot_file, dpi=150, bbox_inches='tight')
#     plt.close()
#     print(f"  Saved: {plot_file}")

def plot_combined(df: pd.DataFrame, stats: dict):
    """Plot combined figure: regime distribution + N_active_TF vs regime."""
    fig = plt.figure(figsize=(12, 5))

    # ==========================================
    # Subplot 1: Regime Distribution (保持不变)
    # ==========================================
    ax_bar = fig.add_subplot(1, 2, 1)
    regimes = ['Collapse', 'Steady']
    counts = [stats['N_collapse'], stats['N_steady']]
    colors = ['#E04532', '#4C7CB8']

    bars = ax_bar.bar(regimes, counts, color=colors, edgecolor='black', linewidth=1.2, width=0.6)

    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax_bar.annotate(f'{count}\n({count/stats["N_total"]*100:.1f}%)',
                        xy=(bar.get_x() + bar.get_width() / 2, height),
                        xytext=(0, 3),
                        textcoords="offset points",
                        ha='center', va='bottom', fontsize=10, fontweight='bold')

    ax_bar.set_xlabel('Regime', fontsize=11)
    ax_bar.set_ylabel('Number of Worlds', fontsize=11)
    ax_bar.set_title(f'Regime Distribution (N={stats["N_total"]})', fontsize=12, fontweight='bold')
    ax_bar.grid(True, alpha=0.3, axis='y')
    ax_bar.set_ylim(0, max(counts) * 1.25)

    # ==========================================
    # Subplot 2: Boxplot + Swarmplot 叠加 (全新修改)
    # ==========================================
    ax_box = fig.add_subplot(1, 2, 2)
    
    # 规范化列名，方便 seaborn 绘图
    df_plot = df.copy()
    df_plot['Regime_Cap'] = df_plot['regime'].str.capitalize()
    order = ['Collapse', 'Steady']
    palette = {'Collapse': '#E04532', 'Steady': '#4C7CB8'}

    # 1. 绘制底层的箱线图（调高透明度，隐藏异常值以免与散点重复）
    sns.boxplot(
        data=df_plot, 
        x='Regime_Cap', 
        y='N_active_TF', 
        order=order, 
        palette=palette,
        ax=ax_box,
        width=0.4,
        boxprops=dict(alpha=0.3, edgecolor='black'),
        whiskerprops=dict(alpha=0.5, color='black'),
        capprops=dict(alpha=0.5, color='black'),
        medianprops=dict(color='orange', linewidth=2),
        showfliers=False, # 重点：关闭原生的离群点绘制
        zorder=1
    )

    # 2. 叠加蜂群图（Swarmplot），展现每一个真实的数据点
    sns.swarmplot(
        data=df_plot, 
        x='Regime_Cap', 
        y='N_active_TF', 
        order=order, 
        palette=palette,
        size=6,           # 点的大小
        alpha=0.8,        # 点的透明度
        edgecolor='white',# 白色描边让点更清晰
        linewidth=1,
        ax=ax_box,
        zorder=2
    )

    ax_box.set_xlabel('Regime', fontsize=11)
    ax_box.set_ylabel('N_active_TF', fontsize=11)
    ax_box.set_title('N_active_TF by Regime', fontsize=12, fontweight='bold')
    ax_box.grid(True, alpha=0.3, axis='y')

    # 计算均值并添加标注
    collapse_data = df_plot[df_plot['Regime_Cap'] == 'Collapse']['N_active_TF'].values
    steady_data = df_plot[df_plot['Regime_Cap'] == 'Steady']['N_active_TF'].values
    collapse_mean = np.mean(collapse_data) if len(collapse_data) > 0 else 0
    steady_mean = np.mean(steady_data) if len(steady_data) > 0 else 0

    y_max = max(collapse_data.max() if len(collapse_data) > 0 else 1,
                steady_data.max() if len(steady_data) > 0 else 1)
    y_upper = y_max * 1.15 if y_max > 0 else 1.15
    ax_box.set_ylim(-0.5, y_upper)

    # 注意：Seaborn 分类轴的 x 坐标是从 0 开始的整数 (0 为 Collapse, 1 为 Steady)
    ax_box.annotate(f'mean: {collapse_mean:.2f}', xy=(0, collapse_mean),
                    xytext=(0.15, min(collapse_mean * 1.1, y_upper * 0.95)), fontsize=9, va='top')
    ax_box.annotate(f'mean: {steady_mean:.2f}', xy=(1, steady_mean),
                    xytext=(1.15, min(steady_mean * 1.1, y_upper * 0.95)), fontsize=9, va='top')

    plt.tight_layout()

    plot_file = os.path.join(PLOTS_DIR, 'sampling_extension_combined.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_file}")


def main():
    print("=" * 60)
    print("Run 004 - Sampling Extension")
    print("=" * 60)

    sample_traj_file = os.path.join(OUTPUT_BASE, f'seed_0_traj.pt')

    if os.path.exists(sample_traj_file):
        print("\n" + "-" * 40)
        print("Loading existing data")
        print("-" * 40)
        df = pd.read_csv(os.path.join(DATA_DIR, 'world_regime_distribution.tsv'), sep='\t')
        results = df.to_dict('records')
        stats = compute_regime_statistics(results)
        print(f"  Loaded {len(results)} worlds from cached trajectory data")
    else:
        print("\n" + "-" * 40)
        print(f"Generating {N_WORLDS} worlds with T={T}")
        print("-" * 40)

        results = generate_worlds_and_simulate(N_WORLDS)

        df = pd.DataFrame(results)
        df.to_csv(os.path.join(DATA_DIR, 'world_regime_distribution.tsv'), sep='\t', index=False)
        print(f"  Saved: {os.path.join(DATA_DIR, 'world_regime_distribution.tsv')}")

        stats = compute_regime_statistics(results)

    print(f"\n  Results:")
    print(f"    Total: {stats['N_total']}, Collapse: {stats['N_collapse']}, Steady: {stats['N_steady']}")
    print(f"    Collapse Rate: {stats['collapse_rate']*100:.1f}%")

    print("\n" + "-" * 40)
    print("Step 2: Save regime summary")
    print("-" * 40)

    regime_summary = pd.DataFrame([stats])
    summary_file = os.path.join(RESULTS_DIR, 'regime_summary.tsv')
    regime_summary.to_csv(summary_file, sep='\t', index=False)
    print(f"  Saved: {summary_file}")

    print("\n" + "-" * 40)
    print("Step 3: Generate P_traj plots from saved data")
    print("-" * 40)

    for world_seed in range(N_WORLDS):
        traj_file = os.path.join(OUTPUT_BASE, f'seed_{world_seed}_traj.pt')
        plot_file = os.path.join(OUTPUT_BASE, f'seed_{world_seed}_P_traj.png')
        if os.path.exists(traj_file) and not os.path.exists(plot_file):
            data = torch.load(traj_file)
            P_traj = data['P_traj'].numpy()
            time = np.arange(1, P_traj.shape[0])
            seed_info = next((r for r in results if r['seed'] == world_seed), None)
            regime = seed_info['regime'] if seed_info else None
            plot_P_traj_by_category(world_seed, P_traj[1:], time, regime)
            print(f"  Generated P_traj plot for seed {world_seed}", end='\r')
    print(f"  P_traj plots generation complete")

    print("\n" + "-" * 40)
    print("Step 4: Generate summary visualizations")
    print("-" * 40)

    plot_combined(df, stats)

    print("\n" + "=" * 60)
    print("Run 004 Sampling Extension completed successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()
