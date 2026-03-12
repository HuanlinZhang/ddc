"""
Phase 0 Regime Analysis Script
===============================

This script performs three analysis tasks to understand the mechanism
behind collapse vs steady-state regime formation in Phase0 DDC model.

Task 1: Parameter Comparison
    - Classify worlds based on simulation trajectories
    - Compare parameter distributions between collapse and steady-state worlds
    
Task 2: TF Regulatory Topology Analysis
    - Analyze TF → TF network structure
    - Compute network metrics (degree, SCC)
    
Task 3: TF Expression Dynamics Analysis
    - Analyze TF gene expression trajectories
    - Identify key driver genes

Author: zhanghl
Date: 2026-03-11
"""

import os
import sys
import torch
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import networkx as nx
from ddc import TF_GENES

DATA_INPUT_DIR = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_1'
OUTPUT_DIR = '/home/zhanghl/projects/ddc_github/runs/run_004'
RESULTS_DIR = os.path.join(OUTPUT_DIR, 'results')
PLOTS_DIR = os.path.join(OUTPUT_DIR, 'plots')
DATA_DIR = os.path.join(OUTPUT_DIR, 'data')

os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PLOTS_DIR, exist_ok=True)
os.makedirs(DATA_DIR, exist_ok=True)

THRESHOLD = 0.1

def extract_seeds_from_directory():
    """Extract seeds from trajectory files in the input data directory."""
    seeds = []
    if not os.path.exists(DATA_INPUT_DIR):
        print(f"Warning: Input data directory {DATA_INPUT_DIR} does not exist")
        return seeds
    
    for filename in os.listdir(DATA_INPUT_DIR):
        if filename.startswith('seed_') and filename.endswith('_traj.pt'):
            try:
                seed_str = filename.replace('seed_', '').replace('_traj.pt', '')
                seed = int(seed_str)
                seeds.append(seed)
            except ValueError:
                print(f"Warning: Could not extract seed from filename: {filename}")
    
    return sorted(seeds)

def load_world_parameters(seed):
    """Load world parameters from JSON file."""
    world_path = os.path.join(DATA_INPUT_DIR, f'seed_{seed}_world.json')
    if not os.path.exists(world_path):
        return None
    
    with open(world_path, 'r') as f:
        world_data = json.load(f)
    
    params = world_data.get('parameters', {})
    return {
        'seed': seed,
        'alpha': params.get('alpha', []),
        'rho': params.get('rho', []),
        'K': params.get('K', []),
        'n': params.get('n', []),
        'delta_x': params.get('delta_x', []),
        'delta_p': params.get('delta_p', []),
        'a_ij': params.get('a_ij', {}),
        'beta_ij': params.get('beta_ij', {})
    }

def compute_regime_statistics(seed):
    """Compute regime statistics for a given seed."""
    traj_path = os.path.join(DATA_INPUT_DIR, f'seed_{seed}_traj.pt')
    if not os.path.exists(traj_path):
        return None
    
    saved_data = torch.load(traj_path)
    X_traj = saved_data['X_traj']
    Z_traj = saved_data['Z_traj']
    
    X_final = X_traj[-1, :].numpy()
    Z_final = Z_traj[-1, :].numpy()
    
    N_active = np.sum(X_final > THRESHOLD)
    
    tf_mask = np.array([i in TF_GENES for i in range(len(X_final))])
    tf_expression = X_final[tf_mask]
    N_active_TF = np.sum(tf_expression > THRESHOLD)
    
    mean_X_T = np.mean(X_final)
    max_X_T = np.max(X_final)
    mean_Z_T = np.mean(Z_final)
    
    return {
        'seed': seed,
        'N_active': N_active,
        'mean_X_T': mean_X_T,
        'max_X_T': max_X_T,
        'N_active_TF': N_active_TF,
        'mean_Z_T': mean_Z_T
    }

def classify_regime(seed):
    """Classify world as collapse or steady based on convergence data."""
    convergence_file = os.path.join(DATA_INPUT_DIR, 'world_parameters_11_seeds.tsv')
    
    if not os.path.exists(convergence_file):
        return "collapse"
    
    try:
        df_convergence = pd.read_csv(convergence_file, sep='\t')
        seed_data = df_convergence[df_convergence['seed'] == seed]
        if not seed_data.empty:
            convergence_status = seed_data['convergence'].iloc[0]
            return "steady" if convergence_status else "collapse"
        else:
            return "collapse"
    except Exception:
        return "collapse"

def task1_parameter_comparison():
    """
    Task 1: Parameter Comparison
    =============================
    Compare parameter distributions between collapse and steady-state worlds.
    
    Outputs:
        - results/world_regime_classification.tsv
        - results/parameter_summary.tsv
        - plots/parameter_distribution.png
    """
    print("\n" + "=" * 60)
    print("Task 1: Parameter Comparison")
    print("=" * 60)
    
    SEEDS = extract_seeds_from_directory()
    print(f"Found {len(SEEDS)} seeds: {SEEDS}")
    
    results = []
    for seed in SEEDS:
        stats = compute_regime_statistics(seed)
        if stats:
            stats['regime'] = classify_regime(seed)
            results.append(stats)
    
    df_summary = pd.DataFrame(results)
    df_summary = df_summary.sort_values('seed').reset_index(drop=True)
    
    df_regime = df_summary[['seed', 'regime']]
    
    regime_file = os.path.join(RESULTS_DIR, 'world_regime_classification.tsv')
    df_regime.to_csv(regime_file, sep='\t', index=False)
    print(f"Saved: {regime_file}")
    
    seeds_collapse = df_summary[df_summary['regime'] == 'collapse']['seed'].tolist()
    seeds_steady = df_summary[df_summary['regime'] == 'steady']['seed'].tolist()
    
    print(f"\nCollapse worlds: {seeds_collapse}")
    print(f"Steady worlds: {seeds_steady}")
    
    params_to_analyze = ['rho', 'delta_x', 'K']
    param_summary = []
    
    for param_name in params_to_analyze:
        collapse_values = []
        steady_values = []
        
        for seed in seeds_collapse:
            params = load_world_parameters(seed)
            if params and (param_name in params):
                collapse_values.extend(params[param_name])
        
        for seed in seeds_steady:
            params = load_world_parameters(seed)
            if params and (param_name in params):
                steady_values.extend(params[param_name])
        
        param_summary.append({
            'parameter': param_name,
            'mean_collapse': np.mean(collapse_values) if collapse_values else 0,
            'std_collapse': np.std(collapse_values) if collapse_values else 0,
            'mean_steady': np.mean(steady_values) if steady_values else 0,
            'std_steady': np.std(steady_values) if steady_values else 0,
            # 'n_collapse': len(collapse_values),
            # 'n_steady': len(steady_values)
        })
    
    df_param_summary = pd.DataFrame(param_summary)
    df_param_summary = df_param_summary.set_index('parameter').T
    df_param_summary = df_param_summary.reset_index().rename(columns={'index': 'statistic'})
    
    param_file = os.path.join(RESULTS_DIR, 'parameter_summary.tsv')
    df_param_summary.to_csv(param_file, sep='\t', index=False)
    print(f"Saved: {param_file}")
    
    print("\nParameter Summary:")
    print(df_param_summary.to_string(index=False))
    
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    
    for idx, param_name in enumerate(params_to_analyze):
        collapse_vals = []
        steady_vals = []
        
        for seed in seeds_collapse:
            params = load_world_parameters(seed)
            if params and param_name in params:
                collapse_vals.extend(params[param_name])
        
        for seed in seeds_steady:
            params = load_world_parameters(seed)
            if params and param_name in params:
                steady_vals.extend(params[param_name])
        
        ax_hist1 = axes[0, idx]
        ax_hist1.hist(collapse_vals, bins=50, label='Collapse', color='#E04532', edgecolor='black')
        ax_hist1.set_title(f'{param_name} - Collapse')
        ax_hist1.set_xlabel(param_name)
        ax_hist1.set_ylabel('Count')
        ax_hist1.grid(True, alpha=0.3)
        
        ax_hist2 = axes[1, idx]
        ax_hist2.hist(steady_vals, bins=50, label='Steady', color='#4C7CB8', edgecolor='black')
        ax_hist2.set_title(f'{param_name} - Steady')
        ax_hist2.set_xlabel(param_name)
        ax_hist2.set_ylabel('Count')
        ax_hist2.grid(True, alpha=0.3)
        
        ax_box = axes[2, idx]
        bp = ax_box.boxplot([collapse_vals, steady_vals], tick_labels=['Collapse', 'Steady'], patch_artist=True)
        bp['boxes'][0].set_facecolor('#E04532')
        bp['boxes'][1].set_facecolor('#4C7CB8')
        
        collapse_mean = np.mean(collapse_vals)
        steady_mean = np.mean(steady_vals)
        
        ax_box.annotate(f'mean:{collapse_mean:.4f}', xy=(1, collapse_mean), xytext=(1.08, collapse_mean),
                       fontsize=9, va='center')
        ax_box.annotate(f'mean:{steady_mean:.4f}', xy=(2, steady_mean), xytext=(2.08, steady_mean),
                       fontsize=9, va='center')
        
        ax_box.set_title(f'{param_name} Boxplot')
        ax_box.set_ylabel(param_name)
        ax_box.grid(True, alpha=0.3)
    
    plt.suptitle('Parameter Distribution: Collapse vs Steady Worlds', fontsize=16)
    plt.tight_layout()
    
    plot_file = os.path.join(PLOTS_DIR, 'parameter_distribution.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {plot_file}")
    
    return df_summary

def task2_step1_extract_tf_gene_list():
    """Step 1: Extract TF gene list"""
    tf_gene_data = []
    for gene_id in range(50):
        tf_gene_data.append({
            'gene_id': gene_id,
            'is_TF': gene_id in TF_GENES
        })
    
    df_tf_genes = pd.DataFrame(tf_gene_data)
    tf_gene_file = os.path.join(DATA_DIR, 'TF_gene_list.tsv')
    df_tf_genes.to_csv(tf_gene_file, sep='\t', index=False)
    print(f"  Saved: {tf_gene_file}")
    print(f"  Total genes: 50, TF genes: {len(TF_GENES)}")
    
    return tf_gene_file

def task2_step2_extract_tf_network_edges(SEEDS):
    """Step 2: Extract TF→TF edges from all worlds"""
    all_edges = []
    
    for seed in SEEDS:
        params = load_world_parameters(seed)
        if not params:
            continue
        
        a_ij = params.get('a_ij', {})
        
        # Extract TF→TF edges
        for source_str, targets in a_ij.items():
            source = int(source_str)
            if source not in TF_GENES:
                continue
            
            for target_str, weight in targets.items():
                target = int(target_str)
                if target in TF_GENES:
                    all_edges.append({
                        'seed': seed,
                        'tf_source': source,
                        'tf_target': target,
                        'weight': weight
                    })
    
    # Save to file
    df_edges = pd.DataFrame(all_edges)
    tf_edges_file = os.path.join(DATA_DIR, 'TF_network_edges.tsv')
    df_edges.to_csv(tf_edges_file, sep='\t', index=False)
    print(f"  Saved: {tf_edges_file}")
    print(f"  Total TF→TF edges: {len(all_edges)}")
    
    return df_edges

def task2_step3_compute_network_metrics(SEEDS):
    """Step 3: Compute network metrics for each world"""
    all_metrics = []
    
    for seed in SEEDS:
        params = load_world_parameters(seed)
        if not params:
            continue
        
        a_ij = params.get('a_ij', {})
        
        # Build TF network
        G = nx.DiGraph()
        G.add_nodes_from(TF_GENES)
        
        # Add TF→TF edges
        for source_str, targets in a_ij.items():
            source = int(source_str)
            if source not in TF_GENES:
                continue
            
            for target_str, weight in targets.items():
                target = int(target_str)
                if target in TF_GENES:
                    G.add_edge(source, target, weight=weight)
        
        # Compute degrees
        in_degrees = dict(G.in_degree())
        out_degrees = dict(G.out_degree())
        total_degrees = dict(G.degree())
        
        # Compute SCC
        scc = list(nx.strongly_connected_components(G))
        largest_scc = max(scc, key=len) if scc else set()
        scc_id_map = {}
        for idx, component in enumerate(scc):
            for node in component:
                scc_id_map[node] = idx
        
        # Detect cycles
        try:
            cycles = list(nx.simple_cycles(G))
            n_cycles = len(cycles)
        except:
            cycles = []
            n_cycles = len(list(nx.simple_cycles(G))) if G.number_of_edges() < 20 else 0
        
        # Collect metrics for each TF node
        for tf_node in TF_GENES:
            all_metrics.append({
                'seed': seed,
                'regime': classify_regime(seed),
                'tf_node': tf_node,
                'in_degree': in_degrees.get(tf_node, 0),
                'out_degree': out_degrees.get(tf_node, 0),
                'total_degree': total_degrees.get(tf_node, 0),
                'SCC_id': scc_id_map.get(tf_node, -1),
                'in_largest_SCC': tf_node in largest_scc,
                'largest_SCC_size': len(largest_scc),
                'n_cycles': n_cycles
            })
    
    # Save to file
    df_metrics = pd.DataFrame(all_metrics)
    metrics_file = os.path.join(RESULTS_DIR, 'TF_network_metrics.tsv')
    df_metrics.to_csv(metrics_file, sep='\t', index=False)
    print(f"  Saved: {metrics_file}")
    
    # Print summary
    for regime in ['collapse', 'steady']:
        regime_data = df_metrics[df_metrics['regime'] == regime]
        seeds_regime = regime_data['seed'].unique()
        print(f"\n  {regime.upper()} worlds ({len(seeds_regime)} seeds):")
        for seed in seeds_regime:
            seed_data = regime_data[regime_data['seed'] == seed]
            largest_scc = seed_data['largest_SCC_size'].iloc[0]
            n_cycles = seed_data['n_cycles'].iloc[0]
            avg_degree = seed_data['total_degree'].mean()
            print(f"    Seed {seed}: SCC_size={largest_scc}, Cycles={n_cycles}, Avg_degree={avg_degree:.2f}")
    
    return df_metrics

def task2_step4_visualize_tf_networks(SEEDS):
    """Step 4: Visualize TF networks for each world"""
    # Sort seeds by regime: collapse first, then steady
    seeds_collapse = []
    seeds_steady = []
    for seed in SEEDS:
        regime = classify_regime(seed)
        if regime == 'collapse':
            seeds_collapse.append(seed)
        else:
            seeds_steady.append(seed)
    
    n_collapse = len(seeds_collapse)
    n_steady = len(seeds_steady)
    n_cols = 4
    
    # Determine rows: collapse in first rows, steady in last row
    n_collapse_rows = (n_collapse + n_cols - 1) // n_cols
    n_rows = n_collapse_rows + 1  # +1 for steady row
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    axes = axes.flatten() if n_rows > 1 else [axes]
    
    # Helper function to plot a single network
    def plot_network(ax, seed):
        params = load_world_parameters(seed)
        if not params:
            return
        
        a_ij = params.get('a_ij', {})
        regime = classify_regime(seed)
        
        # Build TF network
        G = nx.DiGraph()
        G.add_nodes_from(TF_GENES)
        
        for source_str, targets in a_ij.items():
            source = int(source_str)
            if source not in TF_GENES:
                continue
            
            for target_str, weight in targets.items():
                target = int(target_str)
                if target in TF_GENES:
                    G.add_edge(source, target, weight=weight)
        
        # Compute metrics for visualization
        scc = list(nx.strongly_connected_components(G))
        largest_scc = max(scc, key=len) if scc else set()
        scc_size = len(largest_scc)
        
        try:
            cycles = list(nx.simple_cycles(G))
            n_cycles = len(cycles)
        except:
            n_cycles = 0
        
        total_degrees = dict(G.degree())
        
        # Layout
        try:
            pos = nx.circular_layout(G, scale=0.8)
        except:
            pos = nx.random_layout(G, seed=42)
        
        # Draw nodes
        node_colors = []
        node_sizes = []
        for node in G.nodes():
            if node in largest_scc:
                node_colors.append('#4C7CB8')  # In largest SCC - blue
            else:
                node_colors.append('#CCCCCC')  # Not in SCC - gray
            
            node_size = 300 + total_degrees.get(node, 0) * 100
            node_sizes.append(node_size)
        
        # Draw edges
        edge_colors = []
        for u, v in G.edges():
            weight = G[u][v].get('weight', 0)
            if weight > 0:
                edge_colors.append('#E04532')  # Activation - red
            else:
                edge_colors.append('#4C7CB8')  # Repression - blue
        
        # Draw network
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors, 
                              width=1.5, alpha=0.6, arrows=True,
                              arrowsize=10, connectionstyle="arc3,rad=0.1")
        
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                              node_size=node_sizes, alpha=0.8)
        
        nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, 
                               font_weight='bold')
        
        # Title
        ax.set_title(f'{regime.upper()} - Seed {seed}\n'
                    f'SCC: {scc_size} | Cycles: {n_cycles}', fontsize=10)
        ax.axis('off')
        ax.grid(False)
    
    # Plot collapse worlds first
    for idx, seed in enumerate(seeds_collapse):
        plot_network(axes[idx], seed)
    
    # Plot steady worlds in last row
    steady_start_idx = n_collapse_rows * n_cols
    for idx, seed in enumerate(seeds_steady):
        plot_network(axes[steady_start_idx + idx], seed)
    
    # Hide unused axes (between collapse and steady)
    for idx in range(n_collapse, steady_start_idx):
        axes[idx].axis('off')
    
    # Hide unused axes after steady
    for idx in range(steady_start_idx + n_steady, len(axes)):
        axes[idx].axis('off')
    
    plt.suptitle('TF Regulatory Network Topology\n(Collapse worlds in first rows, Steady worlds in last row)', 
                fontsize=16, y=1.02)
    plt.tight_layout()
    
    plot_file = os.path.join(PLOTS_DIR, 'TF_network_topology.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_file}")
    
    return plot_file

def task2_tf_topology_analysis():
    """
    Task 2: TF Regulatory Topology Analysis
    ========================================
    Analyze TF → TF network structure.
    
    Outputs:
        - data/TF_gene_list.tsv
        - data/TF_network_edges.tsv
        - results/TF_network_metrics.tsv
        - plots/TF_network_topology.png
    """
    print("\n" + "=" * 60)
    print("Task 2: TF Regulatory Topology Analysis")
    print("=" * 60)
    
    # Get seeds
    SEEDS = extract_seeds_from_directory()
    print(f"Found {len(SEEDS)} seeds: {SEEDS}")
    
    # Step 1: Extract TF gene list
    print("\nStep 1: Extracting TF gene list...")
    task2_step1_extract_tf_gene_list()
    
    # Step 2: Extract TF network edges
    print("\nStep 2: Extracting TF network edges...")
    task2_step2_extract_tf_network_edges(SEEDS)
    
    # Step 3: Compute network metrics
    print("\nStep 3: Computing network metrics...")
    task2_step3_compute_network_metrics(SEEDS)
    
    # Step 4: Visualize networks
    print("\nStep 4: Visualizing TF networks...")
    task2_step4_visualize_tf_networks(SEEDS)
    
    print("\nTask 2 completed successfully!")



def task3_step1_extract_tf_expression_trajectories(SEEDS):
    """
    Step 1: Extract TF gene expression trajectories from all worlds.
    
    Outputs:
        - data/TF_expression_trajectories.tsv
    """
    all_trajectories = []
    
    for seed in SEEDS:
        traj_path = os.path.join(DATA_INPUT_DIR, f'seed_{seed}_traj.pt')
        if not os.path.exists(traj_path):
            continue
        
        saved_data = torch.load(traj_path)
        P_traj = saved_data['P_traj']  # shape: (time_steps, n_genes)
        
        regime = classify_regime(seed)
        
        # Extract TF gene expression trajectories
        # P_traj shape: (time_steps, n_genes)
        # For each TF gene, get its expression over time
        for tf_gene_id in TF_GENES:
            tf_expression = P_traj[:, tf_gene_id].numpy()  # shape: (time_steps,)
            
            # Store trajectory data with time points
            for t, expr_val in enumerate(tf_expression):
                all_trajectories.append({
                    'seed': seed,
                    'regime': regime,
                    'gene_id': tf_gene_id,
                    'time': t,
                    'expression': expr_val
                })
    
    # Save to file
    df_trajectories = pd.DataFrame(all_trajectories)
    trajectory_file = os.path.join(DATA_DIR, 'TF_expression_trajectories.tsv')
    df_trajectories.to_csv(trajectory_file, sep='\t', index=False)
    print(f"  Saved: {trajectory_file}")
    print(f"  Total trajectory records: {len(all_trajectories)}")
    print(f"  TF genes tracked: {len(TF_GENES)}")
    print(f"  Seeds processed: {len(SEEDS)}")
    
    return df_trajectories

def task3_step2_visualize_tf_dynamics():
    """
    Step 2: Visualize TF expression dynamics.
    
    Layout: 3 rows x 4 cols
        - Row 0-1: Collapse worlds (7 seeds: 4 in row 0, 3 in row 1)
        - Row 2: Steady worlds (4 seeds)
    
    Outputs:
        - plots/TF_dynamics.png
    """
    # Read trajectory data
    trajectory_file = os.path.join(DATA_DIR, 'TF_expression_trajectories.tsv')
    if not os.path.exists(trajectory_file):
        print(f"  Error: {trajectory_file} not found. Run Step 1 first.")
        return None
    
    df = pd.read_csv(trajectory_file, sep='\t')
    
    # Group seeds by regime
    seeds_collapse = sorted(df[df['regime'] == 'collapse']['seed'].unique())
    seeds_steady = sorted(df[df['regime'] == 'steady']['seed'].unique())
    
    n_collapse = len(seeds_collapse)
    n_steady = len(seeds_steady)
    
    print(f"  Collapse worlds: {n_collapse} seeds")
    print(f"  Steady worlds: {n_steady} seeds")
    
    # Create figure: 3 rows x 4 cols
    n_rows = 3
    n_cols = 4
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 12))
    
    # Color palette for TF genes
    tf_colors = plt.cm.tab10(np.linspace(0, 1, len(TF_GENES)))
    
    # Helper function to plot single world
    def plot_single_world(ax, seed, regime):
        world_data = df[df['seed'] == seed].sort_values('time')
        
        for color_idx, gene_id in enumerate(sorted(TF_GENES)):
            gene_data = world_data[world_data['gene_id'] == gene_id]
            ax.plot(gene_data['time'], gene_data['expression'],
                   color=tf_colors[color_idx], linewidth=1.5,
                   label=f'TF{gene_id}')
        
        if regime == 'collapse':
            ax.set_title(f'COLLAPSE - Seed {seed}', fontsize=10, color='#E04532', fontweight='bold')
        else:
            ax.set_title(f'STEADY - Seed {seed}', fontsize=10, color='#4C7CB8', fontweight='bold')
        
        ax.set_xlabel('Time Step')
        ax.set_ylabel('Protein Expression')
        ax.legend(fontsize=7, loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
    
    # Plot collapse worlds in row 0 and row 1
    for idx, seed in enumerate(seeds_collapse):
        row = idx // n_cols
        col = idx % n_cols
        plot_single_world(axes[row, col], seed, 'collapse')
    
    # Hide unused axes in row 1 (for collapse)
    n_collapse_rows = (n_collapse + n_cols - 1) // n_cols
    for col in range(n_collapse % n_cols, n_cols):
        if n_collapse % n_cols != 0:  # Only if not exactly filled
            axes[n_collapse_rows - 1, col].axis('off')
    
    # Plot steady worlds in row 2
    for idx, seed in enumerate(seeds_steady):
        plot_single_world(axes[2, idx], seed, 'steady')
    
    # Hide unused axes in row 2 (for steady)
    for col in range(n_steady, n_cols):
        axes[2, col].axis('off')
    
    plt.suptitle('TF Expression Dynamics: Collapse vs Steady Worlds\n(Rows 1-2: Collapse, Row 3: Steady)',
                fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save plot
    plot_file = os.path.join(PLOTS_DIR, 'TF_dynamics.png')
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {plot_file}")
    
    return plot_file

def task3_tf_dynamics_analysis():
    """
    Task 3: TF Expression Dynamics Analysis
    ========================================
    Analyze TF gene expression trajectories.
    
    Outputs:
        - data/TF_expression_trajectories.tsv
        - plots/TF_dynamics.png
    """
    print("\n" + "=" * 60)
    print("Task 3: TF Expression Dynamics Analysis")
    print("=" * 60)
    
    # Get seeds
    SEEDS = extract_seeds_from_directory()
    print(f"Found {len(SEEDS)} seeds: {SEEDS}")
    
    # Step 1: Extract TF expression trajectories
    print("\nStep 1: Extracting TF expression trajectories...")
    task3_step1_extract_tf_expression_trajectories(SEEDS)
    
    # Step 2: Visualize TF dynamics
    print("\nStep 2: Visualizing TF dynamics...")
    task3_step2_visualize_tf_dynamics()
    
    print("\nTask 3 Step 1-2 completed successfully!")

def main():
    print("Phase 0 Regime Analysis")
    print("=" * 60)
    
    # df_summary = task1_parameter_comparison()
    
    # task2_tf_topology_analysis()
    
    task3_tf_dynamics_analysis()
    
    print("\n" + "=" * 60)
    print("Analysis Complete")
    print("=" * 60)

if __name__ == "__main__":
    main()
