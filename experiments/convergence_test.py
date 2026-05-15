"""
Convergence Test Script
=======================
Test if different random seeds produce convergent (non-zero stable) states.
Uses run_simulation from ddc.py.
"""

import os
import sys
import torch
import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')

OUTPUT_DIR = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_1'
os.makedirs(OUTPUT_DIR, exist_ok=True)

from ddc import run_simulation

np.random.seed(12345)
SEEDS = np.random.randint(0, 10000, size=10).tolist()

SEEDS.append(2026)

print(f"Testing seeds: {SEEDS}")

results = []

for idx, seed in enumerate(SEEDS):
    print(f"\n[{idx+1}/{len(SEEDS)}] Testing seed {seed}...")

    traj_path = os.path.join(OUTPUT_DIR, f'seed_{seed}_traj.pt')
    traj = run_simulation(seed, save_path=traj_path)

    saved_data = torch.load(traj_path)
    world_dict = saved_data.get('world', {})

    world_json_path = os.path.join(OUTPUT_DIR, f'seed_{seed}_world.json')
    with open(world_json_path, 'w') as f:
        json.dump(world_dict, f, indent=2)

    print(f"  Saved: {traj_path}")
    print(f"  Saved: {world_json_path}")

    X_data = traj['X_traj'].numpy()
    P_data = traj['P_traj'].numpy()
    Z_data = traj['Z_traj'].numpy()

    gene_to_macro = world_dict.get('gene_annotation', {}).get('to_macro', {})
    gene_to_micro = world_dict.get('gene_annotation', {}).get('to_micro', {})

    T_steps, G_count = X_data.shape
    t_axis = np.arange(T_steps)

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax1, ax2 = axes[0, 0], axes[0, 1]
    ax3, ax4 = axes[1, 0], axes[1, 1]

    TF_GENES = list(range(0, 6))
    EPI_GENES = list(range(17, 20))
    
    added_labels = set()
    
    for i in range(G_count):
        if i in TF_GENES:
            current_color = plt.cm.Set1(0)
            label = "TF"
            alpha, linewidth, zorder = 0.85, 1.8, 2
        elif i in EPI_GENES:
            current_color = plt.cm.Set1(3)
            label = "Epigenetics"
            alpha, linewidth, zorder = 0.85, 1.8, 2
        else:
            current_color = 'gray'
            label = "Others"
            alpha, linewidth, zorder = 0.5, 0.8, 1
        
        plot_label = label if label not in added_labels else None
        if label:
            added_labels.add(label)
        
        ax1.plot(t_axis, X_data[:, i], color=current_color, alpha=alpha, linewidth=linewidth, zorder=zorder, label=plot_label)
        ax2.plot(t_axis, P_data[:, i], color=current_color, alpha=alpha, linewidth=linewidth, zorder=zorder)
        ax3.plot(t_axis, Z_data[:, i], color=current_color, alpha=alpha, linewidth=linewidth, zorder=zorder)

    ax1.set_xlabel('Time (T)')
    ax1.set_ylabel('mRNA Expression')
    ax1.set_title(f'X_traj (mRNA) - Seed {seed}')
    ax1.grid(True, linestyle='--', alpha=0.6)

    ax2.set_xlabel('Time (T)')
    ax2.set_ylabel('Protein Level')
    ax2.set_title(f'P_traj (Protein) - Seed {seed}')
    ax2.grid(True, linestyle='--', alpha=0.6)

    ax3.set_xlabel('Time (T)')
    ax3.set_ylabel('Chromatin State')
    ax3.set_title(f'Z_traj (Chromatin) - Seed {seed}')
    ax3.grid(True, linestyle='--', alpha=0.6)

    ax4.axis('off')

    handles, labels = ax1.get_legend_handles_labels()
    
    label_order = {"TF": 0, "Epigenetics": 1, "others": 2}
    sorted_pairs = sorted(zip(handles, labels), key=lambda x: label_order.get(x[1], 99))
    handles_sorted, labels_sorted = zip(*sorted_pairs) if sorted_pairs else ([], [])
    
    ax4.legend(handles_sorted, labels_sorted, loc='center left', title="Gene Categories", fontsize=14, title_fontsize=16)

    plt.tight_layout()

    fig_path = os.path.join(OUTPUT_DIR, f'seed_{seed}.png')
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    plt.close()

    X_final = X_data[-1, :].mean()
    P_final = P_data[-1, :].mean()
    convergence = (X_final > 0.01 and P_final > 0.01)

    ### 提取world参数
    world_json_path = os.path.join(OUTPUT_DIR, f'seed_{seed}_world.json')
    with open(world_json_path, 'r') as f:
        world_data = json.load(f)
    
    result = {
        'seed': seed,
        'X_final_mean': X_final,
        'P_final_mean': P_final,
        'convergence': convergence,
    }
    
    results.append(result)

    print(f"  Final X mean: {X_final:.4f}, Final P mean: {P_final:.4f}, Convergence: {convergence}")

print(f"\nAll results saved to: {OUTPUT_DIR}")

df_results = pd.DataFrame(results)
df_results = df_results.set_index('seed')
df_results['X_final_mean'] = df_results['X_final_mean'].round(4)
df_results['P_final_mean'] = df_results['P_final_mean'].round(4)

print("\n" + "=" * 60)
print("Summary Results")
print("=" * 60)
print(df_results.to_string())

seeds_conv = [s for s in df_results.index if df_results.loc[s, 'convergence'] == True]
seeds_collapsed = [s for s in df_results.index if df_results.loc[s, 'convergence'] == False]

print("\n" + "=" * 60)
print("Seeds by Convergence Status")
print("=" * 60)
print(f"Converged: {seeds_conv}")
print(f"Collapsed: {seeds_collapsed}")

seed_convergence = df_results['convergence'].to_dict()

world_params_list = []

for seed in SEEDS:
    convergence = seed_convergence.get(seed, False)
    
    world_json_path = os.path.join(OUTPUT_DIR, f'seed_{seed}_world.json')
    with open(world_json_path, 'r') as f:
        world_data = json.load(f)
    
    params = world_data.get('parameters', {})
    n = params.get('n', [2.0] * 50)
    gamma = params.get('gamma', [1.0] * 50)
    a_ij = params.get('a_ij', {})
    beta_ij = params.get('beta_ij', {})
    
    alpha = params.get('alpha', [])
    rho = params.get('rho', [])
    K = params.get('K', [])
    delta_x = params.get('delta_x', [])
    delta_p = params.get('delta_p', [])
    
    for gene_idx in range(len(alpha)):
        P_regulators = [int(j) for j in a_ij.get(str(gene_idx), {}).keys()]
        P_a_ij_dict = a_ij.get(str(gene_idx), {})
        P_a_ij_str = str({int(k): v for k, v in P_a_ij_dict.items()})
        
        E_regulators = [int(j) for j in beta_ij.get(str(gene_idx), {}).keys()]
        E_beta_ij_dict = beta_ij.get(str(gene_idx), {})
        E_beta_ij_str = str({int(k): v for k, v in E_beta_ij_dict.items()})
        
        world_params_list.append({
            'seed': seed,
            'convergence': convergence,
            'gene_index': gene_idx,
            'alpha': alpha[gene_idx],
            'rho': rho[gene_idx],
            'K': K[gene_idx],
            'n': n[gene_idx] if isinstance(n, list) else n,
            'delta_x': delta_x[gene_idx],
            'delta_p': delta_p[gene_idx],
            'gamma': gamma[gene_idx] if isinstance(gamma, list) else gamma,
            'P_regulators': str(P_regulators),
            'P_a_ij': P_a_ij_str,
            'E_regulators': str(E_regulators),
            'E_beta_ij': E_beta_ij_str,
        })

df_world_params = pd.DataFrame(world_params_list)

for seed in SEEDS:
    seed_df = df_world_params[df_world_params['seed'] == seed]
    seed_df.to_csv(os.path.join(OUTPUT_DIR, f'world_parameters_seed{seed}.tsv'), sep='\t', index=False)

df_world_params = df_world_params.sort_values(by=['seed', 'gene_index']).reset_index(drop=True)
df_world_params.to_csv(os.path.join(OUTPUT_DIR, 'world_parameters_11_seeds.tsv'), sep='\t', index=False)
print(f"\nWorld parameters saved to: {OUTPUT_DIR}/world_parameters_11_seeds.tsv")
