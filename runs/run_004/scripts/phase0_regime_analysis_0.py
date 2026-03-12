"""
Phase 0 Regime Analysis Script -- 0
===============================

Based on the analysis requirements from "Phase0 Regime Analysis and Next Steps" document.
Computes summary statistics for distinguishing collapse vs steady-state worlds.

Author: zhanghl
Date: 2026-03-11
"""

import os
import sys
import torch
import numpy as np
import pandas as pd
import json

sys.path.insert(0, '/home/zhanghl/projects/ddc_github/src')

from ddc import TF_GENES

# Configuration
OUTPUT_DIR = '/home/zhanghl/projects/ddc_github/test_convergence/v1_1_1'
THRESHOLD = 0.1  # ε = 0.1 for active gene definition

def extract_seeds_from_directory():
    """Extract seeds from trajectory files in the output directory."""
    seeds = []
    
    if not os.path.exists(OUTPUT_DIR):
        print(f"Warning: Output directory {OUTPUT_DIR} does not exist")
        return seeds
    
    for filename in os.listdir(OUTPUT_DIR):
        if filename.startswith('seed_') and filename.endswith('_traj.pt'):
            # Extract seed number from filename: seed_123_traj.pt → 123
            try:
                seed_str = filename.replace('seed_', '').replace('_traj.pt', '')
                seed = int(seed_str)
                seeds.append(seed)
            except ValueError:
                print(f"Warning: Could not extract seed from filename: {filename}")
    
    return sorted(seeds)

# Extract seeds from existing files
SEEDS = extract_seeds_from_directory()

def compute_regime_statistics(seed):
    """
    Compute regime statistics for a given seed.
    
    Returns:
        dict with keys: seed, N_active, mean_X_T, max_X_T, N_active_TF, mean_Z_T
    """
    
    # Load trajectory data
    traj_path = os.path.join(OUTPUT_DIR, f'seed_{seed}_traj.pt')
    if not os.path.exists(traj_path):
        print(f"Warning: Trajectory file not found for seed {seed}")
        return None
    
    saved_data = torch.load(traj_path)
    
    # Extract final time step data
    X_traj = saved_data['X_traj']
    Z_traj = saved_data['Z_traj']
    
    # Final time step (T = 200)
    X_final = X_traj[-1, :].numpy()  # Shape: (G,)
    Z_final = Z_traj[-1, :].numpy()  # Shape: (G,)
    
    # 1. Active gene count: N_active = #{X_i(T) > ε}
    N_active = np.sum(X_final > THRESHOLD)
    
    # 2. Active TF gene count: N_active_TF = #{TF genes with X_i(T) > ε}
    tf_mask = np.array([i in TF_GENES for i in range(len(X_final))])
    tf_expression = X_final[tf_mask]
    N_active_TF = np.sum(tf_expression > THRESHOLD)
    
    # 3. Mean terminal expression: mean_X_T = mean_i X_i(T)
    mean_X_T = np.mean(X_final)
    
    # 4. Maximum expression: max_X_T = max_i X_i(T)
    max_X_T = np.max(X_final)
    
    # 5. Mean chromatin openness: mean_Z_T = mean_i Z_i(T)
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
    """
    Classify world as collapse or steady-state based on convergence data.
    Uses convergence column from world_parameters_11_seeds.tsv.
    """
    convergence_file = os.path.join(OUTPUT_DIR, 'world_parameters_11_seeds.tsv')
    
    df_convergence = pd.read_csv(convergence_file, sep='\t')
    # Get convergence status for this seed
    seed_data = df_convergence[df_convergence['seed'] == seed]
    if not seed_data.empty:
        convergence_status = seed_data['convergence'].iloc[0]
        return "steady-state" if convergence_status else "collapse"


def main():
    """Main analysis function."""
    
    print("Phase 0 Regime Analysis")
    print("=" * 60)
    
    # Compute statistics for all seeds
    results = []
    for seed in SEEDS:
        print(f"Analyzing seed {seed}...")
        stats = compute_regime_statistics(seed)
        if stats:
            stats['regime'] = classify_regime(seed)
            results.append(stats)
    
    # Create summary table
    df_summary = pd.DataFrame(results)
    
    # Format output
    df_summary['mean_X_T'] = df_summary['mean_X_T'].round(4)
    df_summary['max_X_T'] = df_summary['max_X_T'].round(4)
    df_summary['mean_Z_T'] = df_summary['mean_Z_T'].round(4)
    
    # Sort by seed
    df_summary = df_summary.sort_values('seed').reset_index(drop=True)
    
    print("\n" + "=" * 60)
    print("Summary Statistics Table")
    print("=" * 60)
    print(df_summary.to_string(index=False))
    
    # Separate seeds by regime
    seeds_collapse = df_summary[df_summary['regime'] == 'collapse']['seed'].tolist()
    seeds_steady = df_summary[df_summary['regime'] == 'steady-state']['seed'].tolist()
    
    print("\n" + "=" * 60)
    print("Regime Classification")
    print("=" * 60)
    print(f"Collapse worlds: {seeds_collapse}")
    print(f"Steady-state worlds: {seeds_steady}")
    
    # Save results
    output_path = os.path.join(OUTPUT_DIR, 'phase0_regime_analysis.tsv')
    df_summary.to_csv(output_path, sep='\t', index=False)
    print(f"\nResults saved to: {output_path}")
    
    return df_summary

if __name__ == "__main__":
    main()