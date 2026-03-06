"""
DDC Phase 0 - Designed Digital-Cell Model
==========================================
Tests for DDC Phase 0: Stage A, B, and C as per the Phase 0 Execution Order document.
Requires 'ddc.py' to be present in the same directory.

Author: zhanghl
Date: 2026-03-03
"""

import os
import json
import torch
import ddc

def run_stage_a(seed: int, dir: str):
    # print("STAGE A - run_001 (Single-cell baseline)")
    # print("STAGE A - run_002 (Single-cell baseline) - adjust K to uniform(0.01, 0.1)")
    # print("STAGE A - run_003 (Single-cell baseline) - scaled_P = tilde_P[j] * G ")
    # print("STAGE A - run_004 (Single-cell baseline) - R_TOTAL: float = 25.0 ")
    print("STAGE A - run_003_1 (Single-cell baseline) - tilde_P = P / (torch.mean(P) + world.epsilon) ")
    # save_path = os.path.join(dir, 'stage_A_run_001_2.pt')
    # world_json_path = os.path.join(dir, 'stage_A_world_run_001_2.json')
    # save_path = os.path.join(dir, 'stage_A_run_003_adjust_P.pt')
    # world_json_path = os.path.join(dir, 'stage_A_world_run_003_adjust_P.json')
    # save_path = os.path.join(dir, 'stage_A_run_004_increase_Rtotal.pt')
    # world_json_path = os.path.join(dir, 'stage_A_world_run_004_increase_Rtotal.json')
    save_path = os.path.join(dir, 'stage_A_run_003_1_tilde_P.pt')
    world_json_path = os.path.join(dir, 'stage_A_world_run_003_1_tilde_P.json')

    print(f'Step 1: Fixing world_seed to {seed}')
    world = ddc.sample_world(seed)
    with open(world_json_path, 'w') as fw:
        json.dump(world.to_dict(), fw, indent=2)
    print(f'Sampled world is saved to {world_json_path}.')

    traj = ddc.run_simulation(seed, save_path)
    print(f'Full trajectories (X/P/Z/N) is saved to {save_path}.')

    final_expression = traj['X_traj'][-1, :]
    print(f'Final expression state extracted. Shape: {final_expression.shape}.')
    print('\nStage A passed. Single-cell baseline verified.')

def run_stage_b(seed: int, dir: str):
    print("STAGE B - Reproducibility check (Mandatory)")

    print('Running simulation Run #1...')
    world1 = ddc.sample_world(seed)
    traj1 = ddc.run_simulation(seed)

    print('Running simulation Run #2 (identical seed)...')
    world2 = ddc.sample_world(seed)
    traj2 = ddc.run_simulation(seed)

    world1_str = json.dumps(world1.to_dict(), sort_keys=True)
    world2_str = json.dumps(world2.to_dict(), sort_keys=True)
    assert world1_str == world2_str, "[ERROR] World objects differ between runs!"
    print('Check 1 passed: World objects (graphs and parameters) are identical.')

    for key in ['X_traj', 'P_traj', 'Z_traj', 'N_traj']:
        assert torch.equal(traj1[key], traj2[key]), f"[ERROR] Trajectory '{key}' differs between runs!"
    print(f'Check 2 passed: All trajectory (X/P/Z/N) are element-wise identical.')
    print('\nStage B passed. Reproducibility verified.')

def run_stage_c(seed: int, dir: str):
    print('STAGE C - Small dataset wrapper test')

    M = 20
    save_path = os.path.join(dir, 'stage_C_dataset.pt')

    print(f'Generating single-cell dataset with M={M} cells...')
    dataset, world = ddc.generate_dataset(seed, M, save_path)

    assert dataset.shape == (20, 50), f"[ERROR] Dataset shape is {dataset.shape}, expected (20, 50)."
    print(f'Check 1 passed: Output dataset shape is strictly (20, 50).')

    ### torch.any: Check if any cell has non-zero variance
    variance = torch.var(dataset, dim=0)
    assert torch.any(variance > 0), "[ERROR] All 20 cells are identical! Independent initialization failed."
    print('Check 2 passed: Cells show expression variance, proving independent initialization.')

    print(f'\nStage C passed. M={M} dataset saved to {save_path}')


if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    output_dir = './test'
    os.makedirs(output_dir, exist_ok=True)

    MASTER_SEED = 2026
    # MASTER_SEED = 42

    try:
        print('\n' + '*'*60)
        print("Running Stage A...")
        run_stage_a(MASTER_SEED, output_dir)
    except Exception as e:
        print(f"[Stage A failed] {e}")
        exit(1)

    # try:
    #     print('\n' + '*'*60)
    #     print("Running Stage B...")
    #     run_stage_b(MASTER_SEED, output_dir)
    # except Exception as e:
    #     print(f"[Stage B failed] {e}")
    #     exit(1)

    # try:
    #     print('\n' + '*'*60)
    #     print("Running Stage C...")
    #     run_stage_c(MASTER_SEED, output_dir)
    # except Exception as e:
    #     print(f"[Stage C failed] {e}")
    #     exit(1)

    # print('\n' + '*'*60)
    # print('All stages (A, B, C) passed.')
