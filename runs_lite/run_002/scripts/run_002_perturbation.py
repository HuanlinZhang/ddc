"""Perturbation analysis for run_002: transient single-gene knockdown at t=500.

For each selected world/cell, the gene with the highest expression at t=500
in the original trajectory is knocked down (set to 0) for one time step.
"""

import json
import os
import sys
import torch

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT, 'src'))

import ddc_model_a as dma

RUN_DIR = os.path.join(ROOT, 'runs_lite', 'run_002')
METADATA_DIR = os.path.join(RUN_DIR, 'world_metadata')
TRAJ_DIR = os.path.join(RUN_DIR, 'trajectories')
OUT_DIR = os.path.join(RUN_DIR, 'perturbations')
os.makedirs(OUT_DIR, exist_ok=True)
SUMMARY_PATH = os.path.join(RUN_DIR, 'analysis', 'perturbation_summary.json')

CONVERGENCE_WINDOW = 50
EPSILON = 1e-4
SLOW_CONVERGENCE_THRESHOLD = 200
DIVERGENCE_THRESHOLD = 1e3
COLLAPSE_THRESHOLD = 1e-3

SELECTION = [
    ('stress_balanced_t005',                        'Type D', 0),
    ('chen_stress_balanced_t008',                   'Type C', 0),
    ('chen_stress_repression_biased_t001',          'Type C', 0),
    ('stress_repression_biased_t006',               'Type B', 1),
    ('chen_moderate_balanced_t007',                 'Type B', 0),
    ('stress_repression_biased_t007',               'Type E', 0),
    ('stress_balanced_t000',                        'Type E', 0),
    ('chen_stress_balanced_t005',                   'Type A', 0),
    ('chen_stress_repression_biased_t007',          'Type A', 0),
]


def detect_equilibrium(X_traj):
    t_steps = X_traj.shape[0] - 1
    for t in range(t_steps - CONVERGENCE_WINDOW + 1):
        all_consecutive = True
        for w in range(CONVERGENCE_WINDOW):
            diff = float(torch.norm(X_traj[t + w + 1] - X_traj[t + w]).item())
            if diff >= EPSILON:
                all_consecutive = False
                break
        if all_consecutive:
            return True, t
    return False, -1


def detect_oscillation(X_traj, converged, conv_time):
    if converged:
        end = conv_time if conv_time >= 0 else X_traj.shape[0] - 1
    else:
        end = X_traj.shape[0] - 1
    if end < 50:
        return 'none'
    N = X_traj.shape[1]
    last_half = X_traj[end // 2:end + 1]
    maxs = torch.max(last_half, dim=0).values
    mins = torch.min(last_half, dim=0).values
    ranges = maxs - mins
    n_osc = int((ranges > 0.05).sum().item())
    if n_osc < 1:
        return 'none'
    if n_osc >= N * 0.3:
        return 'sustained'
    return 'damped'


def classify(result_dict):
    eq_conv = result_dict['equilibrium_converged']
    eq_time = result_dict['equilibrium_time']
    osc = result_dict['oscillation_type']
    divergent = result_dict['divergent']
    collapsed = result_dict['collapsed']

    if divergent:
        return 'Type E'
    if collapsed:
        return 'Type F'
    if osc == 'sustained':
        return 'Type D'
    if osc == 'damped':
        return 'Type C'
    if eq_conv:
        if eq_time >= 0 and eq_time <= SLOW_CONVERGENCE_THRESHOLD:
            return 'Type A'
        return 'Type B'
    return 'Type G'


def main():
    summary = []

    for wid, original_type, cell_idx in SELECTION:
        cell_label = f'cell{cell_idx:02d}'
        print(f'\n=== {wid} / {cell_label} (original: {original_type}) ===')

        meta_path = os.path.join(METADATA_DIR, f'{wid}.json')
        with open(meta_path) as f:
            meta = json.load(f)

        w = dma.World(meta['seed'])
        w.from_dict(meta)

        traj_path = os.path.join(TRAJ_DIR, f'{wid}_{cell_label}.pt')
        orig_traj = torch.load(traj_path, weights_only=False)
        X_at_500 = orig_traj['X_traj'][500]
        knockdown_gene = int(torch.argmax(X_at_500).item())
        max_expr_before = float(X_at_500[knockdown_gene].item())
        print(f'  Highest expr at t=500: gene {knockdown_gene} = {max_expr_before:.4f}')

        cell_seed = meta['seed'] + cell_idx + 1
        X0 = dma.sample_initial_state(cell_seed)

        result = dma.simulate_single_cell(
            w, X0, dma.T_SIM,
            intervention_time=500,
            intervention_config={'knockdown': [knockdown_gene]},
        )

        X_traj = result['X_traj']

        eq_conv, eq_time = detect_equilibrium(X_traj)
        osc_type = detect_oscillation(X_traj, eq_conv, eq_time)
        max_expr = float(X_traj.max().item())
        final_mean = float(X_traj[-1].mean().item())
        divergent = max_expr >= DIVERGENCE_THRESHOLD
        collapsed = final_mean < COLLAPSE_THRESHOLD

        result_dict = {
            'equilibrium_converged': eq_conv,
            'equilibrium_time': eq_time,
            'oscillation_type': osc_type,
            'divergent': divergent,
            'collapsed': collapsed,
        }
        new_type = classify(result_dict)

        max_after_500 = float(X_traj[501:].max().item())
        recovered = abs(max_after_500 - max_expr_before) < 0.05 if not divergent else False

        out_path = os.path.join(OUT_DIR, f'{wid}_{cell_label}_perturb.pt')
        torch.save({
            'X_traj': result['X_traj'],
            'clip_count': result['clip_count'],
            'total_clips': result['total_clips'],
            'world_seed': meta['seed'],
            'cell_seed': cell_seed,
            'world_id': wid,
            'original_attractor': original_type,
            'knockdown_gene': knockdown_gene,
            'knockdown_gene_expr_at_500': max_expr_before,
            'knockdown_time': 500,
            'perturbation_type': 'knockdown',
        }, out_path)
        print(f'  Attractor: {original_type} -> {new_type}')
        print(f'  Saved: {os.path.basename(out_path)}')

        summary.append({
            'world_id': wid,
            'cell': cell_label,
            'original_attractor': original_type,
            'perturbed_attractor': new_type,
            'attractor_changed': new_type != original_type,
            'knockdown_gene': knockdown_gene,
            'knockdown_gene_expr_at_500': max_expr_before,
            'equilibrium_converged': eq_conv,
            'equilibrium_time': eq_time,
            'oscillation_type': osc_type,
            'divergent': divergent,
            'collapsed': collapsed,
            'max_expression': max_expr,
            'final_mean_expression': final_mean,
        })

    os.makedirs(os.path.dirname(SUMMARY_PATH), exist_ok=True)
    with open(SUMMARY_PATH, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f'\nSummary saved to {SUMMARY_PATH}')

    print('\n=== Summary ===')
    unchanged = sum(1 for s in summary if not s['attractor_changed'])
    changed = sum(1 for s in summary if s['attractor_changed'])
    print(f'  Unchanged: {unchanged}/{len(summary)}')
    print(f'  Changed:   {changed}/{len(summary)}')
    for s in summary:
        marker = '  ***' if s['attractor_changed'] else ''
        print(f'  {s["world_id"]:45s} {s["cell"]}  {s["original_attractor"]} -> {s["perturbed_attractor"]}{marker}')


if __name__ == '__main__':
    main()
