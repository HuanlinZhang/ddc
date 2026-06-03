"""
Run 001 — Equilibrium Characterization
========================================

Phase A: fixed graph + parameter sweep

Steps:
  1. sample Level 0 worlds (fixed graph, parameter sweep)
  2. sample initial states per world
  3. run trajectories (T=200)
  4. characterize convergence / stability
  5. aggregate stability statistics
  6. parameter regime analysis

Per: docs/01_DDC_Lite_Curriculum/run_001/
"""
import copy
import os
import sys
import csv
import inspect
import json
import re
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
import ddc_lite as dl


N_WORLDS: int = 20
M_CELLS_PER_WORLD: int = 3
T_SIM: int = 200
BASE_SEED: int = 1000
EPS: float = 0.01
WINDOW: int = 10
COLLAPSE_THR: float = 1e-3
DIVERGENCE_THR: float = 1e3

RUN_DIR: str = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _parse_default_ranges():
    source = inspect.getsource(dl.sample_world)
    m_rho = re.search(r'\.rho\s*=.*?uniform_\(([^,]+),\s*([^,)]+)', source)
    m_K = re.search(r'\.K\s*=.*?uniform_\(([^,]+),\s*([^,)]+)', source)
    m_dx = re.search(r'\.delta_x\s*=.*?uniform_\(([^,]+),\s*([^,)]+)', source)
    m_dp = re.search(r'\.delta_p\s*=.*?uniform_\(([^,]+),\s*([^,)]+)', source)
    return {
        'rho': (float(m_rho.group(1)), float(m_rho.group(2))),
        'K': (float(m_K.group(1)), float(m_K.group(2))),
        'delta_x': (float(m_dx.group(1)), float(m_dx.group(2))),
        'delta_p': (float(m_dp.group(1)), float(m_dp.group(2))),
    }


_DDC_RANGES = _parse_default_ranges()

_RHO_MIN_DEFAULT = _DDC_RANGES['rho'][0]
_RHO_MAX_DEFAULT = _DDC_RANGES['rho'][1]
_K_MIN, _K_MAX = _DDC_RANGES['K']
_DX_MIN, _DX_MAX = _DDC_RANGES['delta_x']
_DP_MIN_DEFAULT = _DDC_RANGES['delta_p'][0]
_DP_MAX_DEFAULT = _DDC_RANGES['delta_p'][1]


def _resample_parameters(w: dl.World, seed: int,
                         rho_min: float = None, rho_max: float = None,
                         dp_min: float = None, dp_max: float = None) -> None:
    _rmin = rho_min if rho_min is not None else _RHO_MIN_DEFAULT
    _rmax = rho_max if rho_max is not None else _RHO_MAX_DEFAULT
    _dpmin = dp_min if dp_min is not None else _DP_MIN_DEFAULT
    _dpmax = dp_max if dp_max is not None else _DP_MAX_DEFAULT
    rng: torch.Generator = torch.Generator()
    rng.manual_seed(seed)
    w.rho = torch.empty(dl.G, dtype=dl.DTYPE).uniform_(_rmin, _rmax, generator=rng)
    w.K = torch.empty(dl.G, dtype=dl.DTYPE).uniform_(_K_MIN, _K_MAX, generator=rng)
    w.delta_x = torch.empty(dl.G, dtype=dl.DTYPE).uniform_(_DX_MIN, _DX_MAX, generator=rng)
    w.delta_p = torch.empty(dl.G, dtype=dl.DTYPE).uniform_(_dpmin, _dpmax, generator=rng)


def run_sweep(sweep_label: str,
              rho_min: float = None, rho_max: float = None,
              dp_min: float = None, dp_max: float = None):
    data_dir = os.path.join(RUN_DIR, 'data', sweep_label)
    worlds_dir = os.path.join(data_dir, 'worlds')
    traj_dir = os.path.join(data_dir, 'trajectories')
    metrics_dir = os.path.join(data_dir, 'metrics')
    results_dir = os.path.join(RUN_DIR, 'results', sweep_label)
    agg_dir = os.path.join(results_dir, 'aggregated')
    plots_dir = os.path.join(RUN_DIR, 'plots', sweep_label)

    for d in [worlds_dir, traj_dir, metrics_dir, agg_dir, plots_dir]:
        os.makedirs(d, exist_ok=True)

    world_seeds: list[int] = [BASE_SEED + i for i in range(N_WORLDS)]

    _rmin_print = rho_min if rho_min is not None else _RHO_MIN_DEFAULT
    _rmax_print = rho_max if rho_max is not None else _RHO_MAX_DEFAULT
    _dpmin_print = dp_min if dp_min is not None else _DP_MIN_DEFAULT
    _dpmax_print = dp_max if dp_max is not None else _DP_MAX_DEFAULT
    print(f'\n===== Sweep: {sweep_label} (rho {_rmin_print}~{_rmax_print}, δp {_dpmin_print}~{_dpmax_print}) =====')
    print(f'Step 1: sampling {N_WORLDS} worlds (fixed graph, parameter sweep) ...')
    base_world: dl.World = dl.sample_world(BASE_SEED)
    dl.export_ground_truth(base_world,
                           os.path.join(worlds_dir, f'world_{BASE_SEED}_base.json'))
    print(f'  base graph world saved (seed={BASE_SEED})')

    worlds: dict[int, dl.World] = {}
    for ws in world_seeds:
        w = copy.deepcopy(base_world)
        _resample_parameters(w, ws, rho_min, rho_max, dp_min, dp_max)
        worlds[ws] = w
        dl.export_ground_truth(w, os.path.join(worlds_dir, f'world_{ws}.json'))
    print(f'  saved {N_WORLDS} world jsons.')

    print(f'Step 2+3: simulating {N_WORLDS} x {M_CELLS_PER_WORLD} trajectories ...')
    all_metrics: list[dict] = []
    cell_idx: int = 0

    for ws in world_seeds:
        w = worlds[ws]
        for c in range(M_CELLS_PER_WORLD):
            cell_seed = BASE_SEED + N_WORLDS + cell_idx
            cell_idx += 1
            X0, P0 = dl.sample_initial_state(cell_seed, w)
            traj = dl.simulate_single_cell(w, X0, P0, T_SIM)

            tsv_path = os.path.join(traj_dir, f'trajectory_{ws}_{cell_seed}.tsv')
            _save_trajectory_tsv(traj, ws, cell_seed, tsv_path)

            metrics = dl.compute_stability_metrics(
                traj, eps=EPS, window=WINDOW,
                collapse_threshold=COLLAPSE_THR,
                divergence_threshold=DIVERGENCE_THR,
            )
            metrics['world_seed'] = ws
            metrics['cell_seed'] = cell_seed

            act_count = sum(
                1 for i in range(dl.G)
                if w.edge_signs[i][w.P_graph[i][0]] == dl.ACTIVATION
            )
            metrics['activation_count'] = act_count
            metrics['repression_count'] = dl.G - act_count
            metrics['repression_ratio'] = (dl.G - act_count) / dl.G

            all_metrics.append(metrics)

    print(f'  saved {cell_idx} trajectory TSVs.')

    print('Step 4: saving metrics ...')
    _save_metrics_tsv(all_metrics, metrics_dir)

    print('Step 5: aggregating statistics ...')
    per_world_metrics: list[dict] = []
    for ws in world_seeds:
        wm_list = [m for m in all_metrics if m['world_seed'] == ws]
        n_cells = len(wm_list)
        n_conv = sum(1 for m in wm_list if m['converged'])
        n_bounded = sum(1 for m in wm_list if m['bounded'])
        n_collapsed = sum(1 for m in wm_list if m['collapsed'])
        conv_times = [m['convergence_time'] for m in wm_list if m['converged']]
        max_expr = max(m['max_expression'] for m in wm_list)
        final_expr = sum(m['final_mean_expression'] for m in wm_list) / n_cells

        pw = {
            'world_seed': ws,
            'n_cells': n_cells,
            'n_converged': n_conv,
            'n_bounded': n_bounded,
            'n_collapsed': n_collapsed,
            'convergence_rate': n_conv / n_cells,
            'bounded_rate': n_bounded / n_cells,
            'collapse_rate': n_collapsed / n_cells,
            'world_converged': n_conv >= 2,
            'world_bounded': n_bounded == n_cells,
            'world_collapsed': n_collapsed >= 2,
            'converged': n_conv >= 2,
            'bounded': n_bounded == n_cells,
            'collapsed': n_collapsed >= 2,
            'convergence_time': (sum(conv_times) / len(conv_times)) if conv_times else -1,
            'max_expression': max_expr,
            'final_mean_expression': final_expr,
            'activation_count': wm_list[0]['activation_count'],
            'repression_count': wm_list[0]['repression_count'],
            'repression_ratio': wm_list[0]['repression_ratio'],
        }
        per_world_metrics.append(pw)

    agg = dl.aggregate_stability_stats(per_world_metrics)
    agg_path = os.path.join(agg_dir, 'stability_statistics.tsv')
    with open(agg_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(agg.keys())
        writer.writerow(agg.values())
    print(f'  convergence_rate={agg["convergence_rate"]:.3f}  '
          f'bounded_rate={agg["bounded_rate"]:.3f}  '
          f'collapse_rate={agg["collapse_rate"]:.3f}')

    print('Step 6: generating plots ...')
    _plot_all_trajectories(worlds, world_seeds, all_metrics, plots_dir, traj_dir)
    _plot_stability_distribution(per_world_metrics, plots_dir)
    _plot_parameter_regime_scan(per_world_metrics, plots_dir)
    _plot_parameter_convergence_analysis(worlds, per_world_metrics, plots_dir)

    print(f'\n  converged: {agg["n_converged"]}/{agg["total_worlds"]} '
          f'({agg["convergence_rate"]:.1%})')
    print(f'  bounded:   {agg["n_bounded"]}/{agg["total_worlds"]} '
          f'({agg["bounded_rate"]:.1%})')
    print(f'  collapsed: {agg["n_collapsed"]}/{agg["total_worlds"]} '
          f'({agg["collapse_rate"]:.1%})')

    print(f'Sweep {sweep_label} complete.\n')


def main():
    run_sweep('default')
    run_sweep('rho0.2-0.8', 0.2, 0.8)
    run_sweep('dp0.1-0.4', 0.5, 2.0, dp_min=0.10, dp_max=0.4)
    run_sweep('rho0.2-0.8_dp0.1-0.4', 0.2, 0.8, dp_min=0.10, dp_max=0.4)
    print('Generating regulatory graph ...')
    _plot_regulatory_graph()
    print('All sweeps complete.')


def _save_trajectory_tsv(traj: dict, world_seed: int, cell_seed: int, path: str):
    X_traj = traj['X_traj']
    P_traj = traj['P_traj']
    T_steps = X_traj.shape[0]
    G = X_traj.shape[1]

    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['world_seed', 'cell_seed', 'time', 'gene_id', 'X', 'P'])
        for t in range(T_steps):
            for g in range(G):
                writer.writerow([
                    world_seed, cell_seed, t, g,
                    f'{X_traj[t, g].item():.8f}',
                    f'{P_traj[t, g].item():.8f}',
                ])


def _save_metrics_tsv(all_metrics: list[dict], metrics_dir: str):
    path = os.path.join(metrics_dir, 'stability_summary.tsv')
    fields = [
        'world_seed', 'cell_seed',
        'converged', 'convergence_time',
        'bounded', 'collapsed',
        'max_expression', 'final_mean_expression',
        'activation_count', 'repression_count', 'repression_ratio',
    ]
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(fields)
        for m in all_metrics:
            writer.writerow([m[k] for k in fields])


def _load_trajectory_tsv(ws: int, cs: int, traj_dir: str) -> np.ndarray:
    path = os.path.join(traj_dir, f'trajectory_{ws}_{cs}.tsv')
    rows_by_time: dict[int, dict[int, float]] = {}
    with open(path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        for row in reader:
            t = int(row[2])
            g = int(row[3])
            p = float(row[5])
            if t not in rows_by_time:
                rows_by_time[t] = {}
            rows_by_time[t][g] = p
    T_steps = max(rows_by_time.keys()) + 1
    P = np.zeros((T_steps, dl.G))
    for t in range(T_steps):
        for g in range(dl.G):
            P[t, g] = rows_by_time[t].get(g, 0.0)
    return P


TF_COLORS: list[str] = ['#E74C3C', '#2ECC71', '#3498DB', '#F39C12']


def _plot_all_trajectories(worlds: dict[int, dl.World], world_seeds: list[int],
                          all_metrics: list[dict], plots_dir: str, traj_dir: str):
    conv_by_cell = {(m['world_seed'], m['cell_seed']): m['converged'] for m in all_metrics}
    WORLDS_PER_PAGE = 4
    pages = (len(world_seeds) + WORLDS_PER_PAGE - 1) // WORLDS_PER_PAGE

    for page in range(pages):
        page_seeds = world_seeds[page * WORLDS_PER_PAGE:(page + 1) * WORLDS_PER_PAGE]
        n_rows = len(page_seeds)
        n_cols = M_CELLS_PER_WORLD

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 3.5 * n_rows))
        if n_rows == 1:
            axes = axes.reshape(1, -1)

        for row, ws in enumerate(page_seeds):
            w = worlds[ws]
            for col in range(n_cols):
                cell_idx = (ws - BASE_SEED) * M_CELLS_PER_WORLD + col
                cell_seed = BASE_SEED + N_WORLDS + cell_idx

                P = _load_trajectory_tsv(ws, cell_seed, traj_dir)
                t_axis = np.arange(P.shape[0])

                ax = axes[row, col]

                for g in range(dl.N_TF):
                    ax.plot(t_axis, P[:, g],
                            color=TF_COLORS[g], linewidth=1.4,
                            label=f'TF{g}')

                for g in range(dl.N_TF, dl.G):
                    ax.plot(t_axis, P[:, g],
                            color='gray', alpha=0.4, linewidth=0.6)

                ax.legend(loc='upper right', fontsize=9, ncol=2,
                          framealpha=0.6, borderpad=0.3, labelspacing=0.2,
                          handlelength=1.0, handletextpad=0.3)
                status = 'converged' if conv_by_cell.get((ws, cell_seed), False) else 'unconverged'
                ax.set_title(f'W{ws} C{col+1} ({status})', fontsize=12)
                ax.grid(True, alpha=0.15)
                ax.tick_params(labelsize=6)
                if col == 0:
                    ax.set_ylabel('P', fontsize=8)
                if row == n_rows - 1:
                    ax.set_xlabel('t', fontsize=8)

        fig.suptitle(
            f'Protein Trajectories — Page {page+1}/{pages}',
            fontsize=16, y=1.01)
        plt.tight_layout()
        path = os.path.join(plots_dir, f'P_traj_page{page+1}.png')
        plt.savefig(path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f'  saved {path}')


def _plot_stability_distribution(per_world_metrics: list[dict], plots_dir: str):
    conv_times = [m['convergence_time'] for m in per_world_metrics if m['converged']]
    not_conv = sum(1 for m in per_world_metrics if not m['converged'])

    fig, ax = plt.subplots(figsize=(8, 5))
    if conv_times:
        n_bins = min(10, max(1, len(conv_times)))
        ax.hist(conv_times, bins=n_bins,
                edgecolor='black', alpha=0.7)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_xlabel('Convergence Time', fontsize=12)
    ax.set_ylabel('Number of Worlds', fontsize=12)
    title = f'Convergence Time Distribution ({len(conv_times)} converged'
    if not_conv > 0:
        title += f', {not_conv} not converged'
    title += ')'
    ax.set_title(title, fontsize=13)
    ax.grid(True, alpha=0.2)
    plt.tight_layout()
    path = os.path.join(plots_dir, 'stability_distribution.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  saved {path}')


def _plot_parameter_regime_scan(per_world_metrics: list[dict], plots_dir: str):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    rep_ratios = sorted(set(m['repression_ratio'] for m in per_world_metrics))
    labels = [f'{rr:.2f}' for rr in rep_ratios]
    x = np.arange(len(labels))

    conv_counts = [sum(1 for m in per_world_metrics
                      if m['repression_ratio'] == rr and m['converged'])
                   for rr in rep_ratios]
    unconv_counts = [sum(1 for m in per_world_metrics
                         if m['repression_ratio'] == rr and not m['converged'])
                      for rr in rep_ratios]

    width = 0.15
    ax1.bar(x - width/2, conv_counts, width, label='converged',
            color='#27AE60', edgecolor='black', linewidth=0.5)
    ax1.bar(x + width/2, unconv_counts, width, label='unconverged',
            color='#E74C3C', edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Repression Ratio', fontsize=11)
    ax1.set_ylabel('Number of Worlds', fontsize=11)
    ax1.set_title('Repression Ratio → Convergence', fontsize=12)
    ax1.set_xticks(x)
    ax1.set_xticklabels(labels)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.2, axis='y')
    ax1.yaxis.set_major_locator(plt.MaxNLocator(integer=True))

    conv_expr = [m['max_expression'] for m in per_world_metrics if m['converged']]
    unconv_expr = [m['max_expression'] for m in per_world_metrics if not m['converged']]
    bp = ax2.boxplot([conv_expr, unconv_expr], positions=[0, 1],
                     widths=0.3, patch_artist=True,
                     medianprops={'color': 'black', 'linewidth': 1.2})
    bp['boxes'][0].set_facecolor('#27AE60')
    bp['boxes'][0].set_alpha(0.5)
    bp['boxes'][1].set_facecolor('#E74C3C')
    bp['boxes'][1].set_alpha(0.5)
    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(['converged\n(n={})'.format(len(conv_expr)),
                         'unconverged\n(n={})'.format(len(unconv_expr))],
                        fontsize=10)
    ax2.set_ylabel('Max Expression', fontsize=11)
    ax2.set_title('Convergence → Max Expression', fontsize=12)
    ax2.grid(True, alpha=0.2, axis='y')

    fig.suptitle('Parameter Regime Scan', fontsize=13, y=1.01)
    plt.tight_layout()
    path = os.path.join(plots_dir, 'parameter_regime_scan.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  saved {path}')


def _plot_parameter_convergence_analysis(
    worlds: dict[int, dl.World],
    per_world_metrics: list[dict],
    plots_dir: str,
):
    params = ['rho', 'K', 'delta_x', 'delta_p']
    groups = {'converged': {k: [] for k in params},
              'non_converged': {k: [] for k in params}}
    drive_key = 'γ·ρ/(δx·δp)'
    groups['converged'][drive_key] = []
    groups['non_converged'][drive_key] = []

    for m in per_world_metrics:
        ws = m['world_seed']
        w = worlds[ws]
        label = 'converged' if m['converged'] else 'non_converged'
        rho_m = float(w.rho.mean().item())
        dx_m = float(w.delta_x.mean().item())
        dp_m = float(w.delta_p.mean().item())
        gamma_m = float(w.gamma.mean().item())
        for key in params:
            val = float(getattr(w, key).mean().item())
            groups[label][key].append(val)
        groups[label][drive_key].append(gamma_m * rho_m / (dx_m * dp_m))

    fig, axes = plt.subplots(3, 2, figsize=(13, 14))
    param_labels = {
        'rho': 'rho (transcription rate)',
        'K': 'K (Hill constant)',
        'delta_x': 'delta_x (mRNA decay)',
        'delta_p': 'delta_p (protein decay)',
    }

    for idx, (key, label) in enumerate(param_labels.items()):
        ax = axes[idx // 2][idx % 2]
        positions = [1, 2]
        data_conv = groups['converged'][key]
        data_non = groups['non_converged'][key]

        bp = ax.boxplot([data_conv, data_non], positions=positions,
                        widths=0.3, patch_artist=True,
                        medianprops={'color': 'black', 'linewidth': 1.2})
        bp['boxes'][0].set_facecolor('#27AE60')
        bp['boxes'][0].set_alpha(0.35)
        bp['boxes'][1].set_facecolor('#E74C3C')
        bp['boxes'][1].set_alpha(0.35)

        jitter_conv = np.random.default_rng(42).uniform(-0.08, 0.08, len(data_conv))
        jitter_non = np.random.default_rng(42).uniform(-0.08, 0.08, len(data_non))
        ax.scatter([1 + j for j in jitter_conv], data_conv,
                   c='#27AE60', s=50, edgecolors='black', linewidth=0.5, zorder=3)
        ax.scatter([2 + j for j in jitter_non], data_non,
                   c='#E74C3C', s=50, edgecolors='black', linewidth=0.5, zorder=3)

        ax.set_xticks([1, 2])
        ax.set_xticklabels(['converged\n(n={})'.format(len(data_conv)),
                            'non-converged\n(n={})'.format(len(data_non))],
                           fontsize=10)
        ax.set_ylabel(label, fontsize=11)
        ax.set_title(f'{key}', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.2, axis='y')

    ax_drive = axes[2][0]
    d_conv = groups['converged'][drive_key]
    d_non = groups['non_converged'][drive_key]
    bp = ax_drive.boxplot([d_conv, d_non], positions=[1, 2],
                          widths=0.3, patch_artist=True,
                          medianprops={'color': 'black', 'linewidth': 1.2})
    bp['boxes'][0].set_facecolor('#27AE60')
    bp['boxes'][0].set_alpha(0.35)
    bp['boxes'][1].set_facecolor('#E74C3C')
    bp['boxes'][1].set_alpha(0.35)
    j_conv = np.random.default_rng(42).uniform(-0.08, 0.08, len(d_conv))
    j_non = np.random.default_rng(42).uniform(-0.08, 0.08, len(d_non))
    ax_drive.scatter([1 + j for j in j_conv], d_conv,
                     c='#27AE60', s=50, edgecolors='black', linewidth=0.5, zorder=3)
    ax_drive.scatter([2 + j for j in j_non], d_non,
                     c='#E74C3C', s=50, edgecolors='black', linewidth=0.5, zorder=3)
    ax_drive.set_xticks([1, 2])
    ax_drive.set_xticklabels(['converged\n(n={})'.format(len(d_conv)),
                              'non-converged\n(n={})'.format(len(d_non))],
                             fontsize=10)
    ax_drive.set_ylabel('γ·ρ / (δx·δp)', fontsize=11)
    ax_drive.set_title('effective drive', fontsize=13, fontweight='bold')
    ax_drive.grid(True, alpha=0.2, axis='y')

    axes[2][1].set_visible(False)

    fig.suptitle('Parameter Profiles: Converged vs Non-Converged Worlds',
                 fontsize=14, y=1.01)
    plt.tight_layout()
    path = os.path.join(plots_dir, 'parameter_convergence_analysis.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f'  saved {path}')


def _plot_regulatory_graph():
    default_world_file = os.path.join(RUN_DIR, 'data', 'default', 'worlds',
                                      f'world_{BASE_SEED}.json')
    if not os.path.exists(default_world_file):
        return
    with open(default_world_file, 'r') as f:
        raw = json.loads(f.read())
    graph = {int(k): v for k, v in raw['P_graph'].items()}
    signs = {int(k): {int(sk): sv for sk, sv in v.items()}
             for k, v in raw['edge_signs'].items()}

    plots_dir = os.path.join(RUN_DIR, 'plots')
    os.makedirs(plots_dir, exist_ok=True)

    n_tf = dl.N_TF
    n_total = dl.G

    tf_node_ids = list(range(n_tf))
    target_node_ids = list(range(n_tf, n_total))

    angles = np.linspace(0, 2 * np.pi, n_total, endpoint=False)
    tf_angles = np.linspace(0, 2 * np.pi, n_tf, endpoint=False)
    tf_r = 1.5
    target_r = 3.5

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_aspect('equal')
    ax.axis('off')

    tf_pos = {}
    for i, g in enumerate(tf_node_ids):
        angle = tf_angles[i] - np.pi / 2
        x = tf_r * np.cos(angle)
        y = tf_r * np.sin(angle)
        tf_pos[g] = (x, y)
        ax.scatter(x, y, s=500, c=TF_COLORS[i], edgecolors='black',
                   linewidth=1.5, zorder=5)
        # 标签沿径向朝外放置，避免与TF间箭头重叠
        label_r = tf_r + 0.55
        lx = label_r * np.cos(angle)
        ly = label_r * np.sin(angle)
        ax.text(lx, ly, f'TF{g}', ha='center', va='center',
                fontsize=14, fontweight='bold', color=TF_COLORS[i])

    for g in target_node_ids:
        i = target_node_ids.index(g)
        angle = angles[n_tf + i] - np.pi / 2
        x = target_r * np.cos(angle)
        y = target_r * np.sin(angle)
        ax.scatter(x, y, s=80, c='#CCCCCC', edgecolors='#888888',
                   linewidth=0.8, zorder=4)
        ax.text(x, y - 0.32, f'{g}', ha='center', va='top',
                fontsize=12, color='#555555')

    from matplotlib.patches import FancyArrowPatch
    for g in range(n_total):
        reg = graph[g][0]
        sign = signs[g][reg]
        if g in tf_node_ids:
            x1, y1 = tf_pos[reg]
            x2, y2 = tf_pos[g]
        elif reg in tf_node_ids:
            x1, y1 = tf_pos[reg]
            angle = angles[g] - np.pi / 2
            x2 = target_r * np.cos(angle)
            y2 = target_r * np.sin(angle)
        else:
            continue

        dx = x2 - x1
        dy = y2 - y1
        d = np.sqrt(dx**2 + dy**2)
        r1 = 0.38 if reg in tf_node_ids else 0.12
        r2 = 0.38 if g in tf_node_ids else 0.10
        sx = x1 + dx * r1 / d
        sy = y1 + dy * r1 / d
        ex = x2 - dx * r2 / d
        ey = y2 - dy * r2 / d

        color = '#E74C3C' if sign == dl.ACTIVATION else '#3498DB'
        lw = 1.2 if (g in tf_node_ids and reg in tf_node_ids) else 0.4
        alpha = 0.9 if (g in tf_node_ids and reg in tf_node_ids) else 0.3
        if sign == dl.ACTIVATION:
            arrow = FancyArrowPatch((sx, sy), (ex, ey),
                                    arrowstyle='-|>', mutation_scale=28,
                                    color=color, linewidth=lw * 1.5, alpha=alpha,
                                    zorder=2,
                                    connectionstyle='arc3,rad=0.05')
            ax.add_patch(arrow)
        else:
            ax.plot([sx, ex], [sy, ey], color=color,
                    linewidth=lw * 1.5, alpha=alpha, zorder=2)
            dx = ex - sx
            dy = ey - sy
            d = np.sqrt(dx**2 + dy**2)
            if d > 0:
                px = -dy / d
                py = dx / d
                bar_len = 0.16 if (g in tf_node_ids and reg in tf_node_ids) else 0.09
                ax.plot([ex - bar_len * px, ex + bar_len * px],
                        [ey - bar_len * py, ey + bar_len * py],
                        color=color, linewidth=lw * 1.5, alpha=alpha, zorder=2)

    ax.set_xlim(-target_r - 0.8, target_r + 0.8)
    ax.set_ylim(-target_r - 0.8, target_r + 0.8)

    plt.tight_layout()
    path = os.path.join(plots_dir, 'regulatory_graph.png')
    plt.savefig(path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f'  saved {path}')


if __name__ == '__main__':
    main()
