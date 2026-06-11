"""
Run 002 — Model A Dynamical Characterization
==============================================

Per: docs/01_DDC_Lite_Curriculum/run_002/

Steps:
  1. parameter sweep config (4 strength_regimes x 2 sign_ratios)
  2. sample worlds + initial states + simulate trajectories
  3. per-trajectory analysis (equilibrium / stability / oscillation / classify)
  4. cross-world aggregation (parameter-regime mapping)
  5. figures / tables / summary

Author: zhanghl
"""
import json
import os
import sys
import copy
import csv
import torch
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Any

ATTRACTOR_NAMES = {
    'Type A': 'Stable equilibrium',
    'Type B': 'Slow convergence',
    'Type C': 'Damped oscillation',
    'Type D': 'Sustained oscillation',
    'Type E': 'Runaway divergence',
    'Type F': 'Numerical collapse',
    'Type G': 'Others',
}

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
import ddc_model_a as dma


RUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES_DIR = os.path.join(RUN_DIR, 'figures')
TABLES_DIR = os.path.join(RUN_DIR, 'tables')
TRAJECTORIES_DIR = os.path.join(RUN_DIR, 'trajectories')
METADATA_DIR = os.path.join(RUN_DIR, 'world_metadata')
PERTURBATIONS_DIR = os.path.join(RUN_DIR, 'perturbations')
ANALYSIS_DIR = os.path.join(RUN_DIR, 'analysis')
SUMMARY_DIR = os.path.join(RUN_DIR, 'summary')

N_TOPOS_PER_SIGN: int = 10
N_INIT_PER_WORLD: int = 5
T_SIM: int = dma.T_SIM

EPSILON: float = 1e-4
CONVERGENCE_WINDOW: int = 50
COLLAPSE_THRESHOLD: float = 1e-3
DIVERGENCE_THRESHOLD: float = 1e3
CLIPPING_FRAC_THRESHOLD: float = 0.1
SLOW_CONVERGENCE_THRESHOLD: int = 200

PERTURBATION_TIME: int = 500
PERTURBATION_SELECTION = [
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

STRENGTH_REGIMES = {
    'baseline':        (dma.A_MIN_BASELINE, dma.A_MAX_BASELINE, 'conservative baseline'),
    'stress':          (dma.A_MIN_STRESS, dma.A_MAX_STRESS, 'stress-test'),
    'chen_moderate':   (dma.A_MIN_CHEN_MODERATE, dma.A_MAX_CHEN_MODERATE, 'Chen-moderate'),
    'chen_stress':     (dma.A_MIN_CHEN_STRESS, dma.A_MAX_CHEN_STRESS, 'Chen-stress'),
}

SIGN_RATIOS = {
    'balanced':         0.5,
    'repression_biased': 0.333,
}

SIGN_RATIO_LABELS = {
    'balanced': 'Balanced',
    'repression_biased': 'Repression-biased',
}

SIGN_RATIO_RATIOS = {
    'balanced': '1:1',
    'repression_biased': '1:2',
}

BASE_SEED: int = 2000

SIGN_BASE_SEEDS = {
    'balanced': 2000,
    'repression_biased': 6000,
}


def make_world_id(strength_key: str, sign_key: str, topo_idx: int) -> str:
    return f"{strength_key}_{sign_key}_t{topo_idx:03d}"


def make_combo_key(strength_key: str, sign_key: str) -> str:
    return f"{strength_key}_{sign_key}"


def detect_equilibrium(X_traj: torch.Tensor) -> Dict[str, Any]:
    t_steps = X_traj.shape[0] - 1
    converged = False
    conv_time = -1

    for t in range(t_steps - CONVERGENCE_WINDOW + 1):
        all_consecutive = True
        for w in range(CONVERGENCE_WINDOW):
            diff = float(torch.norm(X_traj[t + w + 1] - X_traj[t + w]).item())
            if diff >= EPSILON:
                all_consecutive = False
                break
        if all_consecutive:
            converged = True
            conv_time = t
            break

    X_eq = X_traj[-1]
    eq_magnitude = float(torch.norm(X_eq).item())
    eq_sparsity = float((X_eq < EPSILON).sum().item()) / dma.G

    return {
        'converged': converged,
        'convergence_time': conv_time,
        'equilibrium_magnitude': eq_magnitude,
        'equilibrium_sparsity': eq_sparsity,
    }


def compute_spectral_info(world_dict: Dict[str, Any]) -> Tuple[float, Dict[str, List[float]]]:
    A = np.zeros((dma.G, dma.G))
    delta = world_dict['parameters']['delta']
    for i in range(dma.G):
        A[i, i] = 1.0 - delta[i]
    P_graph = {int(k): v for k, v in world_dict['P_graph'].items()}
    edge_signs = {int(k): {int(sk): sv for sk, sv in v.items()}
                  for k, v in world_dict['edge_signs'].items()}
    edge_strengths_dict = world_dict.get('edge_strengths', {})
    edge_strengths = {}
    if edge_strengths_dict:
        edge_strengths = {int(k): {int(sk): sv for sk, sv in v.items()}
                          for k, v in edge_strengths_dict.items()}

    for i in range(dma.G):
        for j in P_graph.get(i, []):
            sign = edge_signs[i][j]
            strength = edge_strengths.get(i, {}).get(j, 0.0)
            A[i, j] = sign * strength

    eigvals = np.linalg.eigvals(A)
    abs_vals = np.abs(eigvals)
    spectral_radius = float(np.max(abs_vals))
    eig_dict = {
        'real_parts': [float(ev.real) for ev in eigvals],
        'imag_parts': [float(ev.imag) for ev in eigvals],
        'abs_values': [float(abs(ev)) for ev in eigvals],
        'spectral_radius': spectral_radius,
    }
    return spectral_radius, eig_dict


def analyze_stability(X_traj: torch.Tensor, clip_count: torch.Tensor,
                      total_clips: int) -> Dict[str, Any]:
    max_expression = float(X_traj.max().item())
    bounded = max_expression < DIVERGENCE_THRESHOLD
    final_mean = float(X_traj[-1].mean().item())
    collapsed = final_mean < COLLAPSE_THRESHOLD
    total_steps = X_traj.shape[0]
    clipping_dominated = total_clips > total_steps * dma.G * CLIPPING_FRAC_THRESHOLD

    divergence_time = -1
    if not bounded:
        for t in range(X_traj.shape[0]):
            if float(X_traj[t].max().item()) >= DIVERGENCE_THRESHOLD:
                divergence_time = t
                break

    return {
        'bounded': bounded,
        'max_expression': max_expression,
        'final_mean_expression': final_mean,
        'collapsed': collapsed,
        'divergence_existence': not bounded,
        'divergence_time': divergence_time if not bounded else None,
        'numerical_collapse': collapsed,
        'clipping_dominated': clipping_dominated,
        'total_clips': total_clips,
    }


def analyze_oscillation(X_traj: torch.Tensor,
                         converged: bool, conv_time: int) -> Dict[str, Any]:
    BURN_IN = 200
    T_total = X_traj.shape[0]

    if converged and conv_time > BURN_IN + 100:
        X = X_traj[BURN_IN:conv_time].numpy()
    else:
        X = X_traj[BURN_IN:].numpy()

    T, G_dim = X.shape
    if T < 50:
        return {
            'oscillation_exists': False,
            'oscillation_type': 'none',
            'amplitude': 0.0,
            'frequency': 0.0,
            'damping_rate': None,
            'oscillatory_genes': [],
        }

    oscillatory_genes = []
    amplitudes = []
    frequencies = []
    damping_rates = []

    for g in range(G_dim):
        x = X[:, g]
        signal_range = float(x.max() - x.min())
        if signal_range < EPSILON:
            continue

        extrema = np.zeros(T, dtype=np.int8)
        for t in range(1, T - 1):
            if x[t] > x[t - 1] and x[t] > x[t + 1]:
                extrema[t] = 1
            elif x[t] < x[t - 1] and x[t] < x[t + 1]:
                extrema[t] = -1

        ext_idx = np.where(extrema != 0)[0]
        if len(ext_idx) < 3:
            continue

        peak_pairs = []
        for i in range(len(ext_idx) - 1):
            s, e = ext_idx[i], ext_idx[i + 1]
            delta = abs(float(x[e] - x[s]))
            peak_pairs.append(delta)

        if len(peak_pairs) < 2:
            continue

        gene_amplitude = float(np.median(peak_pairs))
        relative_amplitude = gene_amplitude / signal_range if signal_range > EPSILON else 0.0

        MIN_REL_AMP = 0.01
        if relative_amplitude < MIN_REL_AMP:
            continue

        gene_freq = len(ext_idx) / (2.0 * T)

        if len(peak_pairs) >= 4:
            mid = len(peak_pairs) // 2
            early = np.mean(peak_pairs[:mid])
            late = np.mean(peak_pairs[mid:])
            damping = float((early - late) / early) if early > EPSILON else 0.0
        else:
            damping = 0.0

        oscillatory_genes.append(g)
        amplitudes.append(gene_amplitude)
        frequencies.append(gene_freq)
        damping_rates.append(damping)

    if not oscillatory_genes:
        return {
            'oscillation_exists': False,
            'oscillation_type': 'none',
            'amplitude': 0.0,
            'frequency': 0.0,
            'damping_rate': None,
            'oscillatory_genes': [],
        }

    avg_damping = float(np.median(damping_rates))
    osc_type = 'damped' if avg_damping > 0.05 else 'sustained'

    return {
        'oscillation_exists': True,
        'oscillation_type': osc_type,
        'amplitude': float(np.mean(amplitudes)),
        'frequency': float(np.mean(frequencies)),
        'damping_rate': avg_damping,
        'oscillatory_genes': oscillatory_genes,
    }


def classify_attractor(eq: Dict, st: Dict, osc: Dict) -> str:
    if st['divergence_existence']:
        return 'Type E'
    if st['numerical_collapse']:
        return 'Type F'
    if osc['oscillation_type'] == 'sustained':
        return 'Type D'
    if osc['oscillation_type'] == 'damped':
        return 'Type C'
    if eq['converged']:
        conv_time = eq['convergence_time']
        if conv_time >= 0 and conv_time <= SLOW_CONVERGENCE_THRESHOLD:
            return 'Type A'
        return 'Type B'
    return 'Type G'


def analyze_trajectory(X_traj: torch.Tensor, clip_count: torch.Tensor,
                       total_clips: int, world_dict: Dict[str, Any]) -> Dict[str, Any]:
    eq = detect_equilibrium(X_traj)
    st = analyze_stability(X_traj, clip_count, total_clips)
    osc = analyze_oscillation(X_traj, eq['converged'], eq['convergence_time'])
    attractor_type = classify_attractor(eq, st, osc)

    spectral_r, eig_dict = compute_spectral_info(world_dict)

    return {
        'equilibrium': eq,
        'stability': st,
        'oscillation': osc,
        'attractor_type': attractor_type,
        'spectral_radius': spectral_r,
        'eigenvalues': eig_dict,
    }


def simulate_topo_stratified(sign_key: str, sign_ratio: float,
                              save_trajectories: bool = True) -> List[Dict[str, Any]]:
    results = []
    base_seed = SIGN_BASE_SEEDS[sign_key]

    topo_pool = {}
    for ti in range(N_TOPOS_PER_SIGN):
        topo_seed = base_seed + ti * 100
        topo_pool[ti] = dma.sample_topology(topo_seed, sign_ratio=sign_ratio)

    for strength_key, (a_min, a_max, regime_label) in STRENGTH_REGIMES.items():
        combo_key = make_combo_key(strength_key, sign_key)
        print(f"  combo: {combo_key}  "
              f"(strength=[{a_min},{a_max}], sign_ratio={sign_ratio})")

        for ti in range(N_TOPOS_PER_SIGN):
            topo_seed = base_seed + ti * 100
            world_seed = topo_seed + 50 + list(STRENGTH_REGIMES.keys()).index(strength_key) * 10
            world_id = make_world_id(strength_key, sign_key, ti)

            w_topo = topo_pool[ti]
            world = dma.sample_world(world_seed, world=w_topo,
                                     a_min=a_min, a_max=a_max)

            world_dict = world.to_dict()
            world_dict['world_id'] = world_id
            world_dict['topo_idx'] = ti
            world_dict['strength_regime'] = strength_key
            world_dict['sign_ratio_label'] = sign_key
            world_dict['a_min'] = a_min
            world_dict['a_max'] = a_max
            world_dict['sign_ratio'] = sign_ratio
            world_dict['regime_label'] = regime_label
            world_dict['edge_density_r'] = dma.EDGE_DENSITY_R
            world_dict['t_sim'] = T_SIM
            world_dict['runtime_version'] = 'v0.2'

            with open(os.path.join(METADATA_DIR, f'{world_id}.json'), 'w') as f:
                json.dump(world_dict, f, indent=2)

            per_world_analyses = []
            for ci in range(N_INIT_PER_WORLD):
                cell_seed = world_seed + 1 + ci
                X0 = dma.sample_initial_state(cell_seed)
                traj = dma.simulate_single_cell(world, X0, t_steps=T_SIM)

                X_traj = traj['X_traj']
                clip_count = traj['clip_count']
                total_clips = traj['total_clips']

                analysis = analyze_trajectory(X_traj, clip_count, total_clips, world_dict)
                per_world_analyses.append(analysis)

                if save_trajectories:
                    traj_save = {
                        'X_traj': X_traj,
                        'clip_count': clip_count,
                        'total_clips': total_clips,
                        'world_id': world_id,
                        'world_seed': world_seed,
                        'cell_seed': cell_seed,
                        'attractor_type': analysis['attractor_type'],
                    }
                    traj_path = os.path.join(
                        TRAJECTORIES_DIR, f'{world_id}_cell{ci:02d}.pt'
                    )
                    torch.save(traj_save, traj_path)

            analysis_summary = {
                'world_id': world_id,
                'world_seed': world_seed,
                'combo_key': combo_key,
                'topo_idx': ti,
                'strength_key': strength_key,
                'sign_key': sign_key,
                'a_min': a_min,
                'a_max': a_max,
                'sign_ratio': sign_ratio,
                'per_cell': per_world_analyses,
            }

            cell_attractors = [a['attractor_type'] for a in per_world_analyses]
            attractor_counts = Counter(cell_attractors)
            primary_attractor = max(attractor_counts, key=attractor_counts.get)
            analysis_summary['primary_attractor'] = primary_attractor
            analysis_summary['attractor_distribution'] = dict(attractor_counts)
            analysis_summary['spectral_radius'] = per_world_analyses[0]['spectral_radius']
            analysis_summary['eigenvalues'] = per_world_analyses[0]['eigenvalues']
            world_converged = all(
                a['equilibrium']['converged'] for a in per_world_analyses
            )
            analysis_summary['converged'] = world_converged
            analysis_summary['mean_convergence_time'] = float(
                np.mean([a['equilibrium']['convergence_time']
                         for a in per_world_analyses])
            ) if world_converged else None
            analysis_summary['clipping_dominated'] = any(
                a['stability']['clipping_dominated'] for a in per_world_analyses
            )

            results.append(analysis_summary)

        n_worlds_this_combo = N_TOPOS_PER_SIGN
        combo_results = [r for r in results if r['combo_key'] == combo_key]
        types = [r['primary_attractor'] for r in combo_results]
        type_dist = Counter(types)
        print(f"    {n_worlds_this_combo} worlds: {dict(type_dist)}")

    return results


def aggregate_across_combos(all_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    combo_groups = defaultdict(list)
    for r in all_results:
        combo_groups[r['combo_key']].append(r)

    aggregation = {}
    for combo_key, worlds in combo_groups.items():
        n = len(worlds)
        attractor_counts = defaultdict(int)
        for w in worlds:
            attractor_counts[w['primary_attractor']] += 1

        aggregation[combo_key] = {
            'n_worlds': n,
            'attractor_distribution': dict(attractor_counts),
            'attractor_distribution_pct': {
                k: round(v / n * 100, 1) for k, v in attractor_counts.items()
            },
            'mean_spectral_radius': float(np.mean([w['spectral_radius'] for w in worlds])),
            'convergence_rate': float(sum(1 for w in worlds if w['converged']) / n),
            'bounded_rate': float(
                sum(1 for w in worlds
                    if all(a['stability']['bounded'] for a in w['per_cell'])
                ) / n
            ),
            'collapse_rate': float(
                sum(1 for w in worlds
                    if any(a['stability']['collapsed'] for a in w['per_cell'])
                ) / n
            ),
        }

    return aggregation


def build_world_summary_table(all_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    rows = []
    for r in all_results:
        cells = r['per_cell']
        n_cells = len(cells)
        conv_cells = [a for a in cells if a['equilibrium']['converged']]
        rows.append({
            'world_id': r['world_id'],
            'combo_key': r['combo_key'],
            'topo_idx': r.get('topo_idx', ''),
            'a_min': r['a_min'],
            'a_max': r['a_max'],
            'sign_ratio': r['sign_ratio'],
            'edge_density_r': dma.EDGE_DENSITY_R,
            'decay_delta': dma.DEFAULT_DELTA,
            'basal_transcription_b': dma.DEFAULT_B,
            'attractor_type': r['primary_attractor'],
            'attractor_distribution': r.get('attractor_distribution', {}),
            'spectral_radius': round(r['spectral_radius'], 6),
            'converged': r['converged'],
            'convergence_time': r.get('mean_convergence_time'),
            'eq_magnitude': round(np.mean([a['equilibrium']['equilibrium_magnitude'] for a in cells]), 4),
            'eq_sparsity': round(np.mean([a['equilibrium']['equilibrium_sparsity'] for a in cells]), 4),
            'bounded': all(a['stability']['bounded'] for a in cells),
            'max_expression': round(max(a['stability']['max_expression'] for a in cells), 4),
            'collapsed': any(a['stability']['collapsed'] for a in cells),
            'clipping_dominated': any(a['stability']['clipping_dominated'] for a in cells),
            'clipping_freq': round(np.mean([a['stability']['total_clips'] for a in cells]) / (T_SIM * dma.G), 4),
            'oscillation_exists': any(a['oscillation']['oscillation_exists'] for a in cells),
        })
    return rows


def aggregate_by_topo(all_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    topo_groups = defaultdict(list)
    for r in all_results:
        sign_key = r.get('sign_key', '')
        ti = r.get('topo_idx', '')
        key = f"{sign_key}_t{ti:03d}"
        topo_groups[key].append(r)

    topo_summaries = {}
    for topo_key, worlds in topo_groups.items():
        regime_map = {}
        for w in worlds:
            sk = w.get('strength_key', w['combo_key'].split('_')[0])
            regime_map[sk] = w['primary_attractor']

        topo_summaries[topo_key] = {
            'n_worlds': len(worlds),
            'strength_regimes': sorted(regime_map.keys()),
            'attractor_by_regime': regime_map,
            'n_distinct_attractors': len(set(regime_map.values())),
            'attractors': list(regime_map.values()),
        }

    return topo_summaries


def analyze_topo_strength_interaction(all_results: List[Dict[str, Any]]) -> Dict[str, Any]:
    topo_summaries = aggregate_by_topo(all_results)

    transition_matrix = defaultdict(lambda: defaultdict(int))
    n_transitions = 0
    n_same_attractor = 0

    for topo_key, summary in topo_summaries.items():
        regimes = summary['strength_regimes']
        attractors = summary['attractor_by_regime']
        if len(regimes) < 2:
            continue
        for i in range(len(regimes)):
            for j in range(i + 1, len(regimes)):
                a_i = attractors[regimes[i]]
                a_j = attractors[regimes[j]]
                if a_i == a_j:
                    n_same_attractor += 1
                else:
                    n_transitions += 1
                transition_matrix[a_i][a_j] += 1
                if a_i != a_j:
                    transition_matrix[a_j][a_i] += 1

    n_topo_diverse = sum(
        1 for s in topo_summaries.values()
        if s['n_distinct_attractors'] > 1
    )

    return {
        'n_topos': len(topo_summaries),
        'n_topo_diverse': n_topo_diverse,
        'topo_diverse_frac': round(n_topo_diverse / len(topo_summaries), 4) if topo_summaries else 0,
        'n_transitions': n_transitions,
        'n_same_attractor': n_same_attractor,
        'topo_summaries': {
            k: {sk: sv for sk, sv in v.items() if sk != 'attractors'}
            for k, v in topo_summaries.items()
        },
    }


def make_regime_label(combo_key: str) -> str:
    parts = combo_key.split('_', 1)
    return f"{parts[0]} | {parts[1]}"


def plot_canonical_trajectories(all_results: List[Dict[str, Any]]):
    exemplars = {
        'Type A': None, 'Type B': None, 'Type C': None,
        'Type D': None, 'Type E': None, 'Type F': None,
        'Type G': None,
    }
    for r in all_results:
        at = r['primary_attractor']
        if exemplars[at] is None:
            exemplars[at] = r['world_id']

    all_found = {k: v for k, v in exemplars.items() if v is not None}
    n_types = len(all_found)
    if n_types == 0:
        return

    fig, axes = plt.subplots(n_types, 1, figsize=(12, 3 * n_types), squeeze=False)
    for idx, (atype, wid) in enumerate(sorted(all_found.items())):
        ax = axes[idx][0]
        traj_files = sorted(
            [f for f in os.listdir(TRAJECTORIES_DIR) if f.startswith(wid) and f.endswith('.pt')]
        )
        if traj_files:
            data = torch.load(os.path.join(TRAJECTORIES_DIR, traj_files[0]), weights_only=False)
            X = data['X_traj'].numpy()
            for g in range(dma.G):
                ax.plot(X[:, g], alpha=0.9, linewidth=0.8)
            at_name = ATTRACTOR_NAMES.get(atype, atype)
            ax.set_title(f'{atype} — {at_name} — {wid}')
            ax.set_xlabel('t')
            ax.set_ylabel('X')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'canonical_trajectories.png'), dpi=150)
    plt.close(fig)


def plot_all_worlds_trajectories():
    os.makedirs(os.path.join(FIGURES_DIR, 'traj'), exist_ok=True)

    topos_per_fig = 4
    strength_order = list(STRENGTH_REGIMES.keys())
    rows = len(strength_order)

    for sign_key in SIGN_RATIOS:
        n_figs = (N_TOPOS_PER_SIGN + topos_per_fig - 1) // topos_per_fig
        for fig_idx in range(n_figs):
            t_start = fig_idx * topos_per_fig
            t_end = min(t_start + topos_per_fig, N_TOPOS_PER_SIGN)
            topo_batch = list(range(t_start, t_end))
            cols = len(topo_batch)

            fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 3),
                                     squeeze=False)

            for r_idx, topo_idx in enumerate(topo_batch):
                for s_idx, strength_key in enumerate(strength_order):
                    wid = make_world_id(strength_key, sign_key, topo_idx)
                    ax = axes[s_idx, r_idx]

                    traj_files = sorted(
                        [f for f in os.listdir(TRAJECTORIES_DIR)
                         if f.startswith(wid) and f.endswith('.pt')]
                    )
                    if not traj_files:
                        ax.set_title(f'{wid} — no data', fontsize=11)
                        continue

                    data = torch.load(os.path.join(TRAJECTORIES_DIR, traj_files[0]),
                                      weights_only=False)
                    X = data['X_traj'].numpy()
                    at = data['attractor_type']
                    at_name = ATTRACTOR_NAMES.get(at, at)
                    for g in range(dma.G):
                        ax.plot(X[:, g], alpha=0.95, linewidth=1.5)
                    ax.set_title(f'[{at}: {at_name}]', fontsize=10)
                    ax.tick_params(labelsize=7)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    if s_idx == 0:
                        ax.set_title(
                            r'$\bf{topology\ ' + f'{topo_idx:03d}' + r'}$' + f'\n[{at}: {at_name}]',
                            fontsize=10)
                    else:
                        ax.set_title(f'[{at}: {at_name}]', fontsize=10)

            for s_idx, strength_key in enumerate(strength_order):
                a_min, a_max, _ = STRENGTH_REGIMES[strength_key]
                axes[s_idx, 0].set_ylabel(
                    f'{strength_key}\na ∈ [{a_min}, {a_max}]',
                    fontsize=11, rotation=0, ha='right', va='center', labelpad=10)

            fig.suptitle(
                f'World Trajectories (Cell 0) — {SIGN_RATIO_LABELS[sign_key]} '
                f'(activation:repression={SIGN_RATIO_RATIOS[sign_key]})\n'
                f'topo {t_start:03d}-{t_end-1:03d}\n'
                f'(columns = topology, rows = strength regime)',
                fontsize=12, fontweight='bold')
            plt.tight_layout()
            fig.savefig(os.path.join(FIGURES_DIR, 'traj',
                        f'all_trajectories_{sign_key}_{fig_idx + 1:02d}.png'),
                        dpi=150)
            plt.close(fig)


def plot_spectral_radius(all_results: List[Dict[str, Any]]):
    attractor_order = [
        ('Type A', 'Stable equilibrium'),
        ('Type B', 'Slow convergence'),
        ('Type C', 'Damped oscillation'),
        ('Type D', 'Sustained oscillation'),
        ('Type E', 'Runaway divergence'),
        ('Type F', 'Numerical collapse'),
        ('Type G', 'Others'),
    ]

    by_type_raw = defaultdict(list)
    by_type_no_clip = defaultdict(list)
    by_type_clip_world_count = defaultdict(int)
    for r in all_results:
        at = r['primary_attractor']
        by_type_raw[at].append(r['spectral_radius'])
        if not all(a['stability']['clipping_dominated'] for a in r['per_cell']):
            by_type_no_clip[at].append(r['spectral_radius'])
        if any(a['stability']['clipping_dominated'] for a in r['per_cell']):
            by_type_clip_world_count[at] += 1

    labels, data, label_text = [], [], []
    for at, name in attractor_order:
        if by_type_raw[at]:
            n_total = len(by_type_raw[at])
            n_clip = by_type_clip_world_count[at]
            labels.append(at)
            data.append(by_type_raw[at])
            label_text.append(f"{at} - {name}\n(n={n_total}, clip={n_clip}/{n_total})")
        if at == 'Type A' and by_type_no_clip['Type A']:
            n_nc = len(by_type_no_clip['Type A'])
            labels.append('Type A*')
            data.append(by_type_no_clip['Type A'])
            label_text.append(f"Type A* - non-clipping\n(n={n_nc}, clip=0/{n_nc})")

    if not data:
        return

    fig, ax = plt.subplots(figsize=(10, 6))
    positions = list(range(1, len(labels) + 1))
    bp = ax.boxplot(data, positions=positions, tick_labels=label_text, patch_artist=True, widths=0.6)
    for patch in bp['boxes']:
        patch.set_facecolor('#a8c5e8')
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='ρ=1')
    ax.set_xticks(positions)
    ax.set_xticklabels(label_text, fontsize=7)
    ax.set_ylabel('Spectral Radius ρ(A)')
    ax.set_title('ρ(A) by Attractor Type\n(label: n = world count, clip = number of clipping-dominated worlds)\n(* = only worlds with at least one non-clipping cell)')
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'spectral_radius.png'), dpi=150)
    plt.close(fig)


def plot_convergence_distribution(all_results: List[Dict[str, Any]]):
    sign_colors = {'balanced': '#3498db', 'repression_biased': '#9b59b6'}

    by_sign = defaultdict(list)
    for r in all_results:
        sk = r.get('sign_key')
        if sk not in sign_colors:
            continue
        if r.get('converged') and r.get('mean_convergence_time') is not None:
            by_sign[sk].append(r['mean_convergence_time'])

    if not by_sign:
        return

    fig, ax = plt.subplots(figsize=(8, 5))

    for sk, color in sign_colors.items():
        times = sorted(by_sign.get(sk, []))
        if not times:
            continue
        cdf = np.arange(1, len(times) + 1) / len(times)
        ax.plot(times, cdf, color=color, linewidth=2, label=f'{sk} (n={len(times)})',
                marker='o', markersize=3, markevery=max(1, len(times) // 20))

    ax.set_xlabel('Convergence Time (steps)')
    ax.set_ylabel('Cumulative Proportion')
    ax.set_title('Convergence Time Distribution by Sign Ratio')
    ax.legend(fontsize=9)
    ax.set_ylim(0, 1.05)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))
    ax.set_xlim(left=0)

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'convergence_time.png'), dpi=150)
    plt.close(fig)


def plot_clipping_distribution(all_results: List[Dict[str, Any]]):
    attractor_order = [
        ('Type A', 'Stable equilibrium', '#2ecc71'),
        ('Type B', 'Slow convergence', '#3498db'),
        ('Type C', 'Damped oscillation', '#f39c12'),
        ('Type D', 'Sustained oscillation', '#e74c3c'),
        ('Type E', 'Runaway divergence', '#9b59b6'),
        ('Type G', 'Others', '#95a5a6'),
    ]
    sign_keys = list(SIGN_RATIOS.keys())
    sign_colors = {'balanced': '#3498db', 'repression_biased': '#9b59b6'}

    by_type = {at: defaultdict(list) for at, _, _ in attractor_order}
    for r in all_results:
        at = r['primary_attractor']
        sk = r.get('sign_key')
        if at not in by_type or sk not in sign_keys:
            continue
        for a in r['per_cell']:
            freq = a['stability']['total_clips'] / ((T_SIM + 1) * dma.G)
            by_type[at][sk].append(freq)

    fig, ax = plt.subplots(figsize=(10, 6))
    x = np.arange(len(attractor_order))
    width = 0.18
    has_data = False

    for si, sk in enumerate(sign_keys):
        positions = x + (si - (len(sign_keys) - 1) / 2) * width
        data = [by_type[at][sk] for at, _, _ in attractor_order]
        if any(len(d) > 0 for d in data):
            has_data = True
        bp = ax.boxplot(data, positions=positions, widths=width * 0.9,
                        patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor(sign_colors[sk])
            patch.set_alpha(0.7)

    if not has_data:
        plt.close(fig)
        return

    ax.set_xticks(x)
    ax.set_xticklabels([f'{at}\n{name}' for at, name, _ in attractor_order], fontsize=10)
    ax.set_ylabel('Clipping Frequency')
    ax.set_title('Clipping Frequency Distribution by Attractor Type\n(grouped by sign_ratio: balanced = blue, repression_biased = purple)')

    from matplotlib.patches import Patch
    legend_patches = [Patch(color=c, label=sk) for sk, c in sign_colors.items()]
    ax.legend(handles=legend_patches, fontsize=9)

    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'clipping_frequency.png'), dpi=150)
    plt.close(fig)


def write_world_summary_csv(rows: List[Dict[str, Any]]):
    if not rows:
        return
    path = os.path.join(TABLES_DIR, 'world_summary.tsv')
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_regime_summary_csv(aggregation: Dict[str, Any],
                              all_results: List[Dict[str, Any]]):
    combo_meta = {}
    for r in all_results:
        ck = r['combo_key']
        if ck not in combo_meta:
            combo_meta[ck] = (r['a_min'], r['a_max'])

    path = os.path.join(TABLES_DIR, 'regime_summary.tsv')
    attractor_types = ['Type A', 'Type B', 'Type C', 'Type D', 'Type E', 'Type F', 'Type G']
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        header = ['combo_key', 'n_worlds', 'a_min', 'a_max', 'decay_delta',
                   'convergence_rate', 'bounded_rate',
                   'collapse_rate', 'mean_spectral_radius'] + attractor_types
        writer.writerow(header)
        for combo_key, agg in sorted(aggregation.items()):
            a_min, a_max = combo_meta.get(combo_key, ('', ''))
            row = [
                combo_key, agg['n_worlds'],
                a_min, a_max, dma.DEFAULT_DELTA,
                round(agg['convergence_rate'], 4),
                round(agg['bounded_rate'], 4),
                round(agg['collapse_rate'], 4),
                round(agg['mean_spectral_radius'], 6),
            ]
            for at in attractor_types:
                row.append(agg['attractor_distribution_pct'].get(at, 0))
            writer.writerow(row)


def write_analysis_json(all_results: List[Dict[str, Any]],
                         aggregation: Dict[str, Any]):
    path = os.path.join(ANALYSIS_DIR, 'aggregated_analysis.json')
    with open(path, 'w') as f:
        json.dump(aggregation, f, indent=2, default=str)

    path = os.path.join(ANALYSIS_DIR, 'all_world_results.json')
    serializable = []
    for r in all_results:
        item = {k: v for k, v in r.items() if k != 'per_cell'}
        item = {k: (float(v) if isinstance(v, (np.floating,)) else v)
                for k, v in item.items()}
        serializable.append(item)
    with open(path, 'w') as f:
        json.dump(serializable, f, indent=2, default=str)


def plot_topo_strength_interaction(topo_analysis: Dict[str, Any]):
    summaries = topo_analysis.get('topo_summaries', {})
    if not summaries:
        return

    attractor_colors = {
        'Type A': '#2ecc71',
        'Type B': '#3498db',
        'Type C': '#f39c12',
        'Type D': '#e74c3c',
        'Type E': '#9b59b6',
        'Type G': '#95a5a6',
    }
    topo_keys = sorted(summaries.keys())
    n_topos = len(topo_keys)
    if n_topos == 0:
        return

    used_attractors = set()
    for tk in topo_keys:
        for sk in STRENGTH_REGIMES.keys():
            at = summaries[tk]['attractor_by_regime'].get(sk, '?')
            if at in attractor_colors:
                used_attractors.add(at)

    fig, ax = plt.subplots(figsize=(10, 6))

    strength_order = list(STRENGTH_REGIMES.keys())
    yticks = []
    ylabels = []
    cell_w, cell_h = 0.85, 0.75
    y_offset = 1.0

    sign_keys = list(SIGN_RATIOS.keys())
    topo_nums = sorted({tk.split("_t")[-1] for tk in topo_keys},
                       key=lambda s: int(s))
    available_sign_keys = [sk for sk in sign_keys
                           if any(f'{sk}_t{tn}' in summaries for tn in topo_nums)]

    intra_gap = 0.9
    inter_gap = 1.5
    y_positions = []
    ordered_keys = []
    for si, sign_key in enumerate(available_sign_keys):
        if si > 0:
            y_offset += inter_gap
        for tn in topo_nums:
            key = f'{sign_key}_t{tn}'
            if key in summaries:
                y_positions.append(y_offset)
                ordered_keys.append(key)
                y_offset += intra_gap

    for row_idx, (topo_key, y) in enumerate(zip(ordered_keys, y_positions)):
        yticks.append(y)
        topo_num = topo_key.split("_t")[-1]
        sign_key = topo_key.split("_t")[0]
        ylabels.append(f'{sign_key}_topo{topo_num}')

        regime_map = summaries[topo_key]['attractor_by_regime']
        for xi, sk in enumerate(strength_order):
            at = regime_map.get(sk, '?')
            color = attractor_colors.get(at, '#cccccc')
            ax.barh(y, cell_w, left=xi + (1 - cell_w) / 2,
                    height=cell_h, color=color,
                    edgecolor='white', linewidth=0.3)

    ax.invert_yaxis()
    ax.set_yticks(yticks)
    ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xticks([i + 0.5 for i in range(len(strength_order))])
    ax.set_xticklabels([f'{sk}\na ∈ [{STRENGTH_REGIMES[sk][0]}, {STRENGTH_REGIMES[sk][1]}]'
                        for sk in strength_order],
                       rotation=0, fontsize=8, ha='center')
    ax.set_xlim(0, len(strength_order))
    ax.set_ylim(20.5, -1.5)

    sign_desc = '  /  '.join(
        f'{SIGN_RATIO_LABELS[sk]} (activation:repression = {SIGN_RATIO_RATIOS[sk]})'
        for sk in SIGN_RATIOS
    )
    fig.suptitle('Topology × Strength Regime: Attractor Type',
                 fontsize=12, fontweight='bold', y=0.97)
    fig.text(0.5, 0.915, sign_desc, fontsize=9, ha='center')

    for xi in range(1, len(strength_order)):
        ax.axvline(x=xi, color='white', linewidth=0.5)

    from matplotlib.patches import Patch
    legend_patches = [Patch(color=attractor_colors[at], label=at)
                      for at in sorted(used_attractors)]
    ax.legend(handles=legend_patches, fontsize=8, ncol=len(used_attractors),
              loc='upper right', bbox_to_anchor=(0.99, 0.99), framealpha=0.9)

    fig.subplots_adjust(left=0.17, right=0.97, top=0.90, bottom=0.08)
    fig.savefig(os.path.join(FIGURES_DIR, 'topo_strength_interaction.png'), dpi=150)
    plt.close(fig)


def plot_repression_ratio_vs_stability(all_results: List[Dict[str, Any]]):
    attractor_labels = [
        ('Type A', 'Stable'),
        ('Type B', 'Slow'),
        ('Type C', 'Damped'),
        ('Type D', 'Sustained'),
        ('Type E', 'Divergence'),
        ('Type F', 'Collapse'),
        ('Type G', 'Others'),
    ]
    sign_keys = list(SIGN_RATIOS.keys())
    combos = list(STRENGTH_REGIMES.keys())

    counts = {sk: {c: {at: 0 for at, _ in attractor_labels} for c in combos} for sk in sign_keys}
    for r in all_results:
        sk = r.get('sign_key')
        ck = r.get('strength_key')
        at = r['primary_attractor']
        if sk in counts and ck in counts[sk] and at in counts[sk][ck]:
            counts[sk][ck][at] += 1

    fig, ax = plt.subplots(figsize=(10, 6))
    n_signs = len(sign_keys)
    width = 0.18
    x = np.arange(len(combos))
    colors_list = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#9b59b6', '#34495e', '#95a5a6']
    sign_short = {'balanced': 'bal', 'repression_biased': 'rep'}

    for si, sk in enumerate(sign_keys):
        bottoms = np.zeros(len(combos))
        for ai, (at, name) in enumerate(attractor_labels):
            vals = np.array([counts[sk][c][at] for c in combos])
            positions_x = x + (si - (n_signs - 1) / 2) * width * 1.4
            bars = ax.bar(positions_x, vals, width,
                   bottom=bottoms, color=colors_list[ai], label=f'{at} - {name}' if si == 0 else None)
            # 标注各段数量（在bar内部居中）
            for rect, b0, v in zip(bars, bottoms, vals):
                if v > 0:
                    ax.text(rect.get_x() + rect.get_width() / 2, b0 + v / 2, str(int(v)),
                            ha='center', va='center', fontsize=8, color='white', fontweight='bold')
            bottoms += vals

    # 为每个 combo 下的两个 sign_ratio 添加底部 label
    for ci, c in enumerate(combos):
        for si, sk in enumerate(sign_keys):
            pos = ci + (si - (n_signs - 1) / 2) * width * 1.4
            total = sum(counts[sk][c].values())
            ax.text(pos, -0.6, sign_short.get(sk, sk),
                    ha='center', va='top', fontsize=7, color='gray')

    ax.set_xticks(x)
    ax.set_xticklabels(combos)
    ax.set_xlabel('strength regime')
    ax.set_ylabel('world count')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    max_stack = max(sum(counts[sk][c].values()) for sk in sign_keys for c in combos)
    ax.set_ylim(bottom=-1.2, top=max_stack * 1.1 + 1)
    ax.set_title('Repression Ratio vs Stability\n(grouped by sign_ratio: bal=balanced, rep=repression_biased)')
    from matplotlib.patches import Patch
    legend_patches = [Patch(color=c, label=at) for at, c in zip([a for a, _ in attractor_labels], colors_list)]
    ax.legend(handles=legend_patches, fontsize=7, ncol=2, loc='upper left')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'repression_ratio_vs_stability.png'), dpi=150)
    plt.close(fig)


def plot_clipping_dominated_by_regime(all_results: List[Dict[str, Any]]):
    sign_keys = list(SIGN_RATIOS.keys())
    combos = list(STRENGTH_REGIMES.keys())

    counts = {sk: {c: {'clip': 0, 'total': 0} for c in combos} for sk in sign_keys}
    for r in all_results:
        sk = r.get('sign_key')
        ck = r.get('strength_key')
        if sk in counts and ck in counts[sk]:
            counts[sk][ck]['total'] += 1
            if r.get('clipping_dominated'):
                counts[sk][ck]['clip'] += 1

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(combos))
    width = 0.35
    sign_colors = {'balanced': '#3498db', 'repression_biased': '#9b59b6'}
    sign_short = {'balanced': 'balanced (1:1)', 'repression_biased': 'repression_biased (1:2)'}

    for si, sk in enumerate(sign_keys):
        fracs = [counts[sk][c]['clip'] / counts[sk][c]['total']
                 if counts[sk][c]['total'] > 0 else 0.0
                 for c in combos]
        positions_x = x + (si - 0.5) * width
        bars = ax.bar(positions_x, [f * 100 for f in fracs], width,
                      color=sign_colors[sk], label=sign_short[sk],
                      edgecolor='white', linewidth=0.5)
        for rect, frac, c in zip(bars, fracs, combos):
            clip_n = counts[sk][c]['clip']
            total_n = counts[sk][c]['total']
            ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + 1,
                    f'{frac*100:.0f}%\n({clip_n}/{total_n})',
                    ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(combos)
    ax.set_xlabel('strength regime')
    ax.set_ylabel('clipping-dominated worlds (%)')
    ax.set_ylim(0, 115)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y)}%'))
    ax.set_title('Clipping-Dominated Worlds by Regime and Sign Ratio')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'clipping_dominated_by_regime.png'), dpi=150)
    plt.close(fig)


def plot_in_degree_distribution_vs_stability(all_results: List[Dict[str, Any]]):
    attractor_labels = [
        ('Type A', 'Stable'),
        ('Type B', 'Slow'),
        ('Type C', 'Damped'),
        ('Type D', 'Sustained'),
        ('Type E', 'Divergence'),
        ('Type F', 'Collapse'),
        ('Type G', 'Others'),
    ]

    world_data = []
    for r in all_results:
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        if not os.path.exists(meta_path):
            continue
        with open(meta_path) as f:
            meta = json.load(f)
        in_deg = meta.get('in_degrees', {})
        if not in_deg:
            continue
        values = np.array([int(v) for v in in_deg.values()], dtype=float)
        std = float(np.std(values, ddof=0))
        world_data.append((std, r['primary_attractor']))

    if not world_data:
        return

    stds = np.array([w[0] for w in world_data])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])
    bucket_edges = [0.0, q25, q50, q75, np.inf]

    def assign_bucket(v):
        if v < q25:    return 'Q1 (low)'
        if v < q50:    return 'Q2'
        if v < q75:    return 'Q3'
        return 'Q4 (high)'

    degree_buckets = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)']
    counts = {at: {b: 0 for b in degree_buckets} for at, _ in attractor_labels}
    for std, at in world_data:
        b = assign_bucket(std)
        for a, _ in attractor_labels:
            if at == a:
                counts[a][b] += 1
                break

    fig, ax = plt.subplots(figsize=(10, 6))
    n_at = len(attractor_labels)
    width = 0.18
    x = np.arange(len(degree_buckets)) * 1.4
    colors_list = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#9b59b6', '#34495e', '#95a5a6']

    for ai, (at, name) in enumerate(attractor_labels):
        vals = [counts[at][b] for b in degree_buckets]
        positions_x = x + (ai - (n_at - 1) / 2) * width
        bars = ax.bar(positions_x, vals, width,
                color=colors_list[ai], label=f'{at} - {name}')
        for rect, v in zip(bars, vals):
            if v > 0:
                ax.text(rect.get_x() + rect.get_width() / 2, v + 0.3, str(int(v)),
                        ha='center', va='bottom', fontsize=10, color='black')

    ax.set_xticks(x)
    labels_xticks = [
        f"{b}\nstd range\n[{bucket_edges[i]:.2f}, {bucket_edges[i+1]:.2f})" if i < 3
        else f"{b}\nstd ≥ {bucket_edges[i]:.2f}"
        for i, b in enumerate(degree_buckets)
    ]
    ax.set_xticklabels(labels_xticks)
    ax.set_xlabel('std of in-degree distribution (quartile bins)')
    ax.set_ylabel('world count (integer)')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(top=max(max(c.values()) for c in counts.values()) * 1.15 + 1)
    ax.set_title(f'Topology Realization / In-Degree Std vs Stability\n(std = std of per-gene in-degree)')
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'in_degree_distribution_vs_stability.png'), dpi=150)
    plt.close(fig)


def plot_out_degree_distribution_vs_stability(all_results: List[Dict[str, Any]]):
    """Plot out-degree distribution std vs stability, same format as in-degree version."""
    attractor_labels = [
        ('Type A', 'Stable'),
        ('Type B', 'Slow'),
        ('Type C', 'Damped'),
        ('Type D', 'Sustained'),
        ('Type E', 'Divergence'),
        ('Type F', 'Collapse'),
        ('Type G', 'Others'),
    ]

    world_data = []
    for r in all_results:
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        if not os.path.exists(meta_path):
            continue
        with open(meta_path) as f:
            meta = json.load(f)
        out_deg = meta.get('out_degrees', {})
        if not out_deg:
            continue
        values = np.array([int(v) for v in out_deg.values()], dtype=float)
        std = float(np.std(values, ddof=0))
        world_data.append((std, r['primary_attractor']))

    if not world_data:
        return

    stds = np.array([w[0] for w in world_data])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])
    bucket_edges = [0.0, q25, q50, q75, np.inf]

    def assign_bucket(v):
        if v < q25:    return 'Q1 (low)'
        if v < q50:    return 'Q2'
        if v < q75:    return 'Q3'
        return 'Q4 (high)'

    degree_buckets = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)']
    counts = {at: {b: 0 for b in degree_buckets} for at, _ in attractor_labels}
    for std, at in world_data:
        b = assign_bucket(std)
        for a, _ in attractor_labels:
            if at == a:
                counts[a][b] += 1
                break

    fig, ax = plt.subplots(figsize=(10, 6))
    n_at = len(attractor_labels)
    width = 0.18
    x = np.arange(len(degree_buckets)) * 1.4
    colors_list = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#9b59b6', '#34495e', '#95a5a6']

    for ai, (at, name) in enumerate(attractor_labels):
        vals = [counts[at][b] for b in degree_buckets]
        positions_x = x + (ai - (n_at - 1) / 2) * width
        bars = ax.bar(positions_x, vals, width,
                color=colors_list[ai], label=f'{at} - {name}')
        for rect, v in zip(bars, vals):
            if v > 0:
                ax.text(rect.get_x() + rect.get_width() / 2, v + 0.3, str(int(v)),
                        ha='center', va='bottom', fontsize=10, color='black')

    ax.set_xticks(x)
    labels_xticks = [
        f"{b}\nstd range\n[{bucket_edges[i]:.2f}, {bucket_edges[i+1]:.2f})" if i < 3
        else f"{b}\nstd ≥ {bucket_edges[i]:.2f}"
        for i, b in enumerate(degree_buckets)
    ]
    ax.set_xticklabels(labels_xticks)
    ax.set_xlabel('std of out-degree distribution (quartile bins)')
    ax.set_ylabel('world count (integer)')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(top=max(max(c.values()) for c in counts.values()) * 1.15 + 1)
    ax.set_title(f'Topology Realization / Out-Degree Std vs Stability\n(std = std of per-gene out-degree; higher std = more hub-like regulatory topology)')
    ax.legend(fontsize=8)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'out_degree_distribution_vs_stability.png'), dpi=150)
    plt.close(fig)


def plot_topology_exemplars_in_degree(all_results: List[Dict[str, Any]]):
    """Draw directed graph structure for 2 worlds per in-degree std quartile (Q1-Q4)."""
    try:
        import networkx as nx
    except ImportError:
        print('  skip topology exemplars: networkx not installed')
        return

    world_data = []
    for r in all_results:
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        if not os.path.exists(meta_path):
            continue
        with open(meta_path) as f:
            meta = json.load(f)
        in_deg = meta.get('in_degrees', {})
        if not in_deg:
            continue
        values = np.array([int(v) for v in in_deg.values()], dtype=float)
        std = float(np.std(values, ddof=0))
        world_data.append((r['world_id'], std, meta, r.get('primary_attractor', '?')))

    if not world_data:
        return

    stds = np.array([w[1] for w in world_data])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])

    def bucket(s):
        if s < q25:    return 'Q1 (low)'
        if s < q50:    return 'Q2'
        if s < q75:    return 'Q3'
        return 'Q4 (high)'

    quartiles = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)']
    bucket_edges = [0.0, q25, q50, q75, np.inf]

    # For each quartile, pick 2 worlds:
    #   1st: a Type A (representative stable attractor)
    #   2nd: prefer C/D/E (oscillation/divergence), then fall back to any non-A
    # Dedup by topo_idx: each combo (regime, sign_ratio) re-uses the same topology
    # realization, so the underlying graph structure is identical across regimes.
    # Two passes: first try C/D/E for slot 1, then fall back to any non-A.
    selected = {q: [] for q in quartiles}
    for wid, std, meta, atype in world_data:
        b = bucket(std)
        topo_idx = meta.get('topo_idx')
        seen_topo = {m.get('topo_idx') for _, m, _ in selected[b]}
        if topo_idx in seen_topo:
            continue
        if len(selected[b]) == 0 and atype != 'Type A':
            continue
        if len(selected[b]) >= 2:
            continue
        if len(selected[b]) == 1 and atype == 'Type A':
            continue
        # Pass 1 (preferred): slot 1 must be C/D/E
        if len(selected[b]) == 1 and atype not in ('Type C', 'Type D', 'Type E'):
            continue
        selected[b].append((wid, meta, atype))
        if all(len(selected[q]) >= 2 for q in quartiles):
            break
    # Pass 2: relax to any non-A for any quartile still missing slot 1
    if not all(len(selected[q]) >= 2 for q in quartiles):
        for wid, std, meta, atype in world_data:
            b = bucket(std)
            topo_idx = meta.get('topo_idx')
            seen_topo = {m.get('topo_idx') for _, m, _ in selected[b]}
            if topo_idx in seen_topo:
                continue
            if len(selected[b]) >= 2:
                continue
            if len(selected[b]) == 1 and atype == 'Type A':
                continue
            selected[b].append((wid, meta, atype))
            if all(len(selected[q]) >= 2 for q in quartiles):
                break

    fig, axes = plt.subplots(4, 2, figsize=(10, 16))
    for row, q in enumerate(quartiles):
        edges = bucket_edges[row:row + 2]
        if row < 3:
            range_label = f'std ∈ [{edges[0]:.2f}, {edges[1]:.2f})'
        else:
            range_label = f'std ≥ {edges[0]:.2f}'
        for col in range(2):
            ax = axes[row, col]
            if col >= len(selected[q]):
                ax.axis('off')
                continue
            wid, meta, atype = selected[q][col]
            in_deg_meta = meta.get('in_degrees', {})
            std_local = float(np.std([int(v) for v in in_deg_meta.values()], ddof=0))
            G = nx.DiGraph()
            G.add_nodes_from(range(dma.G))
            P_graph = meta.get('P_graph', {})
            edge_signs = meta.get('edge_signs', {})
            for reg_str, targets in P_graph.items():
                reg = int(reg_str)
                for t in targets:
                    s = edge_signs.get(reg_str, {}).get(str(t), 1)
                    G.add_edge(reg, t, sign=s)
            n_edges = G.number_of_edges()
            print(f'  {q} col{col}  {wid}  ({atype})  edges={n_edges}')

            in_deg = meta.get('in_degrees', {})
            in_deg_arr = np.array([int(in_deg.get(str(i), 0)) for i in range(dma.G)])
            sizes = 200 + 80 * in_deg_arr
            pos = nx.circular_layout(G)
            edge_color = []
            for u, v, d in G.edges(data=True):
                edge_color.append('#e74c3c' if d.get('sign', 1) < 0 else '#2c3e50')
            nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes,
                                   node_color='#3498db', edgecolors='white', linewidths=1.2)
            nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_color,
                                   arrows=True, arrowsize=10, width=1.1,
                                   connectionstyle='arc3,rad=0.08')
            nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_color='white',
                                    font_weight='bold')
            ax.set_title(
                f'{q} — {wid} ({atype})\n{range_label}  std={std_local:.2f}  edges={n_edges}',
                fontsize=9)
            ax.axis('off')

    legend_elems = [
        plt.Line2D([0], [0], color='#2c3e50', lw=2, label='activation'),
        plt.Line2D([0], [0], color='#e74c3c', lw=2, label='repression'),
        plt.scatter([], [], s=200, c='#3498db', label='gene (node size ∝ in-degree)'),
    ]
    fig.legend(handles=legend_elems, loc='upper center', ncol=3, fontsize=10,
               bbox_to_anchor=(0.5, 0.995))
    fig.suptitle(
        'Topology Exemplars by In-Degree Std Quartile\n'
        '(each row = one quartile; node size ∝ in-degree; red = repression, dark = activation)',
        fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    fig.savefig(os.path.join(FIGURES_DIR, 'topology_exemplars_in_degree.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  topology exemplars saved: topology_exemplars_in_degree.png')


def plot_topology_exemplars_out_degree(all_results: List[Dict[str, Any]]):
    """Draw directed graph structure for 2 worlds per out-degree std quartile (Q1-Q4)."""
    try:
        import networkx as nx
    except ImportError:
        print('  skip out-degree topology exemplars: networkx not installed')
        return

    world_data = []
    for r in all_results:
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        if not os.path.exists(meta_path):
            continue
        with open(meta_path) as f:
            meta = json.load(f)
        out_deg = meta.get('out_degrees', {})
        if not out_deg:
            continue
        values = np.array([int(v) for v in out_deg.values()], dtype=float)
        std = float(np.std(values, ddof=0))
        world_data.append((r['world_id'], std, meta, r.get('primary_attractor', '?')))

    if not world_data:
        return

    stds = np.array([w[1] for w in world_data])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])

    def bucket(s):
        if s < q25:    return 'Q1 (low)'
        if s < q50:    return 'Q2'
        if s < q75:    return 'Q3'
        return 'Q4 (high)'

    quartiles = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)']
    bucket_edges = [0.0, q25, q50, q75, np.inf]

    selected = {q: [] for q in quartiles}
    for wid, std, meta, atype in world_data:
        b = bucket(std)
        topo_idx = meta.get('topo_idx')
        seen_topo = {m.get('topo_idx') for _, m, _ in selected[b]}
        if topo_idx in seen_topo:
            continue
        if len(selected[b]) == 0 and atype != 'Type A':
            continue
        if len(selected[b]) >= 2:
            continue
        if len(selected[b]) == 1 and atype == 'Type A':
            continue
        if len(selected[b]) == 1 and atype not in ('Type C', 'Type D', 'Type E'):
            continue
        selected[b].append((wid, meta, atype))
        if all(len(selected[q]) >= 2 for q in quartiles):
            break
    if not all(len(selected[q]) >= 2 for q in quartiles):
        for wid, std, meta, atype in world_data:
            b = bucket(std)
            topo_idx = meta.get('topo_idx')
            seen_topo = {m.get('topo_idx') for _, m, _ in selected[b]}
            if topo_idx in seen_topo:
                continue
            if len(selected[b]) >= 2:
                continue
            if len(selected[b]) == 1 and atype == 'Type A':
                continue
            selected[b].append((wid, meta, atype))
            if all(len(selected[q]) >= 2 for q in quartiles):
                break

    fig, axes = plt.subplots(4, 2, figsize=(10, 16))
    for row, q in enumerate(quartiles):
        edges = bucket_edges[row:row + 2]
        if row < 3:
            range_label = f'std ∈ [{edges[0]:.2f}, {edges[1]:.2f})'
        else:
            range_label = f'std ≥ {edges[0]:.2f}'
        for col in range(2):
            ax = axes[row, col]
            if col >= len(selected[q]):
                ax.axis('off')
                continue
            wid, meta, atype = selected[q][col]
            out_deg_meta = meta.get('out_degrees', {})
            std_local = float(np.std([int(v) for v in out_deg_meta.values()], ddof=0))
            G = nx.DiGraph()
            G.add_nodes_from(range(dma.G))
            P_graph = meta.get('P_graph', {})
            edge_signs = meta.get('edge_signs', {})
            for reg_str, targets in P_graph.items():
                reg = int(reg_str)
                for t in targets:
                    s = edge_signs.get(reg_str, {}).get(str(t), 1)
                    G.add_edge(reg, t, sign=s)
            n_edges = G.number_of_edges()
            print(f'  {q} col{col}  {wid}  ({atype})  edges={n_edges}')

            out_deg = meta.get('out_degrees', {})
            out_deg_arr = np.array([int(out_deg.get(str(i), 0)) for i in range(dma.G)])
            sizes = 200 + 80 * out_deg_arr
            pos = nx.circular_layout(G)
            edge_color = []
            for u, v, d in G.edges(data=True):
                edge_color.append('#e74c3c' if d.get('sign', 1) < 0 else '#2c3e50')
            nx.draw_networkx_nodes(G, pos, ax=ax, node_size=sizes,
                                   node_color='#3498db', edgecolors='white', linewidths=1.2)
            nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_color,
                                   arrows=True, arrowsize=10, width=1.1,
                                   connectionstyle='arc3,rad=0.08')
            nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_color='white',
                                    font_weight='bold')
            ax.set_title(
                f'{q} — {wid} ({atype})\n{range_label}  std={std_local:.2f}  edges={n_edges}',
                fontsize=9)
            ax.axis('off')

    legend_elems = [
        plt.Line2D([0], [0], color='#2c3e50', lw=2, label='activation'),
        plt.Line2D([0], [0], color='#e74c3c', lw=2, label='repression'),
        plt.scatter([], [], s=200, c='#3498db', label='gene (node size ∝ out-degree)'),
    ]
    fig.legend(handles=legend_elems, loc='upper center', ncol=3, fontsize=10,
               bbox_to_anchor=(0.5, 0.995))
    fig.suptitle(
        'Topology Exemplars by Out-Degree Std Quartile\n'
        '(each row = one quartile; node size ∝ out-degree; red = repression, dark = activation)',
        fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout()
    plt.subplots_adjust(top=0.95)
    fig.savefig(os.path.join(FIGURES_DIR, 'topology_exemplars_out_degree.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('  topology exemplars saved: topology_exemplars_out_degree.png')


def write_topo_analysis_json(topo_analysis: Dict[str, Any]):
    path = os.path.join(ANALYSIS_DIR, 'topo_strength_analysis.json')
    with open(path, 'w') as f:
        json.dump(topo_analysis, f, indent=2, default=str)


def run_perturbation_analysis():
    os.makedirs(PERTURBATIONS_DIR, exist_ok=True)
    summary = []

    for wid, original_type, cell_idx in PERTURBATION_SELECTION:
        cell_label = f'cell{cell_idx:02d}'
        print(f'\n  perturbation: {wid} / {cell_label} (original: {original_type})')

        meta_path = os.path.join(METADATA_DIR, f'{wid}.json')
        with open(meta_path) as f:
            meta = json.load(f)

        w = dma.World(meta['seed'])
        w.from_dict(meta)

        traj_path = os.path.join(TRAJECTORIES_DIR, f'{wid}_{cell_label}.pt')
        orig_traj = torch.load(traj_path, weights_only=False)
        X_at_intervention = orig_traj['X_traj'][PERTURBATION_TIME]
        knockdown_gene = int(torch.argmax(X_at_intervention).item())
        max_expr_before = float(X_at_intervention[knockdown_gene].item())
        print(f'    highest expr at t={PERTURBATION_TIME}: '
              f'gene {knockdown_gene} = {max_expr_before:.4f}')

        cell_seed = meta['seed'] + cell_idx + 1
        X0 = dma.sample_initial_state(cell_seed)

        result = dma.simulate_single_cell(
            w, X0, dma.T_SIM,
            intervention_time=PERTURBATION_TIME,
            intervention_config={'knockdown': [knockdown_gene]},
        )

        X_traj_full = result['X_traj']
        clip_count_full = result['clip_count']

        # Only analyze post-perturbation trajectory (t >= PERTURBATION_TIME)
        X_traj = X_traj_full[PERTURBATION_TIME:]
        clip_count = clip_count_full[PERTURBATION_TIME:]
        total_clips = int(clip_count.sum().item())

        world_dict = meta.copy()
        world_dict['world_id'] = wid

        analysis = analyze_trajectory(X_traj, clip_count, total_clips, world_dict)
        new_type = analysis['attractor_type']
        print(f'    attractor: {original_type} -> {new_type}')

        perturb_history = {
            'time': PERTURBATION_TIME,
            'type': 'knockdown',
            'knockdown_genes': [knockdown_gene],
        }

        eq_conv = analysis['equilibrium']['converged']
        eq_time = analysis['equilibrium']['convergence_time']
        recovery_time = PERTURBATION_TIME + eq_time if eq_conv else None
        recovery_failure = not eq_conv

        out_path = os.path.join(PERTURBATIONS_DIR, f'{wid}_{cell_label}_perturb.pt')
        torch.save({
            'X_traj': result['X_traj'],
            'clip_count': result['clip_count'],
            'total_clips': result['total_clips'],
            'world_seed': meta['seed'],
            'cell_seed': cell_seed,
            'world_id': wid,
            'original_attractor': original_type,
            'perturbed_attractor': new_type,
            'perturbation_history': perturb_history,
            'knockdown_gene': knockdown_gene,
            'knockdown_gene_expr_at_intervention': max_expr_before,
            'intervention_time': PERTURBATION_TIME,
            'perturbation_type': 'knockdown',
        }, out_path)

        summary.append({
            'world_id': wid,
            'cell': cell_label,
            'original_attractor': original_type,
            'perturbed_attractor': new_type,
            'attractor_changed': new_type != original_type,
            'perturbation_history': perturb_history,
            'knockdown_gene': knockdown_gene,
            'knockdown_gene_expr_at_intervention': max_expr_before,
            'intervention_time': PERTURBATION_TIME,
            'recovery_time': recovery_time,
            'recovery_failure': recovery_failure,
            'equilibrium_converged': eq_conv,
            'equilibrium_time': eq_time,
            'oscillation_type': analysis['oscillation']['oscillation_type'],
            'divergent': analysis['stability']['divergence_existence'],
            'collapsed': analysis['stability']['numerical_collapse'],
            'max_expression': analysis['stability']['max_expression'],
            'final_mean_expression': analysis['stability']['final_mean_expression'],
        })

    summary_path = os.path.join(ANALYSIS_DIR, 'perturbation_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)

    changed_count = sum(1 for s in summary if s['attractor_changed'])
    print(f'\n  perturbation done: {changed_count}/{len(summary)} attractor types changed')


def plot_perturbation_trajectories():
    pert_fig_dir = os.path.join(FIGURES_DIR, 'perturbation')
    os.makedirs(pert_fig_dir, exist_ok=True)
    grouped = {}
    for wid, atype, cell_idx in PERTURBATION_SELECTION:
        grouped.setdefault(atype, []).append((wid, cell_idx))

    for atype, entries in grouped.items():
        at_name = ATTRACTOR_NAMES.get(atype, atype)
        n_worlds = len(entries)
        rows = n_worlds
        cols = 2

        fig, axes = plt.subplots(rows, cols, figsize=(14, 4 * rows))
        if rows == 1:
            axes = axes.reshape(1, 2)

        for r, (wid, cell_idx) in enumerate(entries):
            cell_label = f'cell{cell_idx:02d}'

            traj_path = os.path.join(TRAJECTORIES_DIR, f'{wid}_{cell_label}.pt')
            orig = torch.load(traj_path, weights_only=False)
            orig_X = orig['X_traj'].numpy()
            orig_at = orig.get('attractor_type', atype)

            pert_path = os.path.join(PERTURBATIONS_DIR, f'{wid}_{cell_label}_perturb.pt')
            pert = torch.load(pert_path, weights_only=False)
            pert_X = pert['X_traj'].numpy()
            pert_at = pert.get('perturbed_attractor', atype)

            ax_left = axes[r, 0]
            for g in range(dma.G):
                ax_left.plot(orig_X[:, g], alpha=0.9, linewidth=0.8)
            ax_left.axvline(x=PERTURBATION_TIME, color='gray', linestyle='--', linewidth=0.8, alpha=0.6)
            ax_left.set_title(
                f'{wid}  [{orig_at}]  —  before', fontsize=10)
            ax_left.tick_params(labelsize=7)
            ax_left.spines['top'].set_visible(False)
            ax_left.spines['right'].set_visible(False)
            if r == 0:
                ax_left.set_title(
                    f'{wid}  [{orig_at}]  —  before\n({atype}: {at_name})',
                    fontsize=10)

            ax_right = axes[r, 1]
            for g in range(dma.G):
                ax_right.plot(pert_X[:, g], alpha=0.9, linewidth=0.8)
            ax_right.axvline(x=PERTURBATION_TIME, color='red', linestyle='--', linewidth=0.8, alpha=0.6)
            ax_right.set_title(
                f'{wid}  [{pert_at}]  —  after knockdown', fontsize=10)
            ax_right.tick_params(labelsize=7)
            ax_right.spines['top'].set_visible(False)
            ax_right.spines['right'].set_visible(False)
            if r == 0:
                kd_gene = pert.get('knockdown_gene', '?')
                ax_right.set_title(
                    f'{wid}  [{pert_at}]  —  after knockdown\n'
                    f'(t={PERTURBATION_TIME} knockdown gene {kd_gene})',
                    fontsize=10)

        fig.suptitle(
            f'Perturbation Comparison — {atype}: {at_name}',
            fontsize=13, fontweight='bold')
        fig.tight_layout()
        fname = f'perturbation_{atype.replace(" ", "_").lower()}.png'
        fig.savefig(os.path.join(pert_fig_dir, fname), dpi=150)
        plt.close(fig)
        print(f'  perturbation figure saved: {fname}')


def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
    os.makedirs(METADATA_DIR, exist_ok=True)
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    os.makedirs(PERTURBATIONS_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)

    all_results: List[Dict[str, Any]] = []

    for sign_key, sign_ratio in SIGN_RATIOS.items():
        print(f"--- sign_ratio: {sign_key} ({sign_ratio}) ---")
        results = simulate_topo_stratified(
            sign_key, sign_ratio, save_trajectories=True,
        )
        all_results.extend(results)

    print(f"\ntotal worlds: {len(all_results)}")

    aggregation = aggregate_across_combos(all_results)
    topo_analysis = analyze_topo_strength_interaction(all_results)
    print(f"topo-strength interaction: {topo_analysis['n_topo_diverse']}/{topo_analysis['n_topos']} "
          f"topos change regime ({topo_analysis['topo_diverse_frac']:.1%})")

    rows = build_world_summary_table(all_results)
    write_world_summary_csv(rows)
    write_regime_summary_csv(aggregation, all_results)
    write_analysis_json(all_results, aggregation)
    write_topo_analysis_json(topo_analysis)

    plot_canonical_trajectories(all_results)
    plot_spectral_radius(all_results)
    plot_convergence_distribution(all_results)
    plot_clipping_distribution(all_results)
    plot_topo_strength_interaction(topo_analysis)
    plot_repression_ratio_vs_stability(all_results)
    plot_in_degree_distribution_vs_stability(all_results)
    plot_out_degree_distribution_vs_stability(all_results)
    plot_topology_exemplars_in_degree(all_results)
    plot_topology_exemplars_out_degree(all_results)
    plot_clipping_dominated_by_regime(all_results)
    plot_all_worlds_trajectories()

    run_perturbation_analysis()

    plot_perturbation_trajectories()

    print("done")
    print(f"  figures: {FIGURES_DIR}")
    print(f"  tables:  {TABLES_DIR}")


if __name__ == '__main__':
    main()
