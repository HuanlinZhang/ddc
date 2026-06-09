"""
Run 003 — Model A TF-Restricted Topology Control
=================================================

Per: docs/01_DDC_Lite_Curriculum/run_003/

Steps:
  1. Inherit run_002 parameter grid (4 strength_regimes x 2 sign_ratios)
  2. Sample worlds with TF-restricted topology (ddc_model_a_tf)
  3. Simulate trajectories + per-trajectory analysis
  4. Figures: 10 base (aligned with run_002) + 6 run_003-specific comparisons
  5. Tables: world_summary.tsv (incl. TF columns), regime_summary.tsv

Author: zhanghl
"""
import json
import os
import sys
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
import ddc_model_a_tf as dma

RUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGURES_DIR = os.path.join(RUN_DIR, 'figures')
TABLES_DIR = os.path.join(RUN_DIR, 'tables')
TRAJECTORIES_DIR = os.path.join(RUN_DIR, 'trajectories')
METADATA_DIR = os.path.join(RUN_DIR, 'world_metadata')
PERTURBATIONS_DIR = os.path.join(RUN_DIR, 'perturbations')
ANALYSIS_DIR = os.path.join(RUN_DIR, 'analysis')
SUMMARY_DIR = os.path.join(RUN_DIR, 'summary')

# Path to run_002 results for comparison
RUN_002_DIR = os.path.join(os.path.dirname(RUN_DIR), 'run_002')
COMPARISON_ENABLED = os.path.isdir(RUN_002_DIR)

N_TOPOS_PER_SIGN: int = 10
N_INIT_PER_WORLD: int = 5
T_SIM: int = dma.T_SIM

EPSILON: float = 1e-4               # convergence detection: ||X(t+1)-X(t)|| < epsilon
SPARSITY_THRESHOLD: float = 1e-3     # equilibrium_sparsity: x_eq_i < 1e-3
CONVERGENCE_WINDOW: int = 50
COLLAPSE_THRESHOLD: float = 1e-3
DIVERGENCE_THRESHOLD: float = 1e3
CLIPPING_FRAC_THRESHOLD: float = 0.1
SLOW_CONVERGENCE_THRESHOLD: int = 200

PERTURBATION_TIME: int = 500
# Will be populated after simulation: one world per attractor type
PERTURBATION_SELECTION: List[Tuple[str, str, int]] = []

STRENGTH_REGIMES = {
    'baseline':         (dma.A_MIN_BASELINE, dma.A_MAX_BASELINE, 'conservative baseline'),
    'chen_moderate':    (dma.A_MIN_CHEN_MODERATE, dma.A_MAX_CHEN_MODERATE, 'Chen-moderate'),
    'stress':           (dma.A_MIN_STRESS, dma.A_MAX_STRESS, 'stress-test'),
    'chen_stress':      (dma.A_MIN_CHEN_STRESS, dma.A_MAX_CHEN_STRESS, 'Chen-stress'),
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

ATTRACTOR_COLORS = {
    'Type A': '#2ecc71',
    'Type B': '#3498db',
    'Type C': '#f39c12',
    'Type D': '#e74c3c',
    'Type E': '#9b59b6',
    'Type F': '#95a5a6',
    'Type G': '#bdc3c7',
}

# ── helpers ──

def make_world_id(strength_key: str, sign_key: str, topo_idx: int) -> str:
    return f"{strength_key}_{sign_key}_t{topo_idx:03d}"


def make_combo_key(strength_key: str, sign_key: str) -> str:
    return f"{strength_key}_{sign_key}"


# ── trajectory analysis (identical to run_002) ──

def detect_equilibrium(X_traj):
    """Detect convergence & compute equilibrium stats.

    Per 01_Analysis_Plan §7.1:
      x_eq = mean X[t_final-window : t_final], window = 50
      equilibrium_magnitude = mean_i(x_eq_i)
      equilibrium_sparsity = fraction of genes with x_eq_i < 1e-3
    """
    t_steps = X_traj.shape[0] - 1
    eq_window = min(CONVERGENCE_WINDOW, X_traj.shape[0])
    for t in range(t_steps - CONVERGENCE_WINDOW + 1):
        ok = True
        for w in range(CONVERGENCE_WINDOW):
            if float(torch.norm(X_traj[t + w + 1] - X_traj[t + w]).item()) >= EPSILON:
                ok = False
                break
        if ok:
            X_eq = X_traj[-eq_window:].mean(dim=0)
            return {
                'converged': True, 'convergence_time': t,
                'equilibrium_magnitude': float(X_eq.mean().item()),
                'equilibrium_sparsity': float((X_eq < SPARSITY_THRESHOLD).sum().item()) / dma.G,
            }
    X_eq = X_traj[-eq_window:].mean(dim=0)
    return {
        'converged': False, 'convergence_time': -1,
        'equilibrium_magnitude': float(X_eq.mean().item()),
        'equilibrium_sparsity': float((X_eq < SPARSITY_THRESHOLD).sum().item()) / dma.G,
    }


def analyze_stability(X_traj, clip_count, total_clips):
    max_expr = float(X_traj.max().item())
    bounded = max_expr < DIVERGENCE_THRESHOLD
    final_mean = float(X_traj[-1].mean().item())
    collapsed = final_mean < COLLAPSE_THRESHOLD
    total_steps = X_traj.shape[0]
    clipping_dominated = total_clips > total_steps * dma.G * CLIPPING_FRAC_THRESHOLD
    div_time = -1
    if not bounded:
        for t in range(X_traj.shape[0]):
            if float(X_traj[t].max().item()) >= DIVERGENCE_THRESHOLD:
                div_time = t
                break
    return {
        'bounded': bounded, 'max_expression': max_expr,
        'final_mean_expression': final_mean, 'numerical_collapse': collapsed,
        'divergence_existence': not bounded,
        'divergence_time': div_time if not bounded else None,
        'clipping_dominated': clipping_dominated,
        'total_clips': total_clips,
    }


def analyze_oscillation(X_traj, converged, conv_time):
    BURN_IN = 200
    if converged and conv_time > BURN_IN + 100:
        X = X_traj[BURN_IN:conv_time].numpy()
    else:
        X = X_traj[BURN_IN:].numpy()
    T, G_dim = X.shape
    if T < 50:
        return {'oscillation_exists': False, 'oscillation_type': 'none',
                'amplitude': 0.0, 'frequency': 0.0, 'damping_rate': None,
                'oscillatory_genes': []}
    oscillatory = []; amplitudes = []; freqs = []; damping_rates = []
    for g in range(G_dim):
        x = X[:, g]; signal_range = float(x.max() - x.min())
        if signal_range < EPSILON: continue
        extrema = np.zeros(T, dtype=np.int8)
        for t in range(1, T - 1):
            if   x[t] > x[t-1] and x[t] > x[t+1]: extrema[t] = 1
            elif x[t] < x[t-1] and x[t] < x[t+1]: extrema[t] = -1
        ext_idx = np.where(extrema != 0)[0]
        if len(ext_idx) < 3: continue
        peak_pairs = [abs(float(x[e2] - x[e1])) for e1, e2 in zip(ext_idx[:-1], ext_idx[1:])]
        if len(peak_pairs) < 2: continue
        gene_amp = float(np.median(peak_pairs))
        if signal_range > EPSILON and gene_amp / signal_range < 0.01: continue
        gene_freq = len(ext_idx) / (2.0 * T)
        if len(peak_pairs) >= 4:
            mid = len(peak_pairs)//2
            early = np.mean(peak_pairs[:mid]); late = np.mean(peak_pairs[mid:])
            damping = float((early - late)/early) if early > EPSILON else 0.0
        else: damping = 0.0
        oscillatory.append(g); amplitudes.append(gene_amp)
        freqs.append(gene_freq); damping_rates.append(damping)
    if not oscillatory:
        return {'oscillation_exists': False, 'oscillation_type': 'none',
                'amplitude': 0.0, 'frequency': 0.0, 'damping_rate': None,
                'oscillatory_genes': []}
    avg_damp = float(np.median(damping_rates))
    osc_type = 'damped' if avg_damp > 0.05 else 'sustained'
    return {
        'oscillation_exists': True, 'oscillation_type': osc_type,
        'amplitude': float(np.mean(amplitudes)),
        'frequency': float(np.mean(freqs)), 'damping_rate': avg_damp,
        'oscillatory_genes': oscillatory,
    }


def classify_attractor(eq, st, osc):
    if st['divergence_existence']: return 'Type E'
    if st['numerical_collapse']:   return 'Type F'
    if osc['oscillation_type'] == 'sustained': return 'Type D'
    if osc['oscillation_type'] == 'damped':    return 'Type C'
    if eq['converged']:
        return 'Type A' if eq['convergence_time'] <= SLOW_CONVERGENCE_THRESHOLD else 'Type B'
    # Not converged yet, but bounded / non-diverging / non-collapsing / non-oscillating
    # → slow convergence toward equilibrium (Type B)
    if st['bounded'] and not st['divergence_existence'] and not st['numerical_collapse'] \
       and not osc['oscillation_exists']:
        return 'Type B'
    return 'Type G'


def compute_spectral_info(world_dict):
    A = np.zeros((dma.G, dma.G))
    delta = world_dict['parameters']['delta']
    for i in range(dma.G): A[i, i] = 1.0 - delta[i]
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
    spectral_radius = float(np.max(np.abs(eigvals)))
    eig_dict = {
        'real_parts': [float(ev.real) for ev in eigvals],
        'imag_parts': [float(ev.imag) for ev in eigvals],
        'abs_values': [float(abs(ev)) for ev in eigvals],
        'spectral_radius': spectral_radius,
    }
    return spectral_radius, eig_dict


def analyze_trajectory(X_traj, clip_count, total_clips, world_dict):
    eq = detect_equilibrium(X_traj)
    st = analyze_stability(X_traj, clip_count, total_clips)
    osc = analyze_oscillation(X_traj, eq['converged'], eq['convergence_time'])
    at = classify_attractor(eq, st, osc)
    sr, eig = compute_spectral_info(world_dict)
    return {
        'equilibrium': eq, 'stability': st, 'oscillation': osc,
        'attractor_type': at, 'spectral_radius': sr, 'eigenvalues': eig,
    }


# ── simulation ──

def simulate_topo_stratified(sign_key, sign_ratio, save_trajectories=True):
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
            world_dict['runtime_version'] = 'v0.1-tf'
            world_dict['paired_run002_world_id'] = world_id  # same combo+topo pairing

            # TF metadata
            world_dict.update(dma.tf_metadata(world))

            os.makedirs(METADATA_DIR, exist_ok=True)
            with open(os.path.join(METADATA_DIR, f'{world_id}.json'), 'w') as f:
                json.dump(world_dict, f, indent=2)
                # ensure float values are serializable
                pass

            per_world = []
            for ci in range(N_INIT_PER_WORLD):
                cell_seed = world_seed + 1 + ci
                X0 = dma.sample_initial_state(cell_seed)
                traj = dma.simulate_single_cell(world, X0, t_steps=T_SIM)
                X_traj = traj['X_traj']; clip_count = traj['clip_count']
                total_clips = traj['total_clips']
                analysis = analyze_trajectory(X_traj, clip_count, total_clips, world_dict)
                per_world.append(analysis)
                if save_trajectories:
                    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
                    torch.save({
                        'X_traj': X_traj, 'clip_count': clip_count,
                        'total_clips': total_clips, 'world_id': world_id,
                        'world_seed': world_seed, 'cell_seed': cell_seed,
                        'attractor_type': analysis['attractor_type'],
                    }, os.path.join(TRAJECTORIES_DIR, f'{world_id}_cell{ci:02d}.pt'))

            cell_attractors = [a['attractor_type'] for a in per_world]
            attractor_counts = Counter(cell_attractors)
            primary = max(attractor_counts, key=attractor_counts.get)

            summary = {
                'world_id': world_id, 'world_seed': world_seed,
                'combo_key': combo_key, 'topo_idx': ti,
                'strength_key': strength_key, 'sign_key': sign_key,
                'a_min': a_min, 'a_max': a_max, 'sign_ratio': sign_ratio,
                'per_cell': per_world,
                'primary_attractor': primary,
                'attractor_distribution': dict(attractor_counts),
                'spectral_radius': per_world[0]['spectral_radius'],
                'eigenvalues': per_world[0]['eigenvalues'],
                'converged': all(a['equilibrium']['converged'] for a in per_world),
                'mean_convergence_time': float(np.mean(
                    [a['equilibrium']['convergence_time'] for a in per_world]
                )) if all(a['equilibrium']['converged'] for a in per_world) else None,
                'clipping_dominated': any(
                    a['stability']['clipping_dominated'] for a in per_world),
            }
            results.append(summary)

        combo_results = [r for r in results if r['combo_key'] == combo_key]
        type_dist = Counter(r['primary_attractor'] for r in combo_results)
        print(f"    {N_TOPOS_PER_SIGN} worlds: {dict(type_dist)}")

    return results


# ── aggregation ──

def aggregate_across_combos(all_results):
    combo_groups = defaultdict(list)
    for r in all_results:
        combo_groups[r['combo_key']].append(r)
    agg = {}
    for ck, worlds in combo_groups.items():
        n = len(worlds)
        ac = Counter(w['primary_attractor'] for w in worlds)
        agg[ck] = {
            'n_worlds': n,
            'attractor_distribution': dict(ac),
            'attractor_distribution_pct': {k: round(v/n*100,1) for k,v in ac.items()},
            'mean_spectral_radius': float(np.mean([w['spectral_radius'] for w in worlds])),
            'convergence_rate': float(sum(1 for w in worlds if w['converged'])/n),
            'bounded_rate': float(sum(1 for w in worlds
                if all(a['stability']['bounded'] for a in w['per_cell']))/n),
            'collapse_rate': float(sum(1 for w in worlds
                if any(a['stability']['numerical_collapse'] for a in w['per_cell']))/n),
        }
    return agg


# ── tables ──

def build_world_summary_table(all_results):
    rows = []
    for r in all_results:
        cells = r['per_cell']
        n_cells = len(cells)
        conv_cells = [a for a in cells if a['equilibrium']['converged']]
        # Load TF metadata
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        tf_meta = {}
        if os.path.exists(meta_path):
            with open(meta_path) as f: tf_meta = json.load(f)
        rows.append({
            # ── doc §7.1 core columns (in spec order) ──
            'world_id': r['world_id'],
            'topology_type': 'TF-restricted',
            'N_TF': tf_meta.get('N_TF', dma.N_TF),
            'TF_candidate_edge_count': tf_meta.get('tf_candidate_pool_size', dma.N_TF_CANDIDATES),
            'sampled_edge_count': tf_meta.get('sampled_edge_count', dma.N_EDGES),
            'edge_density_r': dma.EDGE_DENSITY_R,
            'tf_pool_internal_density': tf_meta.get('tf_pool_internal_density', ''),
            'sign_ratio': r['sign_ratio'],
            'decay_δ': dma.DEFAULT_DELTA,
            'basal_transcription_b': dma.DEFAULT_B,
            'parameter_regime': r['combo_key'],
            'run_002_comparison_label': f'edge-density-matched (N_EDGES={dma.N_EDGES})',
            'paired_run002_world_id': r['world_id'],
            'n_replicates': n_cells,
            'edge_density_matched_to_run002': True,
            'topo_idx': r.get('topo_idx', ''),
            'a_min': r['a_min'],
            'a_max': r['a_max'],
            'primary_attractor': r['primary_attractor'],
            'attractor_distribution': r.get('attractor_distribution', {}),
            'convergence_time': r.get('mean_convergence_time'),
            'oscillation_exists': any(a['oscillation']['oscillation_exists'] for a in cells),
            'clipping_freq': round(np.mean([a['stability']['total_clips']
                                            for a in cells]) / (T_SIM * dma.G), 4),
            'clipping_dominated': any(a['stability']['clipping_dominated'] for a in cells),
            'equilibrium_magnitude': round(np.mean([a['equilibrium']['equilibrium_magnitude']
                                           for a in cells]), 4),
            'equilibrium_sparsity': round(np.mean([a['equilibrium']['equilibrium_sparsity']
                                          for a in cells]), 4),
            'spectral_radius': round(r['spectral_radius'], 6),
            'converged': r['converged'],
            'bounded': all(a['stability']['bounded'] for a in cells),
            'max_expression': round(max(a['stability']['max_expression'] for a in cells), 4),
            'numerical_collapse': any(a['stability']['numerical_collapse'] for a in cells),
        })
    return rows


def write_world_summary_csv(rows):
    if not rows: return
    os.makedirs(TABLES_DIR, exist_ok=True)
    path = os.path.join(TABLES_DIR, 'world_summary.tsv')
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_regime_summary_csv(aggregation, all_results):
    os.makedirs(TABLES_DIR, exist_ok=True)
    combo_meta = {}
    for r in all_results:
        ck = r['combo_key']
        if ck not in combo_meta:
            combo_meta[ck] = (r['a_min'], r['a_max'])
    path = os.path.join(TABLES_DIR, 'regime_summary.tsv')

    def _split_combo_key(ck):
        for sk in SIGN_RATIOS:
            if ck.endswith('_' + sk):
                return ck[:-(len(sk)+1)], sk
        return ck, ''

    atypes = ['Type A','Type B','Type C','Type D','Type E','Type F','Type G']
    with open(path, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['strength_regime','sign_ratio_regime','n_worlds','a_min','a_max','act_rep',
                    'decay_delta',
                    'convergence_rate','bounded_rate','collapse_rate',
                    'mean_spectral_radius'] + atypes)
        for ck, agg in sorted(aggregation.items()):
            amin, amax = combo_meta.get(ck, ('',''))
            strength, sign_r = _split_combo_key(ck)
            act_rep = SIGN_RATIO_RATIOS.get(sign_r, '')
            w.writerow([strength, sign_r, agg['n_worlds'], amin, amax, act_rep,
                        dma.DEFAULT_DELTA,
                        round(agg['convergence_rate'],4), round(agg['bounded_rate'],4),
                        round(agg['collapse_rate'],4),
                        round(agg['mean_spectral_radius'],6)] +
                       [agg['attractor_distribution_pct'].get(a,0) for a in atypes])


# ── base plots (aligned with run_002) ──

def plot_canonical_trajectories(all_results):
    exemplars = {at: None for at in ATTRACTOR_NAMES}
    for r in all_results:
        at = r['primary_attractor']
        if exemplars[at] is None:
            exemplars[at] = r['world_id']
    found = {k: v for k, v in exemplars.items() if v is not None}
    if not found: return

    fig, axes = plt.subplots(len(found), 1, figsize=(12, 3*len(found)), squeeze=False)
    for idx, (atype, wid) in enumerate(sorted(found.items())):
        ax = axes[idx][0]
        traj_files = sorted(
            [f for f in os.listdir(TRAJECTORIES_DIR) if f.startswith(wid) and f.endswith('.pt')])
        if traj_files:
            data = torch.load(os.path.join(TRAJECTORIES_DIR, traj_files[0]), weights_only=False)
            X = data['X_traj'].numpy()
            for g in range(dma.G):
                if g in dma.TF_GENES:
                    ax.plot(X[:, g], color='red', alpha=0.8, linewidth=1.0, label='TF' if g == 0 else '')
                else:
                    ax.plot(X[:, g], color='grey', alpha=0.35, linewidth=0.8, label='non-TF' if g == dma.N_TF else '')
            ax.legend(loc='upper right', fontsize=10)
            ax.set_title(f'{atype} — {ATTRACTOR_NAMES[atype]} — {wid}')
            ax.set_xlabel('t'); ax.set_ylabel('X')
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
            fig, axes = plt.subplots(rows, cols, figsize=(cols*5, rows*3), squeeze=False)
            for r_idx, ti in enumerate(topo_batch):
                for s_idx, sk in enumerate(strength_order):
                    wid = make_world_id(sk, sign_key, ti)
                    ax = axes[s_idx, r_idx]
                    traj_files = sorted(
                        [f for f in os.listdir(TRAJECTORIES_DIR)
                         if f.startswith(wid) and f.endswith('.pt')])
                    if not traj_files:
                        ax.set_title(f'{wid} — no data', fontsize=11); continue
                    data = torch.load(os.path.join(TRAJECTORIES_DIR, traj_files[0]),
                                      weights_only=False)
                    X = data['X_traj'].numpy()
                    at = data['attractor_type']
                    for g in range(dma.G):
                        if g in dma.TF_GENES:
                            ax.plot(X[:, g], color='red', alpha=0.8, linewidth=1.0,
                                    label='TF' if g == 0 else '')
                        else:
                            ax.plot(X[:, g], color='grey', alpha=0.35, linewidth=0.8,
                                    label='non-TF' if g == dma.N_TF else '')
                    ax.legend(loc='upper right', fontsize=9)
                    if s_idx == 0:
                        ax.set_title(
                            r'$\bf{topology\ ' + f'{ti:03d}' + r'}$' +
                            f'\n[{at}: {ATTRACTOR_NAMES.get(at, at)}]', fontsize=10)
                    else:
                        ax.set_title(f'[{at}: {ATTRACTOR_NAMES.get(at, at)}]', fontsize=10)
                    ax.tick_params(labelsize=7)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
            for s_idx, sk in enumerate(strength_order):
                amin, amax, _ = STRENGTH_REGIMES[sk]
                axes[s_idx, 0].set_ylabel(f'{sk}\na ∈ [{amin}, {amax}]',
                    fontsize=11, rotation=0, ha='right', va='center', labelpad=10)
            fig.suptitle(
                f'World Trajectories (Cell 0) — {SIGN_RATIO_LABELS[sign_key]} '
                f'(activation:repression={SIGN_RATIO_RATIOS[sign_key]})\n'
                f'topo {t_start:03d}-{t_end-1:03d}\n'
                f'(columns = topology, rows = strength regime) [TF-restricted]',
                fontsize=12, fontweight='bold')
            plt.tight_layout()
            fig.savefig(os.path.join(FIGURES_DIR, 'traj',
                        f'all_trajectories_{sign_key}_{fig_idx+1:02d}.png'), dpi=150)
            plt.close(fig)


def plot_spectral_radius(all_results):
    attractor_order = [
        ('Type A', 'Stable equilibrium'), ('Type B', 'Slow convergence'),
        ('Type C', 'Damped oscillation'), ('Type D', 'Sustained oscillation'),
        ('Type E', 'Runaway divergence'), ('Type F', 'Numerical collapse'),
        ('Type G', 'Others'),
    ]
    by_type_raw = defaultdict(list)
    by_type_no_clip = defaultdict(list)
    by_type_clip_count = defaultdict(int)
    for r in all_results:
        at = r['primary_attractor']
        by_type_raw[at].append(r['spectral_radius'])
        if not all(a['stability']['clipping_dominated'] for a in r['per_cell']):
            by_type_no_clip[at].append(r['spectral_radius'])
        if any(a['stability']['clipping_dominated'] for a in r['per_cell']):
            by_type_clip_count[at] += 1

    labels, data, lt = [], [], []
    for at, name in attractor_order:
        if by_type_raw[at]:
            n_total = len(by_type_raw[at])
            n_clip = by_type_clip_count[at]
            labels.append(at); data.append(by_type_raw[at])
            lt.append(f"{at} - {name}\n(n={n_total}, clip={n_clip}/{n_total})")
        if at == 'Type A' and by_type_no_clip['Type A']:
            n_nc = len(by_type_no_clip['Type A'])
            labels.append('Type A*'); data.append(by_type_no_clip['Type A'])
            lt.append(f"Type A* - non-clipping\n(n={n_nc}, clip=0/{n_nc})")
    if not data: return

    fig, ax = plt.subplots(figsize=(10, 6))
    pos = list(range(1, len(labels)+1))
    bp = ax.boxplot(data, positions=pos, tick_labels=lt, patch_artist=True, widths=0.6)
    for p in bp['boxes']: p.set_facecolor('#a8c5e8')
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='ρ=1')
    ax.set_xticklabels(lt, fontsize=7)
    ax.set_ylabel('Spectral Radius ρ(A)')
    ax.set_title('ρ(A) by Attractor Type [TF-restricted]\n'
                 '(* = only worlds with at least one non-clipping cell)')
    ax.legend()
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'spectral_radius.png'), dpi=150)
    plt.close(fig)


def plot_convergence_distribution(all_results):
    sign_colors = {'balanced': '#3498db', 'repression_biased': '#9b59b6'}
    by_sign = defaultdict(list)
    for r in all_results:
        sk = r.get('sign_key')
        if sk not in sign_colors: continue
        if r.get('converged') and r.get('mean_convergence_time') is not None:
            by_sign[sk].append(r['mean_convergence_time'])
    if not by_sign: return

    fig, ax = plt.subplots(figsize=(8, 5))
    for sk, color in sign_colors.items():
        times = sorted(by_sign.get(sk, []))
        if not times: continue
        cdf = np.arange(1, len(times)+1)/len(times)
        ax.plot(times, cdf, color=color, linewidth=2,
                label=f'{sk} (n={len(times)})', marker='o', markersize=3,
                markevery=max(1, len(times)//20))
    ax.set_xlabel('Convergence Time (steps)')
    ax.set_ylabel('Cumulative Proportion')
    ax.set_title('Convergence Time Distribution by Sign Ratio [TF-restricted]')
    ax.legend(fontsize=9); ax.set_ylim(0, 1.05)
    ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))
    ax.set_xlim(left=0)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'convergence_time.png'), dpi=150)
    plt.close(fig)


def plot_clipping_distribution(all_results):
    at_order = [
        ('Type A','Stable equilibrium','#2ecc71'), ('Type B','Slow convergence','#3498db'),
        ('Type C','Damped oscillation','#f39c12'), ('Type D','Sustained oscillation','#e74c3c'),
        ('Type E','Runaway divergence','#9b59b6'), ('Type F','Numerical collapse','#95a5a6'),
        ('Type G','Others','#bdc3c7'),
    ]
    sign_keys = list(SIGN_RATIOS.keys())
    sign_colors = {'balanced': '#3498db', 'repression_biased': '#9b59b6'}
    by_type = {at: defaultdict(list) for at, _, _ in at_order}
    for r in all_results:
        at = r['primary_attractor']; sk = r.get('sign_key')
        if at not in by_type or sk not in sign_keys: continue
        for a in r['per_cell']:
            freq = a['stability']['total_clips'] / ((T_SIM+1)*dma.G)
            by_type[at][sk].append(freq)
    # Only show types that appear in the data
    present = [(at, name, color) for at, name, color in at_order
               if any(by_type[at][sk] for sk in sign_keys)]
    if not present: return
    fig, ax = plt.subplots(figsize=(max(6, len(present)*1.6), 6))
    x = np.arange(len(present)); width = 0.18
    for si, sk in enumerate(sign_keys):
        ps = x + (si - (len(sign_keys)-1)/2)*width
        d = [by_type[at][sk] for at, _, _ in present]
        bp = ax.boxplot(d, positions=ps, widths=width*0.9, patch_artist=True)
        for p in bp['boxes']: p.set_facecolor(sign_colors[sk]); p.set_alpha(0.7)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{at}\n{name}' for at, name, _ in present], fontsize=10)
    ax.set_ylabel('Clipping Frequency')
    ax.set_title('Clipping Frequency Distribution by Attractor Type [TF-restricted]')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=sk) for sk, c in sign_colors.items()], fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'clipping_frequency.png'), dpi=150)
    plt.close(fig)


def plot_strength_vs_regime_type(all_results):
    """Bar plot: attractor type count per strength regime (pooled across sign ratios)."""
    at_labels = [
        ('Type A','Stable'), ('Type B','Slow'), ('Type C','Damped'),
        ('Type D','Sustained'), ('Type E','Divergence'), ('Type F','Collapse'),
        ('Type G','Others'),
    ]
    combos = list(STRENGTH_REGIMES.keys())
    counts = {c: {at: 0 for at, _ in at_labels} for c in combos}
    for r in all_results:
        ck = r.get('strength_key'); at = r['primary_attractor']
        if ck in counts and at in counts[ck]:
            counts[ck][at] += 1

    # Only show types that appear
    present = [(at, name) for at, name in at_labels
               if any(counts[c][at] > 0 for c in combos)]
    if not present: return
    colors = ['#2ecc71','#3498db','#f39c12','#e74c3c','#9b59b6','#34495e','#95a5a6']

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(combos)); width = 0.6
    bottoms = np.zeros(len(combos))
    for ai, (at, name) in enumerate(present):
        vals = np.array([counts[c][at] for c in combos])
        bars = ax.bar(x, vals, width, bottom=bottoms,
                      color=colors[ai], label=f'{at} - {name}')
        for rect, b0, v in zip(bars, bottoms, vals):
            if v > 0:
                ax.text(rect.get_x()+rect.get_width()/2, b0+v/2, str(int(v)),
                        ha='center', va='center', fontsize=10, color='white',
                        fontweight='bold')
        bottoms += vals

    ax.set_xticks(x)
    ax.set_xticklabels([f'{c}\n$A \\in [{STRENGTH_REGIMES[c][0]}, {STRENGTH_REGIMES[c][1]}]$'
                        for c in combos], fontsize=9)
    ax.set_ylabel('World Count')
    ax.set_title('Interaction Strength vs Regime Type [TF-restricted]')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    max_s = max(sum(counts[c].values()) for c in combos)
    ax.set_ylim(0, max_s*1.15)
    ax.legend(fontsize=8, loc='upper right', framealpha=0.9)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'strength_vs_regime_type.png'), dpi=150)
    plt.close(fig)


def plot_repression_ratio_vs_stability(all_results):
    at_labels = [
        ('Type A','Stable'), ('Type B','Slow'), ('Type C','Damped'),
        ('Type D','Sustained'), ('Type E','Divergence'), ('Type F','Collapse'),
        ('Type G','Others'),
    ]
    sign_keys = list(SIGN_RATIOS.keys()); combos = list(STRENGTH_REGIMES.keys())
    counts = {sk: {c: {at: 0 for at, _ in at_labels} for c in combos} for sk in sign_keys}
    for r in all_results:
        sk = r.get('sign_key'); ck = r.get('strength_key'); at = r['primary_attractor']
        if sk in counts and ck in counts[sk] and at in counts[sk][ck]:
            counts[sk][ck][at] += 1
    fig, ax = plt.subplots(figsize=(10, 6))
    n_signs = len(sign_keys); width = 0.18; x = np.arange(len(combos))
    colors_list = ['#2ecc71','#3498db','#f39c12','#e74c3c','#9b59b6','#34495e','#95a5a6']
    sign_short = {'balanced': 'bal', 'repression_biased': 'rep'}
    for si, sk in enumerate(sign_keys):
        bottoms = np.zeros(len(combos))
        for ai, (at, name) in enumerate(at_labels):
            vals = np.array([counts[sk][c][at] for c in combos])
            px = x + (si - (n_signs-1)/2)*width*1.4
            bars = ax.bar(px, vals, width, bottom=bottoms, color=colors_list[ai],
                         label=f'{at} - {name}' if si==0 else None)
            for rect, b0, v in zip(bars, bottoms, vals):
                if v > 0:
                    ax.text(rect.get_x()+rect.get_width()/2, b0+v/2, str(int(v)),
                            ha='center', va='center', fontsize=8, color='white',
                            fontweight='bold')
            bottoms += vals
    for ci, c in enumerate(combos):
        for si, sk in enumerate(sign_keys):
            pos_x = ci + (si-(n_signs-1)/2)*width*1.4
            ax.text(pos_x, -0.6, sign_short.get(sk,sk), ha='center', va='top',
                    fontsize=10, color='gray', fontweight='bold')
    ax.set_xticks(x); ax.set_xticklabels(combos)
    ax.set_xlabel('strength regime'); ax.set_ylabel('world count')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    max_s = max(sum(counts[sk][c].values()) for sk in sign_keys for c in combos)
    ax.set_ylim(bottom=-1.2, top=max_s*1.1+1)
    ax.set_title('Repression Ratio vs Stability [TF-restricted]')
    # Legend: only types that appear
    present_at = [(a,n,c) for (a,n),c in zip(at_labels, colors_list)
                  if any(counts[sk][c2][a] > 0 for sk in sign_keys for c2 in combos)]
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=c, label=f'{a}-{n}') for a,n,c in present_at],
              fontsize=8, ncol=len(present_at), loc='upper left')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'repression_ratio_vs_stability.png'), dpi=150)
    plt.close(fig)


def plot_clipping_dominated_by_regime(all_results):
    sign_keys = list(SIGN_RATIOS.keys()); combos = list(STRENGTH_REGIMES.keys())
    counts = {sk: {c: {'clip':0,'total':0} for c in combos} for sk in sign_keys}
    for r in all_results:
        sk = r.get('sign_key'); ck = r.get('strength_key')
        if sk in counts and ck in counts[sk]:
            counts[sk][ck]['total'] += 1
            if r.get('clipping_dominated'): counts[sk][ck]['clip'] += 1
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(combos)); width = 0.35
    sign_colors = {'balanced': '#3498db', 'repression_biased': '#9b59b6'}
    sign_label = {'balanced': 'balanced (1:1)', 'repression_biased': 'repression_biased (1:2)'}
    for si, sk in enumerate(sign_keys):
        fracs = [counts[sk][c]['clip']/counts[sk][c]['total']
                 if counts[sk][c]['total']>0 else 0.0 for c in combos]
        px = x + (si-0.5)*width
        bars = ax.bar(px, [f*100 for f in fracs], width, color=sign_colors[sk],
                      label=sign_label[sk], edgecolor='white', linewidth=0.5)
        for rect, frac, c in zip(bars, fracs, combos):
            cl = counts[sk][c]['clip']; tot = counts[sk][c]['total']
            ax.text(rect.get_x()+rect.get_width()/2, rect.get_height()+1,
                    f'{frac*100:.0f}%\n({cl}/{tot})', ha='center', va='bottom',
                    fontsize=8)
    ax.set_xticks(x); ax.set_xticklabels(combos)
    ax.set_xlabel('strength regime'); ax.set_ylabel('clipping-dominated worlds (%)')
    ax.set_ylim(0, 115)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{int(y)}%'))
    ax.set_title('Clipping-Dominated Worlds by Regime and Sign Ratio [TF-restricted]')
    ax.legend(fontsize=8, loc='lower right')
    ax.grid(axis='y', linestyle='--', alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'clipping_dominated_by_regime.png'), dpi=150)
    plt.close(fig)


def plot_degree_distribution_vs_stability(all_results):
    at_labels = [
        ('Type A','Stable'), ('Type B','Slow'), ('Type C','Damped'),
        ('Type D','Sustained'), ('Type E','Divergence'), ('Type F','Collapse'),
        ('Type G','Others'),
    ]
    world_data = []
    for r in all_results:
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        if not os.path.exists(meta_path): continue
        with open(meta_path) as f: meta = json.load(f)
        in_deg = meta.get('in_degrees', {})
        if not in_deg: continue
        vals = np.array([int(v) for v in in_deg.values()], dtype=float)
        std = float(np.std(vals, ddof=0))
        world_data.append((std, r['primary_attractor']))
    if not world_data: return

    stds = np.array([w[0] for w in world_data])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])
    edges = [0.0, q25, q50, q75, np.inf]
    def bucket(v):
        if v < q25: return 'Q1 (low)'
        if v < q50: return 'Q2'
        if v < q75: return 'Q3'
        return 'Q4 (high)'
    buckets = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)']
    counts = {at: {b:0 for b in buckets} for at, _ in at_labels}
    for std, at in world_data:
        b = bucket(std)
        for a, _ in at_labels:
            if at == a: counts[a][b] += 1; break

    fig, ax = plt.subplots(figsize=(10, 6))
    n_at = len(at_labels); width = 0.18
    x = np.arange(len(buckets))*1.4
    colors_list = ['#2ecc71','#3498db','#f39c12','#e74c3c','#9b59b6','#34495e','#95a5a6']
    # Only show legend entries for types that appear
    present_ai = []
    present_colors = []
    present_labels = []
    for ai, (at, name) in enumerate(at_labels):
        if any(counts[at][b] > 0 for b in buckets):
            present_ai.append(ai)
            present_colors.append(colors_list[ai])
            present_labels.append(f'{at} - {name}')

    for ai in present_ai:
        at, name = at_labels[ai]
        vals = [counts[at][b] for b in buckets]
        px = x + (ai-(n_at-1)/2)*width
        bars = ax.bar(px, vals, width, color=colors_list[ai])
        for rect, v in zip(bars, vals):
            if v > 0:
                ax.text(rect.get_x()+rect.get_width()/2, v+0.3, str(int(v)),
                        ha='center', va='bottom', fontsize=10, color='black')
    ax.legend(handles=[plt.Rectangle((0,0),1,1, color=c) for c in present_colors],
              labels=present_labels, fontsize=9)
    ax.set_xticks(x)
    xlabels = [f"{b}\nstd range\n[{edges[i]:.2f}, {edges[i+1]:.2f})" if i<3
               else f"{b}\nstd ≥ {edges[i]:.2f}"
               for i, b in enumerate(buckets)]
    ax.set_xticklabels(xlabels)
    ax.set_xlabel('std of in-degree distribution (quartile bins)')
    ax.set_ylabel('world count (integer)')
    ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_ylim(top=max(max(c.values()) for c in counts.values())*1.15+1)
    ax.set_title('Degree Std vs Stability [TF-restricted]\n'
                 '(all regulators are TF genes only)')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'degree_distribution_vs_stability.png'), dpi=150)
    plt.close(fig)


def plot_topology_exemplars(all_results):

    try: import networkx as nx
    except ImportError:
        print('  skip topology exemplars: networkx not installed'); return

    world_data = []
    for r in all_results:
        meta_path = os.path.join(METADATA_DIR, f'{r["world_id"]}.json')
        if not os.path.exists(meta_path): continue
        with open(meta_path) as f: meta = json.load(f)
        in_deg = meta.get('in_degrees', {})
        if not in_deg: continue
        vals = np.array([int(v) for v in in_deg.values()], dtype=float)
        std = float(np.std(vals, ddof=0))
        world_data.append((r['world_id'], std, meta, r.get('primary_attractor', '?')))
    if not world_data: return

    stds = np.array([w[1] for w in world_data])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])
    def bucket(s):
        if s < q25: return 'Q1 (low)'
        if s < q50: return 'Q2'
        if s < q75: return 'Q3'
        return 'Q4 (high)'
    quartiles = ['Q1 (low)', 'Q2', 'Q3', 'Q4 (high)']
    bedges = [0.0, q25, q50, q75, np.inf]

    selected = {q: [] for q in quartiles}
    # pass 1: Type A + C/D/E
    for wid, std, meta, atype in world_data:
        b = bucket(std)
        tidx = meta.get('topo_idx')
        seen = {m.get('topo_idx') for _, m, _ in selected[b]}
        if tidx in seen: continue
        if len(selected[b]) >= 2: continue
        if len(selected[b]) == 0 and atype != 'Type A': continue
        if len(selected[b]) == 1 and atype == 'Type A': continue
        if len(selected[b]) == 1 and atype not in ('Type C','Type D','Type E'): continue
        selected[b].append((wid, meta, atype))
    # pass 2: relax
    if not all(len(selected[q])>=2 for q in quartiles):
        for wid, std, meta, atype in world_data:
            b = bucket(std)
            tidx = meta.get('topo_idx')
            seen = {m.get('topo_idx') for _, m, _ in selected[b]}
            if tidx in seen: continue
            if len(selected[b]) >= 2: continue
            if len(selected[b]) == 1 and atype == 'Type A': continue
            selected[b].append((wid, meta, atype))
            if all(len(selected[q])>=2 for q in quartiles): break

    fig, axes = plt.subplots(4, 2, figsize=(10, 16))
    for row, q in enumerate(quartiles):
        edges_arr = bedges[row:row+2]
        if row < 3:
            range_label = f'std ∈ [{edges_arr[0]:.2f}, {edges_arr[1]:.2f})'
        else:
            range_label = f'std ≥ {edges_arr[0]:.2f}'
        for col in range(2):
            ax = axes[row, col]
            if col >= len(selected[q]): ax.axis('off'); continue
            wid, meta, atype = selected[q][col]
            in_deg_m = meta.get('in_degrees', {})
            std_l = float(np.std([int(v) for v in in_deg_m.values()], ddof=0))

            G = nx.DiGraph()
            G.add_nodes_from(range(dma.G))
            P_graph = meta.get('P_graph', {})
            edge_signs = meta.get('edge_signs', {})
            for tgt_str, regs in P_graph.items():
                tgt = int(tgt_str)
                for r in regs:
                    s = edge_signs.get(tgt_str, {}).get(str(r), 1)
                    G.add_edge(int(r), tgt, sign=s)
            n_edges = G.number_of_edges()

            in_deg = meta.get('in_degrees', {})
            in_deg_arr = np.array([int(in_deg.get(str(i),0)) for i in range(dma.G)])
            sizes = 200 + 80 * in_deg_arr
            pos = nx.circular_layout(G)
            tf_nodes = list(range(dma.N_TF))
            non_tf_nodes = list(range(dma.N_TF, dma.G))
            nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=tf_nodes,
                                   node_size=sizes[:dma.N_TF],
                                   node_color='#E67E22', edgecolors='white', linewidths=1.5)
            nx.draw_networkx_nodes(G, pos, ax=ax, nodelist=non_tf_nodes,
                                   node_size=sizes[dma.N_TF:],
                                   node_color='#3498db', edgecolors='white', linewidths=1.2)
            nx.draw_networkx_labels(G, pos, ax=ax, font_size=10, font_color='white',
                                    font_weight='bold')

            # ── edges: activation (red arrow via nx), repression (black line via nx + T-bar) ──
            # Split edges by sign
            act_edges = [(u,v) for u,v,d in G.edges(data=True) if d.get('sign',1)>0]
            rep_edges = [(u,v) for u,v,d in G.edges(data=True) if d.get('sign',1)<=0]
            # Activation: use networkx's built-in arrow drawing (handles node_size correctly)
            if act_edges:
                nx.draw_networkx_edges(G, pos, ax=ax, edgelist=act_edges,
                                       edge_color='#e74c3c', node_size=sizes,
                                       arrows=True, arrowsize=10, width=1.1,
                                       connectionstyle='arc3,rad=0.08')
            # Repression: draw lines via networkx (no arrows), then overlay T-bar
            if rep_edges:
                nx.draw_networkx_edges(G, pos, ax=ax, edgelist=rep_edges,
                                       edge_color='black', node_size=sizes,
                                       arrows=False, width=1.1,
                                       connectionstyle='arc3,rad=0.08')
            # Add T-bars at repression target endpoints
            ax.autoscale_view()
            ww = ax.bbox.width if ax.bbox.width > 0 else 500
            hh = ax.bbox.height if ax.bbox.height > 0 else 500
            dxr = ax.get_xlim()[1] - ax.get_xlim()[0]
            dyr = ax.get_ylim()[1] - ax.get_ylim()[0]
            ppx = ww / dxr; ppy = hh / dyr
            for u, v in rep_edges:
                rB_pts = np.sqrt(sizes[v] / np.pi)  # true node radius in points
                rB_x = rB_pts / ppx; rB_y = rB_pts / ppy
                x1, y1 = pos[u]; x2, y2 = pos[v]
                dx = x2 - x1; dy = y2 - y1
                d = np.sqrt(dx**2 + dy**2)
                if d > 0:
                    ux, uy = dx/d, dy/d
                    rB = np.sqrt((rB_x*ux)**2 + (rB_y*uy)**2)
                    ex = x2 - dx * rB / d
                    ey = y2 - dy * rB / d
                    # T-bar perpendicular in display coords
                    ux_d = dx * ppx; uy_d = dy * ppy
                    dn = np.sqrt(ux_d**2 + uy_d**2)
                    if dn > 0:
                        pdx = -uy_d/dn; pdy = ux_d/dn
                        bar = 3.0 / ppx  # ~3 pts in data coords
                        ax.plot([ex + bar*pdx, ex - bar*pdx],
                                [ey + bar*pdy, ey - bar*pdy],
                                color='black', linewidth=1.1, zorder=6)
            ax.set_title(f'{q} — {wid} ({atype})\n{range_label}  std={std_l:.2f}  edges={n_edges}',
                         fontsize=9)
            ax.axis('off')

    legend_elems = [
        plt.Line2D([0],[0],color='#e74c3c',lw=2,label='activation'),
        plt.Line2D([0],[0],color='black',lw=2,label='repression'),
        plt.scatter([],[],s=200,c='#E67E22',label='TF gene'),
        plt.scatter([],[],s=200,c='#3498db',label='non-TF gene'),
    ]
    fig.legend(handles=legend_elems, loc='upper center', ncol=4, fontsize=10,
               bbox_to_anchor=(0.5,0.995))
    fig.suptitle('Topology Exemplars by In-Degree Std Quartile [TF-restricted]\n'
                 '(only TF genes 0-4 can be regulators)',
                 fontsize=12, fontweight='bold', y=1.02)
    plt.tight_layout(); plt.subplots_adjust(top=0.95)
    fig.savefig(os.path.join(FIGURES_DIR, 'topology_exemplars.png'), dpi=150, bbox_inches='tight')
    plt.close(fig)


def plot_topo_strength_interaction(all_results):
    """Topology × Strength regime attractor heatmap (adapted for TF-restricted)."""
    summaries = {}
    for r in all_results:
        sign_key = r.get('sign_key',''); ti = r.get('topo_idx','')
        key = f"{sign_key}_t{ti:03d}"
        summaries.setdefault(key, {'regime_map': {}, 'n_worlds': 0})
        summaries[key]['n_worlds'] += 1
        summaries[key]['regime_map'][r['strength_key']] = r['primary_attractor']
    if not summaries: return

    topo_keys = sorted(summaries.keys())
    strength_order = list(STRENGTH_REGIMES.keys())
    sign_keys = list(SIGN_RATIOS.keys())
    topo_nums = sorted({tk.split("_t")[-1] for tk in topo_keys}, key=lambda s: int(s))

    fig, ax = plt.subplots(figsize=(10, 6))
    yticks, ylabels = [], []
    y_offset = 1.0; intra_gap = 0.9; inter_gap = 1.5
    cell_w, cell_h = 0.85, 0.75
    ordered = []

    for si, sk in enumerate(sign_keys):
        if si > 0: y_offset += inter_gap
        for tn in topo_nums:
            key = f'{sk}_t{tn}'
            if key in summaries:
                ordered.append(key)
                yticks.append(y_offset)
                ylabels.append(f'{sk}_topo{tn}')
                y_offset += intra_gap

    for topo_key, y in zip(ordered, yticks):
        regime_map = summaries[topo_key]['regime_map']
        for xi, sk in enumerate(strength_order):
            at = regime_map.get(sk, '?')
            color = ATTRACTOR_COLORS.get(at, '#cccccc')
            ax.barh(y, cell_w, left=xi+(1-cell_w)/2, height=cell_h,
                    color=color, edgecolor='white', linewidth=0.3)

    ax.invert_yaxis()
    ax.set_yticks(yticks); ax.set_yticklabels(ylabels, fontsize=7)
    ax.set_xticks([i+0.5 for i in range(len(strength_order))])
    ax.set_xticklabels([f'{sk}\na ∈ [{STRENGTH_REGIMES[sk][0]}, {STRENGTH_REGIMES[sk][1]}]'
                        for sk in strength_order], rotation=0, fontsize=8, ha='center')
    ax.set_xlim(0, len(strength_order))
    used_a = set()
    for tk in topo_keys:
        for sk in strength_order:
            at = summaries[tk]['regime_map'].get(sk,'?')
            if at in ATTRACTOR_COLORS: used_a.add(at)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=ATTRACTOR_COLORS[a], label=a) for a in sorted(used_a)],
              fontsize=8, ncol=len(used_a), loc='upper right',
              bbox_to_anchor=(0.99,0.99), framealpha=0.9)
    for xi in range(1, len(strength_order)):
        ax.axvline(x=xi, color='white', linewidth=0.5)
    fig.suptitle('Topology × Strength Regime: Attractor Type [TF-restricted]',
                 fontsize=12, fontweight='bold', y=0.97)
    fig.subplots_adjust(left=0.17, right=0.97, top=0.92, bottom=0.08)
    fig.savefig(os.path.join(FIGURES_DIR, 'topo_strength_interaction.png'), dpi=150)
    plt.close(fig)


# ── run_002 comparison plots ──

def load_run002_results():
    """Load run_002 aggregated results for comparison."""
    path = os.path.join(RUN_002_DIR, 'analysis', 'all_world_results.json')
    if not os.path.exists(path): return None
    with open(path) as f:
        return json.load(f)


def plot_run002_vs_run003_attractor(all_results, r2_data):
    """Side-by-side stacked bar: attractor distribution per combo, run_002 vs run_003."""
    if r2_data is None: return

    at_order = ['Type A','Type B','Type C','Type D','Type E','Type F','Type G']
    combos = [make_combo_key(sk, srk) for sk in STRENGTH_REGIMES for srk in SIGN_RATIOS]
    at_colors = ['#2ecc71','#3498db','#f39c12','#e74c3c','#9b59b6','#95a5a6','#bdc3c7']

    # run_003
    r3_counts = {c: Counter() for c in combos}
    for r in all_results:
        r3_counts[r['combo_key']][r['primary_attractor']] += 1
    # run_002
    r2_counts = {c: Counter() for c in combos}
    for r in r2_data:
        ck = r.get('combo_key', ''); at = r.get('primary_attractor', '')
        if ck in r2_counts and at in at_order:
            r2_counts[ck][at] += 1

    n_combos = len(combos)
    n_cols = 4
    n_rows = (n_combos + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols*3.2, n_rows*4.5), sharey=True)
    axes = axes.flatten()
    max_count = 10

    # Determine which attractor types actually appear
    used_at_set = set()
    for ck in combos:
        used_at_set.update(r2_counts[ck].keys())
        used_at_set.update(r3_counts[ck].keys())
    used_at = [at for at in at_order if at in used_at_set]
    used_colors = [at_colors[at_order.index(at)] for at in used_at]

    for ci, ck in enumerate(combos):
        ax = axes[ci]
        x = np.arange(2)
        width = 0.5
        bottoms_r2 = np.zeros(2); bottoms_r3 = np.zeros(2)
        for ai, at in enumerate(at_order):
            v2 = r2_counts[ck].get(at, 0); v3 = r3_counts[ck].get(at, 0)
            if v2 == 0 and v3 == 0:
                continue
            # run_002 bar
            ax.bar(x[0], v2, width, bottom=bottoms_r2[0],
                   color=at_colors[ai], edgecolor='white', linewidth=0.5,
                   label=at if ci == 0 else None)
            if v2 > 0:
                ax.text(x[0], bottoms_r2[0]+v2/2, str(int(v2)), ha='center',
                        va='center', fontsize=9, color='white', fontweight='bold')
            bottoms_r2[0] += v2
            # run_003 bar
            ax.bar(x[1], v3, width, bottom=bottoms_r3[1],
                   color=at_colors[ai], edgecolor='white', linewidth=0.5)
            if v3 > 0:
                ax.text(x[1], bottoms_r3[1]+v3/2, str(int(v3)), ha='center',
                        va='center', fontsize=9, color='white', fontweight='bold')
            bottoms_r3[1] += v3
        ax.set_xticks(x); ax.set_xticklabels(['Run 002','Run 003'], fontsize=11)
        # 2-line title: split combo_key into strength_regime \n sign_ratio
        str_line, sign_line = ck, ''
        for sk in sorted(STRENGTH_REGIMES.keys(), key=len, reverse=True):
            if ck.startswith(sk):
                str_line = sk
                sign_line = ck[len(sk)+1:]
                break
        ax.set_title(f'{str_line}\n{sign_line}', fontsize=10,
                     linespacing=1.3)
        ax.set_ylim(0, max_count + 0.5)
        ax.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        if ci == 0 or ci == n_cols:
            ax.set_ylabel('World count', fontsize=11)
    # Hide unused subplots
    for ci in range(n_combos, len(axes)):
        axes[ci].set_visible(False)
    fig.suptitle('Attractor Distribution: run_002 vs run_003 (TF-restricted)\n'
                 '(same parameter grid, different topology generator)',
                 fontsize=12, fontweight='bold')
    leg_labels = [f'{at} — {ATTRACTOR_NAMES.get(at, at)}' for at in used_at]
    leg_handles = [plt.Rectangle((0,0),1,1, color=c) for c in used_colors]
    # Reorder for column-first fill to appear row-first: [0,2,4,...,1,3,5,...]
    n_at = len(used_at)
    reorder = [i for i in range(0, n_at, 2)] + [i for i in range(1, n_at, 2)]
    leg_handles = [leg_handles[i] for i in reorder]
    leg_labels = [leg_labels[i] for i in reorder]
    fig.legend(handles=leg_handles, labels=leg_labels, fontsize=11, ncol=2,
               loc='lower center', bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(os.path.join(FIGURES_DIR, 'run002_vs_run003_attractor.png'), dpi=150,
                bbox_inches='tight')
    plt.close(fig)


def plot_run002_vs_run003_spectral(all_results, r2_data):
    """Rho distribution comparison per combo."""
    if r2_data is None: return

    combos = [make_combo_key(sk, srk) for sk in STRENGTH_REGIMES for srk in SIGN_RATIOS]
    r3_by_combo = {c: [] for c in combos}
    for r in all_results:
        r3_by_combo[r['combo_key']].append(r['spectral_radius'])
    r2_by_combo = {c: [] for c in combos}
    for r in r2_data:
        ck = r.get('combo_key',''); sr = r.get('spectral_radius')
        if ck in r2_by_combo and sr is not None:
            r2_by_combo[ck].append(sr)

    n = len(combos)
    fig, ax = plt.subplots(figsize=(n*1.3+2, 5.5))
    x = np.arange(n)
    width = 0.35
    for ci, ck in enumerate(combos):
        if r2_by_combo[ck]:
            ax.boxplot(r2_by_combo[ck], positions=[x[ci]-width/2], widths=width*0.8,
                       patch_artist=True, boxprops=dict(facecolor='#a8c5e8'),
                       medianprops=dict(color='#2c3e50'))
        if r3_by_combo[ck]:
            ax.boxplot(r3_by_combo[ck], positions=[x[ci]+width/2], widths=width*0.8,
                       patch_artist=True, boxprops=dict(facecolor='#e8a8a8'),
                       medianprops=dict(color='#8b0000'))
    ax.set_xticks(x)
    xt_labels = []
    for c in combos:
        sl, sgl = c, ''
        for sk in sorted(STRENGTH_REGIMES.keys(), key=len, reverse=True):
            if c.startswith(sk):
                sl = sk; sgl = c[len(sk)+1:]; break
        xt_labels.append(f'{sl}\n{sgl}')
    ax.set_xticklabels(xt_labels, fontsize=8)
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='ρ=1')
    ax.set_ylabel('Spectral Radius ρ(A)')
    ax.set_title('ρ(A) Comparison: run_002 vs run_003 (TF-restricted)')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color='#a8c5e8',label='run_002'),
                       Patch(color='#e8a8a8',label='run_003')], fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'run002_vs_run003_spectral_radius.png'), dpi=150)
    plt.close(fig)


def plot_run002_vs_run003_clipping(all_results, r2_data):
    """Clipping frequency comparison per combo."""
    if r2_data is None: return

    combos = [make_combo_key(sk, srk) for sk in STRENGTH_REGIMES for srk in SIGN_RATIOS]
    r3_clip = {c: [] for c in combos}
    for r in all_results:
        for a in r['per_cell']:
            freq = a['stability']['total_clips'] / ((T_SIM+1)*dma.G)
            r3_clip[r['combo_key']].append(freq)
    r2_clip = {c: [] for c in combos}
    # run_002 clip data from trajectories if available
    if os.path.isdir(os.path.join(RUN_002_DIR, 'trajectories')):
        for r in r2_data:
            ck = r.get('combo_key',''); wid = r.get('world_id','')
            if ck not in r2_clip: continue
            for ci in range(N_INIT_PER_WORLD):
                tp = os.path.join(RUN_002_DIR, 'trajectories', f'{wid}_cell{ci:02d}.pt')
                if os.path.exists(tp):
                    d = torch.load(tp, weights_only=False)
                    freq = d['total_clips'] / (d['X_traj'].shape[0] * dma.G)
                    r2_clip[ck].append(freq)

    n = len(combos)
    fig, ax = plt.subplots(figsize=(n*1.3+2, 5.5))
    x = np.arange(n); width = 0.35
    for ci, ck in enumerate(combos):
        if r2_clip[ck]:
            ax.boxplot(r2_clip[ck], positions=[x[ci]-width/2], widths=width*0.8,
                       patch_artist=True, boxprops=dict(facecolor='#a8c5e8'),
                       medianprops=dict(color='#2c3e50'))
        if r3_clip[ck]:
            ax.boxplot(r3_clip[ck], positions=[x[ci]+width/2], widths=width*0.8,
                       patch_artist=True, boxprops=dict(facecolor='#e8a8a8'),
                       medianprops=dict(color='#8b0000'))
    ax.set_xticks(x)
    xt_labels = []
    for c in combos:
        sl, sgl = c, ''
        for sk in sorted(STRENGTH_REGIMES.keys(), key=len, reverse=True):
            if c.startswith(sk):
                sl = sk; sgl = c[len(sk)+1:]; break
        xt_labels.append(f'{sl}\n{sgl}')
    ax.set_xticklabels(xt_labels, fontsize=8)
    ax.axhline(y=0.1, color='red', linestyle='--', alpha=0.7, label='clip=0.1')
    ax.set_ylabel('Clipping Frequency')
    ax.set_title('Clipping Frequency: run_002 vs run_003 (TF-restricted)')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color='#a8c5e8',label='run_002'),
                       Patch(color='#e8a8a8',label='run_003')], fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'run002_vs_run003_clipping.png'), dpi=150)
    plt.close(fig)


def plot_run002_vs_run003_convergence(all_results, r2_data):
    """Convergence time CDF overlay."""
    if r2_data is None: return

    sign_keys = list(SIGN_RATIOS.keys())
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for si, sk in enumerate(sign_keys):
        ax = axes[si]
        for run_label, results, color, ls in [
            ('run_002', r2_data, '#3498db', '--'),
            ('run_003', all_results, '#e74c3c', '-'),
        ]:
            times = []
            for r in results:
                if r.get('sign_key') != sk: continue
                if r.get('converged') and r.get('mean_convergence_time') is not None:
                    times.append(r['mean_convergence_time'])
            times = sorted(times)
            if not times: continue
            cdf = np.arange(1, len(times)+1)/len(times)
            ax.plot(times, cdf, color=color, linestyle=ls, linewidth=2,
                    label=f'{run_label} (n={len(times)})')
        ax.set_title(f'{SIGN_RATIO_LABELS[sk]} ({SIGN_RATIO_RATIOS[sk]})')
        ax.set_xlabel('Convergence Time (steps)'); ax.set_ylabel('CDF')
        ax.legend(fontsize=9); ax.set_ylim(0, 1.05)
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(100))
        ax.set_xlim(left=0)
    fig.suptitle('Convergence Time CDF: run_002 vs run_003 (TF-restricted)',
                 fontsize=12, fontweight='bold')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'run002_vs_run003_convergence.png'), dpi=150)
    plt.close(fig)



def plot_regime_shift_map(all_results, r2_data):
    """Confusion matrix: run_002 attractor → run_003 attractor, per combo."""
    if r2_data is None: return

    at_order = ['Type A','Type B','Type C','Type D','Type E','Type F','Type G']
    combos = [make_combo_key(sk, srk) for sk in STRENGTH_REGIMES for srk in SIGN_RATIOS]

    # pair worlds by combo+topo_idx
    r3_map = {}
    for r in all_results:
        key = (r['combo_key'], r['topo_idx'])
        r3_map[key] = r['primary_attractor']
    r2_map = {}
    for r in r2_data:
        ck = r.get('combo_key', ''); ti = r.get('topo_idx')
        at = r.get('primary_attractor', '?')
        if ck and ti is not None:
            r2_map[(ck, ti)] = at

    # Build confusion matrix across all combos
    cm = {a2: {a3: 0 for a3 in at_order} for a2 in at_order}
    for key, a2 in r2_map.items():
        a3 = r3_map.get(key, '?')
        if a2 in cm and a3 in cm[a2]:
            cm[a2][a3] += 1

    # Check if any data exists
    total = sum(sum(row.values()) for row in cm.values())
    if total == 0: return

    fig, ax = plt.subplots(figsize=(7, 6))
    n = len(at_order)
    matrix = np.zeros((n, n))
    for i, a2 in enumerate(at_order):
        for j, a3 in enumerate(at_order):
            matrix[i, j] = cm[a2][a3]

    im = ax.imshow(matrix, cmap='Blues', aspect='auto')
    plt.colorbar(im, ax=ax, label='world count')
    ax.set_xticks(range(n)); ax.set_yticks(range(n))
    ax.set_xticklabels(at_order, fontsize=9)
    ax.set_yticklabels(at_order, fontsize=9)
    ax.set_xlabel('run_003 (TF-restricted)')
    ax.set_ylabel('run_002 (unrestricted)')
    for i in range(n):
        for j in range(n):
            if matrix[i, j] > 0:
                ax.text(j, i, str(int(matrix[i, j])), ha='center', va='center',
                        fontsize=10, color='white' if matrix[i,j] > total*0.1 else 'black',
                        fontweight='bold')
    ax.set_title(f'Regime Shift: run_002 → run_003\n({total} paired worlds)')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, 'regime_shift_map.png'), dpi=150)
    plt.close(fig)


# ── perturbation analysis ──

def pick_perturbation_targets(all_results):
    """Pick one world per attractor type for perturbation analysis."""
    global PERTURBATION_SELECTION
    selected = {}
    for at in ATTRACTOR_NAMES:
        for r in all_results:
            if r['primary_attractor'] == at:
                wid = r['world_id']
                # pick cell 0
                selected[at] = (wid, at, 0)
                break
    PERTURBATION_SELECTION = list(selected.values())
    print(f'  perturbation targets: {[(w,a) for w,a,_ in PERTURBATION_SELECTION]}')


def run_perturbation_analysis():
    if not PERTURBATION_SELECTION:
        print('  skip perturbation: no targets selected'); return
    os.makedirs(PERTURBATIONS_DIR, exist_ok=True)
    summary = []

    for wid, original_type, cell_idx in PERTURBATION_SELECTION:
        cell_label = f'cell{cell_idx:02d}'
        print(f'\n  perturbation: {wid} / {cell_label} (original: {original_type})')

        meta_path = os.path.join(METADATA_DIR, f'{wid}.json')
        with open(meta_path) as f: meta = json.load(f)

        w = dma.World(meta['seed']); w.from_dict(meta)

        traj_path = os.path.join(TRAJECTORIES_DIR, f'{wid}_{cell_label}.pt')
        if not os.path.exists(traj_path):
            print(f'    skip: trajectory not found'); continue
        orig_traj = torch.load(traj_path, weights_only=False)
        X_at_pert = orig_traj['X_traj'][PERTURBATION_TIME]
        kd_gene = int(torch.argmax(X_at_pert).item())
        max_expr = float(X_at_pert[kd_gene].item())
        kd_gene_role = 'TF' if kd_gene in dma.TF_GENES else 'non-TF'
        print(f'    highest expr at t={PERTURBATION_TIME}: gene {kd_gene} ({kd_gene_role}) = {max_expr:.4f}')

        cell_seed = meta['seed'] + cell_idx + 1
        X0 = dma.sample_initial_state(cell_seed)
        result = dma.simulate_single_cell(
            w, X0, dma.T_SIM,
            intervention_time=PERTURBATION_TIME,
            intervention_config={'knockdown': [kd_gene]},
        )
        X_traj_full = result['X_traj']
        clip_count_full = result['clip_count']

        # post-perturbation analysis
        X_post = X_traj_full[PERTURBATION_TIME:]
        clip_post = clip_count_full[PERTURBATION_TIME:]
        total_clips = int(clip_post.sum().item())

        world_dict = meta.copy(); world_dict['world_id'] = wid
        analysis = analyze_trajectory(X_post, clip_post, total_clips, world_dict)
        new_type = analysis['attractor_type']
        print(f'    attractor: {original_type} -> {new_type}')

        # recovery metrics
        eq_conv = analysis['equilibrium']['converged']
        eq_time = analysis['equilibrium']['convergence_time']
        recovery_time = None; recovery_failure = False
        if eq_conv:
            recovery_time = PERTURBATION_TIME + eq_time - PERTURBATION_TIME
        else:
            recovery_failure = True

        pert_history = {
            'time': PERTURBATION_TIME, 'type': 'knockdown',
            'knockdown_genes': [kd_gene],
        }
        torch.save({
            'X_traj': X_traj_full, 'clip_count': clip_count_full,
            'total_clips': total_clips, 'world_seed': meta['seed'],
            'cell_seed': cell_seed, 'world_id': wid,
            'original_attractor': original_type,
            'perturbed_attractor': new_type,
            'perturbation_history': pert_history,
            'knockdown_gene': kd_gene,
            'knockdown_gene_role': kd_gene_role,
            'knockdown_gene_expr_at_intervention': max_expr,
            'intervention_time': PERTURBATION_TIME,
            'perturbation_type': 'knockdown',
        }, os.path.join(PERTURBATIONS_DIR, f'{wid}_{cell_label}_perturb.pt'))

        summary.append({
            'world_id': wid, 'cell': cell_label,
            'original_attractor': original_type,
            'perturbed_attractor': new_type,
            'attractor_changed': new_type != original_type,
            'perturbation_history': pert_history,
            'knockdown_gene': kd_gene,
            'knockdown_gene_role': kd_gene_role,
            'knockdown_gene_expr_at_intervention': max_expr,
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

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    with open(os.path.join(ANALYSIS_DIR, 'perturbation_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    changed = sum(1 for s in summary if s['attractor_changed'])
    print(f'\n  perturbation done: {changed}/{len(summary)} attractor types changed')


def plot_perturbation_trajectories():
    pert_dir = os.path.join(FIGURES_DIR, 'perturbation')
    os.makedirs(pert_dir, exist_ok=True)
    grouped = {}
    for wid, atype, cell_idx in PERTURBATION_SELECTION:
        grouped.setdefault(atype, []).append((wid, cell_idx))
    for atype, entries in grouped.items():
        at_name = ATTRACTOR_NAMES.get(atype, atype)
        n = len(entries)
        fig, axes = plt.subplots(n, 2, figsize=(14, 4*n))
        if n == 1: axes = axes.reshape(1, 2)
        for r, (wid, cell_idx) in enumerate(entries):
            cl = f'cell{cell_idx:02d}'
            # original
            tp_o = os.path.join(TRAJECTORIES_DIR, f'{wid}_{cl}.pt')
            if os.path.exists(tp_o):
                orig = torch.load(tp_o, weights_only=False)
                for g in range(dma.G):
                    if g in dma.TF_GENES:
                        axes[r,0].plot(orig['X_traj'].numpy()[:,g],
                                       color='red', alpha=0.8, linewidth=1.0,
                                       label='TF' if g == 0 else '')
                    else:
                        axes[r,0].plot(orig['X_traj'].numpy()[:,g],
                                       color='grey', alpha=0.35, linewidth=0.8,
                                       label='non-TF' if g == dma.N_TF else '')
                axes[r,0].legend(loc='upper right', fontsize=9)
                axes[r,0].axvline(x=PERTURBATION_TIME, color='gray', linestyle='--',
                                  linewidth=0.8, alpha=0.6)
                axes[r,0].set_title(f'{wid} — before', fontsize=10)
            # perturbed
            tp_p = os.path.join(PERTURBATIONS_DIR, f'{wid}_{cl}_perturb.pt')
            if os.path.exists(tp_p):
                pert = torch.load(tp_p, weights_only=False)
                for g in range(dma.G):
                    if g in dma.TF_GENES:
                        axes[r,1].plot(pert['X_traj'].numpy()[:,g],
                                       color='red', alpha=0.8, linewidth=1.0,
                                       label='TF' if g == 0 else '')
                    else:
                        axes[r,1].plot(pert['X_traj'].numpy()[:,g],
                                       color='grey', alpha=0.35, linewidth=0.8,
                                       label='non-TF' if g == dma.N_TF else '')
                axes[r,1].legend(loc='upper right', fontsize=9)
                axes[r,1].axvline(x=PERTURBATION_TIME, color='red', linestyle='--',
                                  linewidth=0.8, alpha=0.6)
                axes[r,1].set_title(f'{wid} [{pert.get("perturbed_attractor","?")}] — after KD',
                                    fontsize=10)
            for ax in axes[r]: ax.tick_params(labelsize=7)
        fig.suptitle(f'Perturbation — {atype}: {at_name} [TF-restricted]',
                     fontsize=13, fontweight='bold')
        fig.tight_layout()
        fname = f'perturbation_{atype.replace(" ", "_").lower()}.png'
        fig.savefig(os.path.join(pert_dir, fname), dpi=150)
        plt.close(fig)
        print(f'  perturbation figure saved: {fname}')


# ── JSON output ──

def write_analysis_json(all_results, aggregation):
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    with open(os.path.join(ANALYSIS_DIR, 'aggregated_analysis.json'), 'w') as f:
        json.dump(aggregation, f, indent=2, default=str)
    serializable = []
    for r in all_results:
        item = {k: (float(v) if isinstance(v, (np.floating,)) else v)
                for k, v in r.items() if k != 'per_cell'}
        serializable.append(item)
    with open(os.path.join(ANALYSIS_DIR, 'all_world_results.json'), 'w') as f:
        json.dump(serializable, f, indent=2, default=str)


# ── main ──

def main():
    os.makedirs(FIGURES_DIR, exist_ok=True)
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    os.makedirs(SUMMARY_DIR, exist_ok=True)
    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
    os.makedirs(METADATA_DIR, exist_ok=True)

    all_results: List[Dict[str, Any]] = []

    for sign_key, sign_ratio in SIGN_RATIOS.items():
        print(f"\n=== sign_ratio: {sign_key} ({sign_ratio}) ===")
        results = simulate_topo_stratified(sign_key, sign_ratio, save_trajectories=True)
        all_results.extend(results)

    print(f"\nTotal worlds: {len(all_results)}")

    aggregation = aggregate_across_combos(all_results)
    for ck, agg in sorted(aggregation.items()):
        print(f"  {ck}: {agg['attractor_distribution']}")

    # Tables
    rows = build_world_summary_table(all_results)
    write_world_summary_csv(rows)
    write_regime_summary_csv(aggregation, all_results)
    write_analysis_json(all_results, aggregation)

    # Base figures (aligned with run_002)
    print("\nGenerating base figures...")
    plot_canonical_trajectories(all_results)
    plot_spectral_radius(all_results)
    plot_convergence_distribution(all_results)
    plot_clipping_distribution(all_results)
    plot_strength_vs_regime_type(all_results)
    plot_repression_ratio_vs_stability(all_results)
    plot_clipping_dominated_by_regime(all_results)
    plot_degree_distribution_vs_stability(all_results)
    plot_topology_exemplars(all_results)
    plot_topo_strength_interaction(all_results)
    plot_all_worlds_trajectories()

    # Perturbation
    pick_perturbation_targets(all_results)
    run_perturbation_analysis()
    plot_perturbation_trajectories()

    # run_002 vs run_003 comparison figures
    r2_data = load_run002_results()
    if r2_data is not None:
        print("\nGenerating run_002 vs run_003 comparison figures...")
        plot_run002_vs_run003_attractor(all_results, r2_data)
        plot_run002_vs_run003_spectral(all_results, r2_data)
        plot_run002_vs_run003_clipping(all_results, r2_data)
        plot_run002_vs_run003_convergence(all_results, r2_data)
        plot_regime_shift_map(all_results, r2_data)
    else:
        print(f"\n  skip run_002 comparisons: no data at {RUN_002_DIR}")

    print(f"\nDone. Output in: {RUN_DIR}")
    print(f"  figures: {FIGURES_DIR}")
    print(f"  tables:  {TABLES_DIR}")
    print(f"  analysis: {ANALYSIS_DIR}")


if __name__ == '__main__':
    main()