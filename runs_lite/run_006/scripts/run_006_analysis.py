#!/usr/bin/env python3
"""
Run 006 — Model A Edge Perturbation Sensitivity Experiment
=============================================================

Execution grid: core (12 selected worlds from run_005)
  Fork parallel continuation: baseline branch vs perturbed branch

Core question:
  Does acute regulatory edge deletion reveal hidden sensitivity in Model A worlds?

Pipeline:
  §1 — World selection & data loading (edge table, SCC, spectral info)
  §2 — Initial-condition screening & pre-perturb state classes
  §3 — Edge deletion strategy generation
  §4 — Fork simulation (baseline + perturbed branches)
  §5 — Trajectory-level regime classification
  §6 — Perturbation outcome assignment
  §7 — State landing analysis (within-world clustering)
  §8 — Table generation
  §9 — Figure generation

See: docs/01_DDC_Lite_Curriculum/run_006/  (00_Problem_Definition, 01_Analysis_Plan, 02_Output_Spec)
"""
import json
import os
import sys
import csv as _csv
import glob as _glob
import re
import time
import numpy as np
import pandas as pd
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from sklearn.decomposition import PCA
from collections import Counter, defaultdict
from typing import Any, Dict, List, Optional, Tuple, Set

# ── Add ddc source to path ────────────────────────────────────────
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))
import ddc_model_a as dma
import util

# ── Operating-condition override (aligned with run_005) ───────────
dma.DEFAULT_B     = 10.0
dma.DEFAULT_DELTA = 0.2

# ── Project paths ─────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
PROJ_ROOT   = os.path.dirname(SCRIPT_DIR)                       # run_006/
RUN_005_ROOT = os.path.join(os.path.dirname(PROJ_ROOT), 'run_005')

DATA_DIR    = os.path.join(PROJ_ROOT, 'data')
TRAJ_DIR    = os.path.join(PROJ_ROOT, 'trajectories')
FIGURES_DIR = os.path.join(PROJ_ROOT, 'figures')
SUMMARY_DIR = os.path.join(PROJ_ROOT, 'summary')

for d in [DATA_DIR, TRAJ_DIR, FIGURES_DIR, SUMMARY_DIR]:
    os.makedirs(d, exist_ok=True)

RUN_005_META_DIR = os.path.join(RUN_005_ROOT, 'world_metadata')
RUN_005_TRAJ_DIR = os.path.join(RUN_005_ROOT, 'trajectories')
RUN_005_TABLES_DIR = os.path.join(RUN_005_ROOT, 'tables')


# ═══════════════════════════════════════════════════════════════════
# Constants
# ═══════════════════════════════════════════════════════════════════

EXECUTION_GRID = 'core'

G     = 20
B     = 10.0
DELTA = 0.2
W_DIAG = 1.0 - DELTA  # 0.8

T_FORK  = 500          # perturbation intervention time
T_POST  = 500          # post-fork simulation steps per branch
T_TOTAL = T_FORK + T_POST  # 1000 — full trajectory length
K_DELETE = 5           # number of edges to delete
R_RANDOM = 10          # repeats for Strategy A (random), per doc §8.1 recommendation
N_INIT_SCREEN = 5      # initial conditions per world for screening
LATE_WINDOW  = 100     # last N steps for late-window analysis (per Analysis Plan §13.2)

# ── 12 selected worlds (from 01_Analysis_Plan §4) ─────────────────
SELECTED_WORLDS = [
    # H1 — Clip-free stable anchors
    {'world_id': 'baseline_balanced_r010_t000',
     'tier': 'H1_stable',
     'rationale': 'lowest rho(W) in full dataset, zero clipping; gold-standard stable'},
    {'world_id': 'chen_balanced_r010_t006',
     'tier': 'H1_stable',
     'rationale': 'moderate strength still zero clipping; baseline vs chen robustness'},

    # H2 — Boundary: stress_balanced_r010 gradient
    {'world_id': 'stress_balanced_r010_t003',
     'tier': 'H2_boundary',
     'rationale': 'convergent-side anchor; clip-mediated stable'},
    {'world_id': 'stress_balanced_r010_t001',
     'tier': 'H2_boundary',
     'rationale': 'mixed_regime specimen; deletion may push toward convergent or non-convergent'},
    {'world_id': 'stress_balanced_r010_t007',
     'tier': 'H2_boundary',
     'rationale': 'pure sustained oscillatory; test whether deletion breaks oscillation'},
    {'world_id': 'stress_balanced_r010_t005',
     'tier': 'H2_boundary',
     'rationale': 'pure divergent; test whether deletion rescues'},

    # H2 — Other boundary combo controls
    {'world_id': 'chen_balanced_r025_t005',
     'tier': 'H2_control',
     'rationale': 'rare convergent survivor in same combo'},
    {'world_id': 'chen_balanced_r025_t000',
     'tier': 'H2_control',
     'rationale': 'pure oscillatory; cross-combo comparison with #5'},
    {'world_id': 'stress_repression_biased_r010_t004',
     'tier': 'H2_control',
     'rationale': 'rho(W) just above 1 but pure oscillation; test rho(W) insufficiency'},

    # H3 — Rare specimens and regime outliers
    {'world_id': 'chen_repression_biased_r025_t007',
     'tier': 'H3_rare',
     'rationale': 'run_005 only multi-equilibrium world'},
    {'world_id': 'chen_repression_biased_r025_t001',
     'tier': 'H3_rare',
     'rationale': 'same combo, rho(W) near #10; topology contrast control'},
    {'world_id': 'stress_balanced_r025_t001',
     'tier': 'H3_rare',
     'rationale': 'high-rho divergent control; test whether deletion can rescue'},
]

SELECTED_IDS = [w['world_id'] for w in SELECTED_WORLDS]
SELECTED_LOOKUP = {w['world_id']: w for w in SELECTED_WORLDS}

# Init-type selection (4 types × 1 each + 1 extra)
INIT_TYPES = ['low', 'medium', 'high', 'sparse']
N_PER_TYPE = 1  # pick 1 per type
N_EXTRA    = 1  # pick 1 extra from 'low'


# ═══════════════════════════════════════════════════════════════════
# §1  World selection & data loading
# ═══════════════════════════════════════════════════════════════════

def s1_load_and_prepare():
    """Load 12 selected worlds from run_005, build matrices, compute topology.

    Returns
    -------
    selected_worlds : dict
        world_id → {
            'meta':          world_metadata JSON dict,
            'summary_row':   world_summary.tsv row (dict),
            'W':             (G, G) effective transition matrix W = I-Delta + A_reg (numpy),
            'b_vec':         (G, 1) bias vector (numpy),
            'rho_W':         float spectral radius of W,
            'scc_info':      dict from util.compute_scc_from_pgraph(),
            'n_edges':       int,
            'n_activation':  int,
            'n_repression':  int,
        }
    """
    print(f'\n{"="*70}')
    print(f'§1  World Selection & Data Loading')
    print(f'{"="*70}')

    # ── §1.1  Load world_summary.tsv and filter ──
    ws_path = os.path.join(RUN_005_TABLES_DIR, 'world_summary.tsv')
    df_ws = pd.read_csv(ws_path, sep='\t')
    df_ws.set_index('world_id', inplace=True)
    print(f'[§1.1] Loaded {len(df_ws)} worlds from world_summary.tsv')

    selected = {}
    for wid in SELECTED_IDS:
        if wid not in df_ws.index:
            print(f'  WARNING: {wid} not found in world_summary.tsv, skipping')
            continue
        selected[wid] = {
            'summary_row': df_ws.loc[wid].to_dict(),
        }

    # ── §1.2  Load world_metadata JSONs ──
    for wid in SELECTED_IDS:
        if wid not in selected:
            continue
        meta_path = os.path.join(RUN_005_META_DIR, f'{wid}.json')
        if not os.path.exists(meta_path):
            print(f'  WARNING: metadata JSON not found for {wid}, skipping')
            del selected[wid]
            continue
        with open(meta_path) as f:
            selected[wid]['meta'] = json.load(f)

    print(f'[§1.2] Loaded metadata for {len(selected)}/{len(SELECTED_IDS)} worlds')

    # ── §1.3  Build W matrices, compute rho(W) and SCC ──
    for world_idx, (wid, info) in enumerate(sorted(selected.items())):
        info['_world_idx'] = world_idx  # for stable seeding in §3
        meta = info['meta']
        W, b_vec = util.build_A_from_world_dict(meta, G=G)   # note: named 'A' in util but is W = I-Delta + A_reg
        info['W'] = W
        info['b_vec'] = b_vec

        # Spectral radius of W (already has self-decay 1-delta on diagonal)
        eigvals = np.linalg.eigvals(W)
        info['rho_W'] = float(np.max(np.abs(eigvals)))

        # SCC from P_graph
        pg = {int(k): v for k, v in meta['P_graph'].items()}
        info['scc_info'] = util.compute_scc_from_pgraph(pg, G=G)

        # Edge counts
        n_act = sum(1 for i in range(G)
                    for s in meta.get('edge_signs', {}).get(str(i), {}).values()
                    if s == dma.ACTIVATION)
        n_rep = sum(1 for i in range(G)
                    for s in meta.get('edge_signs', {}).get(str(i), {}).values()
                    if s == dma.REPRESSION)
        info['n_edges'] = n_act + n_rep
        info['n_activation'] = n_act
        info['n_repression'] = n_rep

        tier = SELECTED_LOOKUP[wid]['tier']
        rho = info['rho_W']
        scc_size = info['scc_info']['scc_sizes'][0] if info['scc_info']['scc_sizes'] else 0
        regime = info['summary_row'].get('world_attractor_label', '?')
        print(f'  [{tier}] {wid}: rho(W)={rho:.4f}, n_edges={info["n_edges"]}, '
              f'SCC={scc_size}, regime={regime}')

    # ── §1.4  Build edge table ──
    edge_rows = _build_edge_table(selected)
    print(f'[§1.4] edge_table_pre.csv: {len(edge_rows)} edges')

    # ── §1.5  Build source world selection table ──
    sel_rows = _build_source_world_selection(selected)
    print(f'[§1.5] source_world_selection.csv: {len(sel_rows)} worlds')

    return selected


# ── §1 helpers ────────────────────────────────────────────────────

def _build_edge_table(selected_worlds):
    """Build edge_table_pre.csv rows from selected worlds."""
    rows_written = []
    for wid in sorted(selected_worlds):
        info = selected_worlds[wid]
        meta = info['meta']
        scc_info = info['scc_info']
        comp_labels = scc_info['comp_labels']
        largest_id = scc_info['largest_scc_id']

        pg = {int(k): v for k, v in meta['P_graph'].items()}
        signs = {int(k): {int(sk): sv for sk, sv in v.items()}
                 for k, v in meta.get('edge_signs', {}).items()}
        strengths = {}
        es_dict = meta.get('edge_strengths', {})
        if es_dict:
            strengths = {int(k): {int(sk): sv for sk, sv in v.items()}
                         for k, v in es_dict.items()}

        # Pre-compute out-degree / in-degree
        out_deg = {i: 0 for i in range(G)}
        in_deg = {i: 0 for i in range(G)}
        for tgt, regs in pg.items():
            for src in regs:
                out_deg[src] += 1
                in_deg[tgt] += 1

        edge_idx = 0
        for tgt, regs in pg.items():
            for src in regs:
                sign = signs[tgt][src]
                strength_val = strengths.get(tgt, {}).get(src, 0.0)

                # is_in_largest_scc: both source and target in the same largest SCC
                in_largest = (comp_labels[src] == largest_id
                              and comp_labels[tgt] == largest_id)

                edge_id = f'{wid}_e{edge_idx:04d}'
                rows_written.append({
                    'source_world_id':    wid,
                    'edge_id':            edge_id,
                    'source_gene':        src,
                    'target_gene':        tgt,
                    'sign':               sign,
                    'strength':           round(sign * strength_val, 8),
                    'abs_strength':       round(abs(strength_val), 8),
                    'source_outdegree':   out_deg[src],
                    'target_indegree':    in_deg[tgt],
                    'is_in_largest_scc':  in_largest,
                    'is_cycle_associated': 'NA',
                })
                edge_idx += 1

    # Write CSV
    if rows_written:
        out_path = os.path.join(DATA_DIR, 'edge_table_pre.csv')
        with open(out_path, 'w', newline='') as f:
            w = _csv.DictWriter(f, fieldnames=rows_written[0].keys())
            w.writeheader()
            w.writerows(rows_written)

    return rows_written


def _build_source_world_selection(selected_worlds):
    """Build source_world_selection.csv rows."""
    columns = [
        'source_world_id', 'selection_tier', 'combo_key',
        'G', 'b', 'delta', 'r',
        'strength_regime', 'strength_low', 'strength_high', 'sign_ratio',
        'edge_count_pre', 'activation_count_pre', 'repression_count_pre',
        'rho_W_pre', 'largest_scc_size_pre',
        'source_primary_regime', 'source_secondary_label',
        'source_attractor_label', 'source_equilibrium_count',
        'source_clipping_frequency', 'selection_rationale',
    ]

    # ── Collect per-world secondary_labels from run_005 trajectory_summary ──
    traj_path = os.path.join(RUN_005_TABLES_DIR, 'trajectory_summary.tsv')
    secondary_by_world: Dict[str, Dict[str, int]] = {}  # wid -> {label: count}
    traj_total_by_world: Dict[str, int] = {}             # wid -> total traj count
    if os.path.exists(traj_path):
        with open(traj_path) as f:
            reader = _csv.DictReader(f, delimiter='\t')
            for row in reader:
                wid = row['world_id']
                sl = row.get('secondary_labels', '').strip()
                if sl:
                    sec = secondary_by_world.setdefault(wid, {})
                    sec[sl] = sec.get(sl, 0) + 1
                traj_total_by_world[wid] = traj_total_by_world.get(wid, 0) + 1

    rows = []
    for wid in sorted(selected_worlds):
        info = selected_worlds[wid]
        meta = info['meta']
        sr = info['summary_row']
        lookup = SELECTED_LOOKUP[wid]

        # Extract r from world_id
        r_val = 0.25 if 'r025' in wid else 0.10

        # combo_key: {strength}_{sign_ratio}
        strength_regime = meta.get('strength_regime', '')
        sign_label = meta.get('sign_ratio_label', '')
        combo_key = f'{strength_regime}_{sign_label}'

        scc_sizes = info['scc_info']['scc_sizes']
        largest_scc_size = scc_sizes[0] if scc_sizes else 0

        # source_secondary_label: all distinct secondary_labels across trajectories, with counts
        sl_counts = secondary_by_world.get(wid, {})
        total_traj = traj_total_by_world.get(wid, 0)
        if sl_counts:
            parts = [f'{label}({cnt}/{total_traj})'
                     for label, cnt in sorted(sl_counts.items(), key=lambda x: -x[1])]
            secondary_label = ' | '.join(parts)
        else:
            secondary_label = ''

        rows.append({
            'source_world_id':             wid,
            'selection_tier':              lookup['tier'],
            'combo_key':                   combo_key,
            'G':                           G,
            'b':                           B,
            'delta':                       DELTA,
            'r':                           r_val,
            'strength_regime':             strength_regime,
            'strength_low':                meta.get('a_min', ''),
            'strength_high':               meta.get('a_max', ''),
            'sign_ratio':                  sign_label,
            'edge_count_pre':              info['n_edges'],
            'activation_count_pre':        info['n_activation'],
            'repression_count_pre':        info['n_repression'],
            'rho_W_pre':                   round(info['rho_W'], 4),
            'largest_scc_size_pre':        largest_scc_size,
            'source_primary_regime':       sr.get('world_attractor_label', ''),
            'source_secondary_label':      secondary_label,
            'source_attractor_label':      sr.get('world_attractor_label', ''),
            'source_equilibrium_count':    sr.get('number_of_distinct_equilibria', ''),
            'source_clipping_frequency':   sr.get('mean_clipping_frequency', ''),
            'selection_rationale':         lookup['rationale'],
        })

    if rows:
        out_path = os.path.join(DATA_DIR, 'source_world_selection.csv')
        with open(out_path, 'w', newline='') as f:
            w = _csv.DictWriter(f, fieldnames=columns)
            w.writeheader()
            w.writerows(rows)

    return rows


# ═══════════════════════════════════════════════════════════════════
# §2  Initial-condition screening & pre-perturb state classes
# ═══════════════════════════════════════════════════════════════════

CLUSTER_THRESHOLD = 1e-3       # relative L2 threshold for X(t_fork) clustering (convergent worlds only)


def _pick_init_files(world_id: str) -> List[Dict[str, Any]]:
    """Pick up to N_INIT_SCREEN init files per world from run_005 trajectories.

    Strategy: take 1 per init_type (low, medium, high, sparse), then 1 extra.
    Returns list sorted by init_id.
    """
    pattern = os.path.join(RUN_005_TRAJ_DIR, f'{world_id}_init_*.pt')
    all_files = sorted(_glob.glob(pattern))
    if not all_files:
        return []

    # Parse init_type from filename: ..._init_<type>_<NN>.pt
    pat = re.compile(rf'{re.escape(world_id)}_init_(\w+)_(\d+)\.pt$')

    by_type: Dict[str, List[Tuple[str, str]]] = {}
    for fp in all_files:
        basename = os.path.basename(fp)
        m = pat.match(basename)
        if not m:
            continue
        itype, num = m.group(1), m.group(2)
        if itype not in INIT_TYPES:
            continue  # skip unknown types
        by_type.setdefault(itype, []).append((f'{itype}_{num}', fp))

    picked = []
    for itype in INIT_TYPES:
        entries = sorted(by_type.get(itype, []), key=lambda x: x[0])
        if entries:
            picked.append(entries[0])          # first in type

    # Extra: take second from 'low' (index 1) if available, otherwise from any type
    if len(picked) < N_INIT_SCREEN:
        for itype in INIT_TYPES:  # prefer 'low' first
            entries = sorted(by_type.get(itype, []), key=lambda x: x[0])
            if len(entries) > 1:
                extra = entries[1]
                if extra[0] not in {p[0] for p in picked}:
                    picked.append(extra)
                    break

    # Load and return
    results = []
    for init_id, fp in picked:
        data = torch.load(fp, map_location='cpu')
        # X_traj shape (T, G); t_fork is 500
        X_fork = data['X_traj'][T_FORK].numpy().astype(np.float64)
        X0 = data['X0'].numpy().astype(np.float64)
        results.append({
            'init_id':   init_id,
            'init_type': data.get('init_type', init_id.split('_')[0]),
            'X_fork':    X_fork,
            'X0':        X0,
            'pt_path':   fp,
        })
    return results


def _cluster_preperturb_states(
    init_files: List[Dict[str, Any]],
    world_id: str,
) -> Dict[str, Any]:
    """Cluster X(t_fork) states with relative-L2 average-linkage.

    Returns dict with keys: class_dict, labels, abs_dm, rel_dm, n_classes.
    """
    states = [np.array(if_['X_fork']) for if_ in init_files]
    labels, abs_dm, rel_dm = util.hierarchical_cluster_relative_l2(
        states, threshold=CLUSTER_THRESHOLD)

    # Group indices by label
    class_dict: Dict[int, List[int]] = {}
    for idx, lbl in enumerate(labels):
        class_dict.setdefault(int(lbl), []).append(idx)

    return {
        'class_dict': class_dict,
        'labels':     labels,
        'abs_dm':     abs_dm,
        'rel_dm':     rel_dm,
        'n_classes':  len(class_dict),
    }


def _select_representatives(
    init_files: List[Dict[str, Any]],
    class_dict: Dict[int, List[int]],
) -> List[Dict[str, Any]]:
    """Select representative trajectories per pre-perturb state class.

    For each class, pick the member closest to centroid.
    Used only for convergent worlds (non_convergent worlds skip clustering).
    """
    reps = []
    for cid, members in sorted(class_dict.items()):
        states = np.stack([init_files[m]['X_fork'] for m in members])
        centroid = np.mean(states, axis=0)
        dists = np.linalg.norm(states - centroid, axis=1)
        closest_idx = int(np.argmin(dists))
        closest_member = members[closest_idx]
        reps.append({
            'class_id':  str(cid),
            'init_id':   init_files[closest_member]['init_id'],
            'init_type': init_files[closest_member]['init_type'],
            'X_fork':    init_files[closest_member]['X_fork'],
            'X0':        init_files[closest_member]['X0'],
            'pt_path':   init_files[closest_member]['pt_path'],
            'selection': 'closest_to_centroid',
        })
    return reps


def _make_init_direct_reps(
    init_files: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """For non_convergent worlds: each init becomes its own representative class."""
    reps = []
    for idx, f in enumerate(init_files):
        reps.append({
            'class_id':  str(idx),
            'init_id':   f['init_id'],
            'init_type': f['init_type'],
            'X_fork':    f['X_fork'],
            'X0':        f['X0'],
            'pt_path':   f['pt_path'],
            'selection': 'direct_no_clustering',
        })
    return reps


def _save_preperturb_state_classes(
    screening: Dict[str, Any],
) -> None:
    """Write data/pre_perturb_state_classes.csv."""
    columns = [
        'source_world_id', 'pre_perturb_state_class_id', 't_perturb',
        'pre_primary_regime', 'pre_secondary_label',
        'n_screened_initial_conditions', 'member_initial_condition_ids',
        'representative_initial_condition_id', 'representative_trajectory_id',
        'representative_state_vector_path', 'state_class_assignment_method',
        'within_class_max_distance', 'notes',
    ]
    rows = []
    for wid, scr in screening.items():
        sr = scr['world_sr']
        regime = sr.get('world_attractor_label', '')
        cluster_info = scr['cluster_info']
        class_dict = cluster_info['class_dict']
        rel_dm = cluster_info['rel_dm']
        init_files = scr['init_files']

        # Build class_id → representative mapping
        rep_map = {r['class_id']: r for r in scr['representatives']}

        for cid, members in sorted(class_dict.items()):
            cid_str = str(cid)
            rep = rep_map.get(cid_str)
            member_ids = [init_files[m]['init_id'] for m in members]

            # Max pairwise rel-L2 within this class
            max_dist = 0.0
            for a in members:
                for b in members:
                    if a < b and rel_dm is not None:
                        d = rel_dm[a][b]
                        if d > max_dist:
                            max_dist = d

            is_convergent = scr.get('is_convergent', True)
            method = ('hierarchical_average_linkage_relative_L2'
                      if is_convergent else 'direct_per_init_no_clustering')
            notes = (f'clustering_threshold={CLUSTER_THRESHOLD}'
                     if is_convergent else 'non_convergent_world_each_init_direct')

            rows.append({
                'source_world_id':                   wid,
                'pre_perturb_state_class_id':        cid_str,
                't_perturb':                         T_FORK,
                'pre_primary_regime':                regime,
                'pre_secondary_label':               '',  # filled at §5
                'n_screened_initial_conditions':     len(init_files),
                'member_initial_condition_ids':      ','.join(member_ids),
                'representative_initial_condition_id': rep['init_id'] if rep else '',
                'representative_trajectory_id':       rep['init_id'] if rep else '',
                'representative_state_vector_path':    rep['pt_path'] if rep else '',
                'state_class_assignment_method':      method,
                'within_class_max_distance':          round(max_dist, 6),
                'notes': notes,
            })

    if rows:
        out_path = os.path.join(DATA_DIR, 'pre_perturb_state_classes.csv')
        with open(out_path, 'w', newline='') as f:
            w = _csv.DictWriter(f, fieldnames=columns)
            w.writeheader()
            w.writerows(rows)
        print(f'  Saved pre_perturb_state_classes.csv: {len(rows)} classes')


def _save_representative_xforks(
    screening: Dict[str, Any],
) -> None:
    """Write trajectories/{wid}_rep{n}_Xfork.pt for each representative."""
    count = 0
    for wid, scr in screening.items():
        for rep in scr['representatives']:
            out_path = os.path.join(TRAJ_DIR, f'{wid}_rep{rep["class_id"]}_Xfork.pt')
            torch.save({
                'world_id':  wid,
                'class_id':  rep['class_id'],
                'init_id':   rep['init_id'],
                'init_type': rep['init_type'],
                'X_fork':    torch.from_numpy(rep['X_fork']),
                'X0':        torch.from_numpy(rep['X0']),
                't_fork':    T_FORK,
                'selection': rep['selection'],
            }, out_path)
            count += 1
    print(f'  Saved {count} representative Xfork .pt files')


def s2_screen_init_conditions(selected_worlds):
    """Screen 5 initial conditions per world, cluster X(t_fork) states,
    select representative trajectories.

    Returns
    -------
    screening : dict
        world_id → {
            'init_files':       list of dicts (init_id, init_type, X_fork, ...),
            'cluster_info':     dict from _cluster_preperturb_states(),
            'representatives':  list of rep info dicts,
            'world_sr':         world_summary row for this world,
        }
    """
    print(f'\n{"="*70}')
    print(f'§2  Initial-condition Screening & Pre-perturb State Classes')
    print(f'{"="*70}')
    print(f'  N_screen = {N_INIT_SCREEN} per world, '
          f'clustering threshold = {CLUSTER_THRESHOLD}')

    screening = {}

    for wid, info in sorted(selected_worlds.items()):
        world_sr = info['summary_row']
        regime = world_sr.get('world_attractor_label', '')
        is_convergent = ('convergent' in regime.lower() 
                         and 'non_convergent' not in regime.lower())

        # §2.1 — Pick init files
        init_files = _pick_init_files(wid)
        if len(init_files) == 0:
            print(f'  [{wid}]  WARNING: no init files found, skipping')
            screening[wid] = {
                'init_files': [], 'cluster_info': {}, 'representatives': [],
                'world_sr': world_sr, 'is_convergent': is_convergent,
            }
            continue

        if is_convergent:
            # Convergent worlds: cluster X(t_fork), find multi-basin
            cluster_info = _cluster_preperturb_states(init_files, wid)
            reps = _select_representatives(init_files, cluster_info['class_dict'])
        else:
            # Non-convergent worlds: each init is its own class
            n = len(init_files)
            class_dict = {i: [i] for i in range(n)}
            cluster_info = {
                'class_dict': class_dict,
                'labels':     list(range(n)),
                'abs_dm':     None,
                'rel_dm':     None,
                'n_classes':  n,
            }
            reps = _make_init_direct_reps(init_files)

        screening[wid] = {
            'init_files':      init_files,
            'cluster_info':    cluster_info,
            'representatives': reps,
            'world_sr':        world_sr,
            'is_convergent':   is_convergent,
        }

        n_classes = cluster_info['n_classes']
        n_reps = len(reps)
        classes_str = ','.join(f'c{c}={len(v)}' for c, v in sorted(cluster_info['class_dict'].items()))
        print(f'  [{wid}]  regime={regime:40s}  '
              f'init_files={len(init_files)}, classes={n_classes}({classes_str}), reps={n_reps}')

    # §2.4 — Write pre_perturb_state_classes.csv
    _save_preperturb_state_classes(screening)

    # §2.5 — Write representative Xfork .pt files
    _save_representative_xforks(screening)

    total_reps = sum(len(scr['representatives']) for scr in screening.values())
    print(f'\n[§2]  DONE: {total_reps} total representatives across {len(screening)} worlds')
    return screening


# ═══════════════════════════════════════════════════════════════════
# §3  Edge deletion strategy generation
# ═══════════════════════════════════════════════════════════════════

STRATEGY_CODES = {
    'random_delete_5':                'A',
    'target_strongest_5':             'B',
    'target_hub_outgoing_5':          'C',
    'target_largest_scc_internal_5':  'D',
}


def _build_edge_catalog(wid: str, meta: dict, info: dict) -> List[dict]:
    """Build per-world edge catalog from metadata.

    Returns list of edge dicts with keys matching edge_table_pre.csv columns.
    """
    pg = {int(k): v for k, v in meta['P_graph'].items()}
    signs = {int(k): {int(sk): sv for sk, sv in v.items()}
             for k, v in meta.get('edge_signs', {}).items()}
    es_dict = meta.get('edge_strengths', {})
    strengths = {}
    if es_dict:
        strengths = {int(k): {int(sk): sv for sk, sv in v.items()}
                     for k, v in es_dict.items()}

    comp_labels = info['scc_info']['comp_labels']
    largest_id = info['scc_info']['largest_scc_id']

    # Read degrees from metadata (already computed in run_005)
    in_deg = {int(k): v for k, v in meta.get('in_degrees', {}).items()}
    out_deg = {int(k): v for k, v in meta.get('out_degrees', {}).items()}

    catalog = []
    edge_idx = 0
    for tgt in range(G):
        for src in pg.get(tgt, []):
            sign = signs[tgt][src]
            strength_val = strengths.get(tgt, {}).get(src, 0.0)
            in_largest = (comp_labels[src] == largest_id and
                          comp_labels[tgt] == largest_id)
            # comp_id: which SCC this edge is internal to, or -1 if inter-SCC
            comp_id = (comp_labels[src]
                       if comp_labels[src] == comp_labels[tgt] else -1)
            catalog.append({
                'edge_id':            f'{wid}_e{edge_idx:04d}',
                'source_world_id':    wid,
                'source_gene':        src,
                'target_gene':        tgt,
                'sign':               sign,
                'strength':           round(sign * strength_val, 8),  # signed
                'abs_strength':       round(abs(strength_val), 8),
                'source_outdegree':   out_deg[src],
                'target_indegree':    in_deg[tgt],
                'is_in_largest_scc':  in_largest,
                'comp_id':            comp_id,
                'is_cycle_associated': 'NA',
            })
            edge_idx += 1
    return catalog


def _compute_W_post(W: np.ndarray, deleted_sources: List[int],
                    deleted_targets: List[int], meta: Dict[str, Any] = None,
                    ) -> np.ndarray:
    """Return W' with deleted edges zeroed out.

    Two methods are compared for safety:
      method_1 — zero out off-diagonal entries in the original W
      method_2 — rebuild W from modified metadata (ground truth)
    If both are available they must agree within float tolerance;
    the ground-truth rebuild is returned.
    """
    # method_1: zero entries
    W_post = W.copy()
    for src, tgt in zip(deleted_sources, deleted_targets):
        if src != tgt:  # never zero out diagonal (self-decay)
            W_post[tgt, src] = 0.0

    # method_2: rebuild from metadata (ground truth)
    if meta is not None:
        meta_mod = json.loads(json.dumps(meta))  # deep-copy via json round-trip
        pg = {int(k): v for k, v in meta_mod['P_graph'].items()}
        # pg[tgt] shares the underlying list with meta_mod['P_graph'][str(tgt)],
        # so removing src here mutates meta_mod in-place — no need to write back.
        for src, tgt in zip(deleted_sources, deleted_targets):
            if src != tgt and tgt in pg and src in pg[tgt]:
                pg[tgt].remove(src)
        W_rebuilt, _ = util.build_A_from_world_dict(meta_mod, G=G)
        # Verify consistency with zero-out method
        if not np.allclose(W_post, W_rebuilt, atol=1e-12):
            raise ValueError(
                f'Mismatch between zero-out and rebuild methods: '
                f'max|diff| = {np.max(np.abs(W_post - W_rebuilt))}'
            )
        return W_rebuilt
    return W_post


def _strategy_a_random(catalog: List[dict], world_id: str, rep_class: str,
                       rep_seed: int, n_repeats: int) -> List[dict]:
    """Strategy A: uniform random K=5 edges across all edges.

    Each repeat gets its own independently-seeded rng so that
    random_seed_perturbation = rep_seed + r exactly reproduces that repeat.
    """
    plans = []
    for r in range(n_repeats):
        rng = np.random.default_rng(rep_seed + r)
        chosen = rng.choice(len(catalog), size=min(K_DELETE, len(catalog)),
                            replace=False)
        selected = [catalog[i] for i in chosen]
        plans.append(_make_plan(world_id, 'random_delete_5', rep_class, r,
                                  selected, rep_seed + r))
    return plans


def _strategy_b_strongest(catalog: List[dict], world_id: str, rep_class: str) -> List[dict]:
    """Strategy B: top K=5 edges by |strength| descending (deterministic)."""
    sorted_cat = sorted(catalog, key=lambda e: e['abs_strength'], reverse=True)
    selected = sorted_cat[:K_DELETE]
    return [_make_plan(world_id, 'target_strongest_5', rep_class, 0, selected, None)]


def _strategy_c_hub_outgoing(catalog: List[dict], world_id: str, rep_class: str,
                             rng: np.random.Generator, rep_seed: int) -> List[dict]:
    """Strategy C: outgoing edges of top hub genes, weighted sample by |strength|.

    Rank genes by out-degree descending. Expand hub gene set greedily until
    we have >= K edges; then weighted-sample K by abs_strength.
    """
    out_deg = {}
    for e in catalog:
        out_deg.setdefault(e['source_gene'], e['source_outdegree'])
    ranked_genes = sorted(out_deg.items(), key=lambda x: x[1], reverse=True)

    hub_edges = []
    hub_weights = []
    for gene, _ in ranked_genes:
        for e in catalog:
            if e['source_gene'] == gene:
                hub_edges.append(e)
                hub_weights.append(e['abs_strength'])
        if len(hub_edges) >= K_DELETE:
            break

    weights = np.array(hub_weights, dtype=np.float64)
    weights /= weights.sum()
    chosen = rng.choice(len(hub_edges), size=min(K_DELETE, len(hub_edges)),
                        replace=False, p=weights)
    selected = [hub_edges[i] for i in chosen]
    return [_make_plan(world_id, 'target_hub_outgoing_5', rep_class, 0, selected,
                       rep_seed)]


def _strategy_d_scc_internal(catalog: List[dict], world_id: str, rep_class: str,
                             scc_info: dict, rng: np.random.Generator,
                             rep_seed: int) -> List[dict]:
    """Strategy D: edges inside largest SCC. If < K, expand to next SCC.

    Picks edges internal to the largest SCC first, then expands through
    progressively smaller SCCs (by component size descending) until K edges
    are accumulated, per doc §7.3 Strategy D.
    """
    # Group edges by comp_id (only internal-to-same-SCC edges have comp_id >= 0)
    edges_by_comp = {}
    for e in catalog:
        cid = e['comp_id']
        if cid >= 0:
            edges_by_comp.setdefault(cid, []).append(e)

    # Component IDs sorted by SCC size (descending) = comp_id itself
    comp_order = sorted(edges_by_comp.keys(),
                        key=lambda c: scc_info['scc_sizes'][c],
                        reverse=True)

    selected = []
    for cid in comp_order:
        selected.extend(edges_by_comp[cid])
        if len(selected) >= K_DELETE:
            break

    if len(selected) > K_DELETE:
        chosen = rng.choice(len(selected), size=K_DELETE, replace=False)
        selected = [selected[i] for i in chosen]

    return [_make_plan(world_id, 'target_largest_scc_internal_5', rep_class, 0,
                       selected, rep_seed)]


def _make_plan(world_id: str, mode: str, rep_class: str, repeat_index: int,
               selected: List[dict], seed: Optional[int]) -> dict:
    """Assemble a single perturbation plan dict.

    seed: the rng seed that produced this plan's edge selection,
          or None if the strategy is deterministic.
    """
    pid = f'{world_id}_{STRATEGY_CODES[mode]}_r{repeat_index}_rep{rep_class}'
    return {
        'perturbation_id':                   pid,
        'source_world_id':                   world_id,
        'perturbation_mode':                 mode,
        'delete_k':                          len(selected),
        'repeat_index':                      repeat_index,
        'representative_class_id':           rep_class,
        'deleted_edge_ids':                  ','.join(e['edge_id'] for e in selected),
        'deleted_edge_sources':              ','.join(str(e['source_gene']) for e in selected),
        'deleted_edge_targets':              ','.join(str(e['target_gene']) for e in selected),
        'deleted_edge_signs':                ','.join(str(e['sign']) for e in selected),
        'deleted_edge_strengths':            ','.join(str(e['strength']) for e in selected),
        'deleted_abs_strength_sum':          round(sum(e['abs_strength'] for e in selected), 8),
        'deleted_activation_count':          sum(1 for e in selected if e['sign'] == dma.ACTIVATION),
        'deleted_repression_count':          sum(1 for e in selected if e['sign'] == dma.REPRESSION),
        'deleted_edges_in_largest_scc_count': sum(1 for e in selected if e['is_in_largest_scc']),
        'deleted_cycle_associated_count':    0,  # NA
        'random_seed_perturbation':          seed if seed is not None else '',
    }


def s3_generate_perturbations(selected_worlds, screening):
    """Generate edge deletion plans for all strategies (A–D).

    Returns
    -------
    pert_plan : list of perturbation plan dicts
    """
    print(f'\n{"="*70}')
    print(f'§3  Edge Deletion Strategy Generation')
    print(f'{"="*70}')
    print(f'  R_random = {R_RANDOM}, K_delete = {K_DELETE}')

    all_plans = []
    total_reps = sum(len(s['representatives']) for s in screening.values())

    for wid, info in sorted(selected_worlds.items()):
        meta = info['meta']
        scr = screening[wid]
        reps = scr['representatives']
        if not reps:
            continue

        # Build edge catalog once per world
        catalog = _build_edge_catalog(wid, meta, info)
        W = info['W']
        base_seed = info.get('_world_idx', 0) * 10000  # stable seeding

        for rep in reps:
            class_id = rep['class_id']
            # Seed per (world, rep) pair for reproducibility
            rep_seed = base_seed + int(class_id) * 1000
            rng = np.random.default_rng(rep_seed)

            # Strategy A — random (each repeat independently seeded)
            plans_a = _strategy_a_random(catalog, wid, class_id, rep_seed, R_RANDOM)
            all_plans.extend(plans_a)

            # Strategy B — strongest (deterministic, no seed)
            plans_b = _strategy_b_strongest(catalog, wid, class_id)
            all_plans.extend(plans_b)

            # Strategy C — hub outgoing
            plans_c = _strategy_c_hub_outgoing(catalog, wid, class_id, rng,
                                               rep_seed)
            all_plans.extend(plans_c)

            # Strategy D — largest-SCC internal
            plans_d = _strategy_d_scc_internal(catalog, wid, class_id,
                                               info['scc_info'], rng, rep_seed)
            all_plans.extend(plans_d)

    # Compute rho_W_pre / rho_W_post for each plan
    for plan in all_plans:
        wid = plan['source_world_id']
        W = selected_worlds[wid]['W']
        src_str = plan['deleted_edge_sources']
        tgt_str = plan['deleted_edge_targets']
        deleted_src = [int(x) for x in src_str.split(',')] if src_str else []
        deleted_tgt = [int(x) for x in tgt_str.split(',')] if tgt_str else []
        W_post = _compute_W_post(W, deleted_src, deleted_tgt,
                                meta=selected_worlds[wid]['meta'])
        plan['rho_W_pre'] = round(selected_worlds[wid]['rho_W'], 6)
        plan['rho_W_post'] = round(np.max(np.abs(np.linalg.eigvals(W_post))), 6)
        plan['delta_rho_W'] = round(plan['rho_W_post'] - plan['rho_W_pre'], 6)
        plan['n_edges_pre'] = selected_worlds[wid]['n_edges']
        plan['n_edges_post'] = plan['n_edges_pre'] - len(deleted_src)
        plan['t_perturb'] = T_FORK
        plan['T_post'] = T_POST
        plan['T_total'] = T_FORK + T_POST
        plan['representative_init_id'] = next(
            (r['init_id'] for r in screening[wid]['representatives']
             if r['class_id'] == plan['representative_class_id']), '')

    # Save perturbation_metadata.csv
    if all_plans:
        out_path = os.path.join(DATA_DIR, 'perturbation_metadata.csv')
        columns = [
            'perturbation_id', 'source_world_id', 'perturbation_mode',
            'delete_k', 'repeat_index',
            't_perturb', 'T_post', 'T_total',
            'representative_class_id', 'representative_init_id',
            'deleted_edge_ids', 'deleted_edge_sources', 'deleted_edge_targets',
            'deleted_edge_signs', 'deleted_edge_strengths',
            'deleted_abs_strength_sum',
            'deleted_activation_count', 'deleted_repression_count',
            'deleted_edges_in_largest_scc_count', 'deleted_cycle_associated_count',
            'rho_W_pre', 'rho_W_post', 'delta_rho_W',
            'n_edges_pre', 'n_edges_post',
            'random_seed_perturbation',
        ]
        with open(out_path, 'w', newline='') as f:
            w = _csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
            w.writeheader()
            w.writerows(all_plans)

    print(f'  Generated {len(all_plans)} perturbation plans '
          f'({total_reps} reps × {R_RANDOM + 3} strategies)')
    print(f'  Saved perturbation_metadata.csv')

    return all_plans


# ═══════════════════════════════════════════════════════════════════
# §4  Fork simulation (baseline + perturbed)
# ═══════════════════════════════════════════════════════════════════
#
# Scheme:
#   baseline  — full 1000-step sim from X0 with original W
#   perturbed — 500-step pre-perturb (W) + 500-step post-perturb (W_post)
#               stitched into a single 1001-point trajectory
#   sanity    — compare baseline full trajectory against run_005 X_traj
#
# Late window: last LATE_WINDOW steps of the full trajectory → [T_TOTAL - LATE_WINDOW, T_TOTAL]

def _run_branch(W: np.ndarray, b_vec: np.ndarray, X0: np.ndarray,
                n_steps: int) -> Tuple[np.ndarray, np.ndarray]:
    """Run simulation: X(t+1) = max(W @ X(t) + b, 0).

    Returns
    -------
    X_traj    : (n_steps+1, G)   trajectory including X0 at index 0
    clip_mask : (n_steps+1, G)   bool, which genes were clipped per step
    """
    G = len(X0)
    X_traj = np.zeros((n_steps + 1, G), dtype=np.float64)
    clip_mask = np.zeros((n_steps + 1, G), dtype=bool)
    X_traj[0] = X0
    X = X0.copy()
    for t in range(n_steps):
        raw = W @ X + b_vec.flatten()
        X_next = np.maximum(raw, 0.0)
        clip_mask[t + 1] = raw < 0
        X_traj[t + 1] = X_next
        X = X_next
    return X_traj, clip_mask


def _stitch_perturbed_traj(
    pre_traj: np.ndarray,   # (T_FORK+1, G)  — indices 0..T_FORK
    pre_mask: np.ndarray,   # (T_FORK+1, G)
    post_traj: np.ndarray,  # (T_POST+1, G)  — indices 0..T_POST
    post_mask: np.ndarray,  # (T_POST+1, G)
) -> Tuple[np.ndarray, np.ndarray]:
    """Stitch pre- and post-perturb segments into one full trajectory."""
    full_traj = np.concatenate([pre_traj, post_traj[1:]], axis=0)    # (1001, G)
    full_mask = np.concatenate([pre_mask, post_mask[1:]], axis=0)    # (1001, G)
    return full_traj, full_mask


def s4_run_fork_simulations(selected_worlds, screening, pert_plan):
    """Run fork simulations with baseline sanity check.

    For each perturbation plan:
      baseline  — simulates full T_TOTAL=1000 steps from X0 with W
                  and compares against run_005 source trajectory.
      perturbed — stitched two-segment:
                    t=0..500 with W, then t=500..1000 with W_post

    Returns
    -------
    fork_results : dict
        perturbation_id → {
            perturbation_id, source_world_id, X0, X_fork,
            baseline_traj, perturbed_traj,
            baseline_clip_mask, perturbed_clip_mask,
            late_window_start, late_window_end,
            sanity_passed,
        }
    """
    print(f'\n{"="*70}')
    print(f'§4  Fork Simulation (full-trajectory + baseline sanity)')
    print(f'{"="*70}')
    print(f'  T_total = {T_TOTAL}, T_fork = {T_FORK}, T_post = {T_POST}')
    print(f'  Late window = last {LATE_WINDOW} steps (indices [{T_TOTAL - LATE_WINDOW}, {T_TOTAL}])')

    # Build lookup: (wid, rep_class) → (X0, pt_path)
    rep_lookup: Dict[tuple, Tuple[np.ndarray, str]] = {}
    for wid, scr in screening.items():
        for rep in scr['representatives']:
            rep_lookup[(wid, rep['class_id'])] = (rep['X0'], rep['pt_path'])

    fork_results: Dict[str, dict] = {}
    total = len(pert_plan)
    sanity_fails = 0
    T_SIM = T_FORK + T_POST  # 1000

    for idx, plan in enumerate(pert_plan):
        wid = plan['source_world_id']
        rep_class = plan['representative_class_id']

        # Locate X0 and source trajectory path
        rep_info = rep_lookup.get((wid, rep_class))
        if rep_info is None:
            print(f'  WARNING: no rep info for {wid} rep{rep_class}, skipping')
            continue
        X0, pt_path = rep_info

        W = selected_worlds[wid]['W']
        b_vec = selected_worlds[wid]['b_vec']
        meta = selected_worlds[wid]['meta']

        # Parse deleted edges
        src_str = plan['deleted_edge_sources']
        tgt_str = plan['deleted_edge_targets']
        deleted_src = [int(x) for x in src_str.split(',')] if src_str else []
        deleted_tgt = [int(x) for x in tgt_str.split(',')] if tgt_str else []

        # Compute W_post (ground-truth rebuild with consistency check)
        W_post = _compute_W_post(W, deleted_src, deleted_tgt, meta=meta)

        # ── Baseline: full 1000-step simulation ──
        base_traj, base_mask = _run_branch(W, b_vec, X0, T_SIM)
        X_fork_val = base_traj[T_FORK]  # state at t=500

        # ── Sanity check: compare baseline against run_005 source ──
        if not pt_path or not os.path.exists(pt_path):
            raise FileNotFoundError(
                f'Source trajectory not found: {pt_path} (wid={wid}, rep={rep_class})')
        source_data = torch.load(pt_path, map_location='cpu')
        source_traj = source_data['X_traj'].numpy().astype(np.float64)
        if source_traj.shape != (T_SIM + 1, G):
            raise ValueError(
                f'Source trajectory shape mismatch: expected ({T_SIM+1},{G}), '
                f'got {source_traj.shape} (wid={wid}, rep={rep_class})')
        diff = np.max(np.abs(base_traj - source_traj))
        sanity_ok = np.allclose(base_traj, source_traj, atol=1e-12)
        if not sanity_ok:
            sanity_fails += 1
            print(f'  SANITY FAIL: {wid} rep{rep_class} max|diff|={diff:.2e}')

        # ── Perturbed: two-segment stitched ──
        pre_traj, pre_mask = _run_branch(W, b_vec, X0, T_FORK)       # t=0..500  (501 pts)
        post_traj, post_mask = _run_branch(W_post, b_vec, pre_traj[-1], T_POST)  # t=500..1000
        pert_traj, pert_mask = _stitch_perturbed_traj(
            pre_traj, pre_mask, post_traj, post_mask)

        pid = plan['perturbation_id']
        fork_results[pid] = {
            'perturbation_id':         pid,
            'source_world_id':         wid,
            'X0':                      X0,
            'X_fork':                  X_fork_val,
            'baseline_traj':           base_traj,
            'perturbed_traj':          pert_traj,
            'baseline_clip_mask':      base_mask,
            'perturbed_clip_mask':     pert_mask,
            'late_window_start':       T_SIM - LATE_WINDOW,
            'late_window_end':         T_SIM,
            'sanity_passed':           sanity_ok,
        }

        # Progress indicator
        if (idx + 1) % 50 == 0 or idx + 1 == total:
            print(f'  [{idx+1}/{total}] {pid}')

    # Summary
    print(f'\n  Sanity check: {sanity_fails}/{total} baseline trajectories differ from run_005')
    if sanity_fails == 0:
        print(f'  All baseline re-simulations match run_005 source trajectories')

    # Save each fork result as .pt file
    print(f'  Saving {len(fork_results)} fork trajectories to {TRAJ_DIR}/...')
    for pid, result in fork_results.items():
        out_path = os.path.join(TRAJ_DIR, f'{pid}_fork.pt')
        save_dict = {k: (torch.from_numpy(v) if isinstance(v, np.ndarray) else v)
                     for k, v in result.items()}
        torch.save(save_dict, out_path)

    print(f'  DONE: {len(fork_results)} fork pairs simulated '
          f'(×2 branches = {len(fork_results)*2} trajectories)')
    return fork_results


# ═══════════════════════════════════════════════════════════════════
# §5  Trajectory-level regime classification
# ═══════════════════════════════════════════════════════════════════
#
# Delegates to util.py primitives (shared with run_grn_perturbation):
#   detect_equilibrium, analyze_stability, analyze_oscillation, classify_attractor
#
# clip_count / total_clips are derived from clip_mask stored in §4.

# ── Type A-G → document regime label mapping ──
_TYPE_TO_PRIMARY: Dict[str, str] = {
    'Type A': 'Convergent',
    'Type B': 'Convergent',
    'Type C': 'Convergent',
    'Type D': 'Sustained oscillatory',
    'Type E': 'Divergent',
    'Type F': 'Collapse',
    'Type G': 'Ambiguous',
}

_TYPE_TO_SECONDARY: Dict[str, str] = {
    'Type A': 'fast convergence',
    'Type B': 'slow convergence',
    'Type C': 'damped oscillatory transient',
    'Type D': '',
    'Type E': 'runaway divergence',
    'Type F': 'numerical collapse',
    'Type G': '',
}


def _compute_equilibrium(X_traj: np.ndarray, converged: bool,
                         conv_time: int) -> np.ndarray:
    """Equilibrium vector: mean over convergence window or late window.

    - Converged: mean of X_traj[conv_time : conv_time + CONVERGENCE_WINDOW]
    - Not converged: mean of last CONVERGENCE_WINDOW steps
    """
    if converged and conv_time >= 0:
        end = min(conv_time + util.CONVERGENCE_WINDOW, X_traj.shape[0])
        return X_traj[conv_time:end].mean(axis=0)
    return X_traj[-util.CONVERGENCE_WINDOW:].mean(axis=0)


def _classify_one_branch(X_traj: np.ndarray, clip_mask: np.ndarray,
                         t_fork: int = 0
                         ) -> Dict[str, Any]:
    """Run full detection + classification chain for a single trajectory branch.

    Parameters
    ----------
    X_traj   : (T_total+1, G) full trajectory
    clip_mask: (T_total+1, G) clip mask
    t_fork   : intervention time. If > 0, slice to [t_fork:] before analysis
               (so post-perturb equilibrium is computed from post-fork segment only).
    """
    if t_fork > 0:
        X_traj    = X_traj[t_fork:]
        clip_mask = clip_mask[t_fork:]
    # conv_time from detect_equilibrium is 0-indexed within X_traj (sliced if t_fork>0).

    X_t = torch.from_numpy(X_traj.astype(np.float64))
    clip_count = torch.from_numpy(clip_mask.sum(axis=1).astype(np.int64))
    total_clips = int(clip_mask.sum())

    eq = util.detect_equilibrium(X_t)
    st = util.analyze_stability(X_t, clip_count, total_clips, G=G)
    osc = util.analyze_oscillation(X_t, eq['converged'], eq['convergence_time'])
    atype = util.classify_attractor(eq, st, osc)

    # conv_time comes back indexed to X_t = post-fork segment if sliced.
    eq_vec = _compute_equilibrium(X_traj, eq['converged'], eq['convergence_time'])
    # ↑ still operates on the (possibly sliced) X_traj.

    return {
        'attractor_type':        atype,
        'primary_regime':        _TYPE_TO_PRIMARY[atype],
        'secondary_label':       _TYPE_TO_SECONDARY[atype],
        'converged':             eq['converged'],
        'convergence_time':      eq['convergence_time'],
        'equilibrium':           eq_vec,
        'equilibrium_magnitude': eq['equilibrium_magnitude'],
        'equilibrium_sparsity':  eq['equilibrium_sparsity'],
        'divergence_existence':  st['divergence_existence'],
        'divergence_time':       st.get('divergence_time'),
        'max_expression':        st['max_expression'],
        'numerical_collapse':    st['numerical_collapse'],
        'clipping_dominated':    st['clipping_dominated'],
        'total_clips':           st['total_clips'],
        'oscillation_type':      osc['oscillation_type'],
    }


def s5_classify_regimes(fork_results):
    """Classify primary regime for each baseline and perturbed branch.

    Returns
    -------
    regime_results : dict
        perturbation_id → {
            perturbation_id, source_world_id,
            base_primary_regime,  base_secondary_label,  base_converged, ...
            pert_primary_regime,  pert_secondary_label,  pert_converged, ...
        }
    """
    print(f'\n{"="*70}')
    print(f'§5  Trajectory-level Regime Classification')
    print(f'{"="*70}')
    print(f'  Thresholds: EPSILON={util.EPSILON}, CONVERGENCE_WINDOW={util.CONVERGENCE_WINDOW}, '
          f'DIVERGENCE={util.DIVERGENCE_THRESHOLD}, COLLAPSE={util.COLLAPSE_THRESHOLD}')
    print(f'  Delegating to util.py: detect_equilibrium, analyze_stability, '
          f'analyze_oscillation, classify_attractor')

    regime_results: Dict[str, Dict] = {}
    total = len(fork_results)
    counts: Dict[str, int] = {}

    for idx, (pid, fr) in enumerate(fork_results.items()):
        wid = fr['source_world_id']

        # baseline: full trajectory
        base = _classify_one_branch(fr['baseline_traj'], fr['baseline_clip_mask'])
        # perturbed: restrict to POST-FORK segment so equilibrium reflects post-perturb dynamics
        pert = _classify_one_branch(fr['perturbed_traj'], fr['perturbed_clip_mask'],
                                    t_fork=T_FORK)

        regime_results[pid] = {
            'perturbation_id':         pid,
            'source_world_id':         wid,
            # baseline
            'base_primary_regime':     base['primary_regime'],
            'base_secondary_label':    base['secondary_label'],
            'base_attractor_type':     base['attractor_type'],
            'base_converged':          base['converged'],
            'base_convergence_time':   base['convergence_time'],
            'base_equilibrium':        base['equilibrium'],
            'base_eq_magnitude':       base['equilibrium_magnitude'],
            'base_eq_sparsity':        base['equilibrium_sparsity'],
            'base_divergence':         base['divergence_existence'],
            'base_divergence_time':    base['divergence_time'],
            'base_collapse':           base['numerical_collapse'],
            'base_clipping_dominated': base['clipping_dominated'],
            'base_total_clips':        base['total_clips'],
            'base_oscillation_type':   base['oscillation_type'],
            # perturbed
            'pert_primary_regime':     pert['primary_regime'],
            'pert_secondary_label':    pert['secondary_label'],
            'pert_attractor_type':     pert['attractor_type'],
            'pert_converged':          pert['converged'],
            'pert_convergence_time':   pert['convergence_time'],
            'pert_equilibrium':        pert['equilibrium'],
            'pert_eq_magnitude':       pert['equilibrium_magnitude'],
            'pert_eq_sparsity':        pert['equilibrium_sparsity'],
            'pert_divergence':         pert['divergence_existence'],
            'pert_divergence_time':    pert['divergence_time'],
            'pert_collapse':           pert['numerical_collapse'],
            'pert_clipping_dominated': pert['clipping_dominated'],
            'pert_total_clips':        pert['total_clips'],
            'pert_oscillation_type':   pert['oscillation_type'],
        }

        # Tally baseline attractor types
        at = base['attractor_type']
        counts[at] = counts.get(at, 0) + 1

        if (idx + 1) % 100 == 0 or idx + 1 == total:
            print(f'  [{idx+1}/{total}] classified')

    print(f'\n  Baseline attractor type distribution:')
    for t in ['Type A', 'Type B', 'Type C', 'Type D', 'Type E', 'Type F', 'Type G']:
        if t in counts:
            print(f'    {t}: {counts[t]}')

    print(f'  DONE: {len(regime_results)} trajectories classified')
    return regime_results


# ═══════════════════════════════════════════════════════════════════
# §6  Perturbation outcome assignment
# ═══════════════════════════════════════════════════════════════════

LANDING_SAME_THRESHOLD = 0.05   # normalized L2: robust vs landing_shift
EPSILON_NORM            = 1e-8  # numerical floor for norm denominators


def _late_window_state(X_traj: np.ndarray) -> np.ndarray:
    """Mean state vector over the late-window (last LATE_WINDOW steps)."""
    return X_traj[-LATE_WINDOW:].mean(axis=0)


def _compute_state_distances(v_pre: np.ndarray, v_post: np.ndarray
                             ) -> Dict[str, float]:
    """Normalized L2, cosine, and correlation distances between two vectors."""
    norm_pre = float(np.linalg.norm(v_pre))
    norm_post = float(np.linalg.norm(v_post))
    l2_raw = float(np.linalg.norm(v_post - v_pre))
    l2_normed = l2_raw / max(norm_pre, EPSILON_NORM)

    dot = float(np.dot(v_pre, v_post))
    mag = norm_pre * max(norm_post, EPSILON_NORM)
    # cosine distance: 0 when parallel, 2 when opposite
    cosine_dist = 1.0 - dot / mag

    # Pearson correlation distance: ∈ [0,1], 0 when perfectly correlated
    corr_dist = 0.5 * (1.0 - float(np.corrcoef(v_pre, v_post)[0, 1]))

    return {
        'pre_reference_state_norm':     norm_pre,
        'post_final_state_norm':        norm_post,
        'pre_post_normalized_L2':       l2_normed,
        'pre_post_cosine_distance':     cosine_dist,
        'pre_post_correlation_distance': corr_dist,
    }


def _decide_outcome(base_prim: str, pert_prim: str,
                    dists: Dict[str, float]) -> str:
    """Map (baseline regime, perturbed regime, distances) → outcome label."""
    if base_prim == 'Ambiguous' or pert_prim == 'Ambiguous':
        return 'ambiguous'
    if pert_prim == 'Collapse':
        return 'collapse'

    if base_prim == pert_prim:
        # same regime — distinguish robust vs landing_shift by distance
        if dists['pre_post_normalized_L2'] < LANDING_SAME_THRESHOLD:
            return 'robust'
        return 'landing_shift'
    else:
        return 'regime_transition'


def s6_assign_outcomes(traj_results, fork_results):
    """Compare baseline vs perturbed branch outcomes, assign perturbation_outcome.

    Parameters
    ----------
    traj_results : dict  (from s5_classify_regimes)
        perturbation_id → {base_*, pert_*, ...}
    fork_results : dict  (from s4_run_fork_simulations)
        perturbation_id → {perturbed_clip_mask, ...}

    Returns
    -------
    outcomes : dict
        perturbation_id → { outcome fields }
    """
    print(f'\n{"="*70}')
    print(f'§6  Perturbation Outcome Assignment')
    print(f'{"="*70}')
    print(f'  landing_same_threshold = {LANDING_SAME_THRESHOLD}')
    print(f'  epsilon_norm           = {EPSILON_NORM}')
    print(f'  late_window            = last {LATE_WINDOW} steps')

    outcomes: Dict[str, Dict] = {}
    tally: Dict[str, int] = {}

    for idx, (pid, rr) in enumerate(traj_results.items()):
        base_conv   = rr['base_converged']
        pert_conv   = rr['pert_converged']
        base_conv_t = rr['base_convergence_time']
        pert_conv_t = rr['pert_convergence_time']

        dists = _compute_state_distances(
            rr['base_equilibrium'], rr['pert_equilibrium'])

        result = _decide_outcome(
            rr['base_primary_regime'], rr['pert_primary_regime'], dists)

        regime_transition = (rr['base_primary_regime'] != rr['pert_primary_regime'])

        # clipping frequencies: baseline & perturbed, full post-fork period (t=501..1000)
        base_mask = fork_results[pid]['baseline_clip_mask']
        pert_mask = fork_results[pid]['perturbed_clip_mask']
        base_clip_freq = base_mask[T_FORK + 1:].sum() / (T_POST * G)
        post_clip_freq = pert_mask[T_FORK + 1:].sum() / (T_POST * G)

        outcomes[pid] = {
            'perturbation_id':             pid,
            'source_world_id':             rr['source_world_id'],
            'perturbation_outcome':        result,
            'regime_transition':           regime_transition,
            'base_primary_regime':         rr['base_primary_regime'],
            'base_secondary_label':        rr['base_secondary_label'],
            'base_attractor_type':         rr['base_attractor_type'],
            'pert_primary_regime':         rr['pert_primary_regime'],
            'pert_secondary_label':        rr['pert_secondary_label'],
            'pert_attractor_type':         rr['pert_attractor_type'],
            # convergence
            'converged_pre':               base_conv,
            'converged_post':              pert_conv,
            'convergence_time_pre':        base_conv_t,
            'convergence_time_post':       pert_conv_t,
            'divergence_time_post':        rr['pert_divergence_time'],
            # distances
            **dists,
            # clipping
            'base_clipping_frequency':     base_clip_freq,
            'post_clipping_frequency':     post_clip_freq,
            # oscillation
            'base_oscillation_type':       rr['base_oscillation_type'],
            'pert_oscillation_type':       rr['pert_oscillation_type'],
            # flags
            'base_divergence':             rr['base_divergence'],
            'pert_divergence':             rr['pert_divergence'],
            'base_collapse':               rr['base_collapse'],
            'pert_collapse':               rr['pert_collapse'],
            # equilibrium vectors (for §7 clustering)
            'base_equilibrium':            rr['base_equilibrium'],
            'pert_equilibrium':            rr['pert_equilibrium'],
        }

        tally[result] = tally.get(result, 0) + 1

        if (idx + 1) % 100 == 0 or idx + 1 == len(traj_results):
            print(f'  [{idx+1}/{len(traj_results)}] assigned')

    print(f'\n  Outcome distribution:')
    for label in ['robust', 'landing_shift', 'regime_transition', 'collapse', 'ambiguous']:
        if label in tally:
            print(f'    {label}: {tally[label]}')
    print(f'  DONE: {len(outcomes)} perturbation outcomes assigned')
    return outcomes


# ═══════════════════════════════════════════════════════════════════
# §7  State landing analysis (within-world clustering)
# ═══════════════════════════════════════════════════════════════════

# ── Active-set thresholds (run_006 Analysis Plan §13.2) ──
ZERO_THRESHOLD                   = 1e-8   # abundance below this counts as "zero"
ACTIVE_ZERO_FREQ_THRESHOLD       = 0.5    # gene is zero in >50% of late-window steps
ACTIVE_CLIPPING_FREQ_THRESHOLD   = 0.1    # gene clips in >10% of late-window steps


def _encode_gene_set(gene_indices: np.ndarray) -> str:
    """Encode a set of gene indices as a compact string for CSV storage."""
    if len(gene_indices) == 0:
        return ''
    return ','.join(str(int(g)) for g in sorted(gene_indices))


def _compute_active_sets(pert_traj: np.ndarray, clip_mask: np.ndarray
                         ) -> Tuple[str, str]:
    """Compute active_zero_set_id and active_clipping_set_id from perturbed branch.

    Per run_006 §13.2: late-window (~50-100 steps) frequency-based.
      - active_zero_set: abundance < 1e-8 in >50% of late-window steps
      - active_clipping_set: clipped in >10% of late-window steps
    """
    late_traj = pert_traj[-LATE_WINDOW:]    # (LATE_WINDOW, G)
    late_mask = clip_mask[-LATE_WINDOW:]    # (LATE_WINDOW, G)

    # genes that are near-zero in >= 50% of late-window steps
    zero_freq = (late_traj < ZERO_THRESHOLD).mean(axis=0)   # (G,)
    zero_genes = np.where(zero_freq >= ACTIVE_ZERO_FREQ_THRESHOLD)[0]

    # genes that trigger clipping in >= 10% of late-window steps
    clip_freq = late_mask.mean(axis=0)                      # (G,)
    clip_genes = np.where(clip_freq >= ACTIVE_CLIPPING_FREQ_THRESHOLD)[0]

    return _encode_gene_set(zero_genes), _encode_gene_set(clip_genes)


# ── Supplementary metrics for non-convergent regimes (§12.2-12.4) ──

def _compute_oscillatory_metrics(pert_traj: np.ndarray) -> Dict[str, object]:
    """Compute oscillatory summary metrics from perturbed trajectory late-window.

    Returns
    -------
    dict with keys:
      - oscillatory_amplitude_per_gene : (G,)  peak-to-peak amplitude per gene
      - oscillatory_mean_amplitude     : float  mean amplitude across genes
      - oscillatory_period             : int    dominant period in steps, -1 if undetectable
    """
    late = pert_traj[-LATE_WINDOW:]                          # (LATE_WINDOW, G)

    # Peak-to-peak amplitude per gene
    amp_per_gene = late.max(axis=0) - late.min(axis=0)       # (G,)

    # Mean state for autocorrelation
    mean_curve = late.mean(axis=1)                           # (LATE_WINDOW,)
    mean_curve -= mean_curve.mean()
    if np.std(mean_curve) < 1e-12:
        period = -1
    else:
        # Autocorrelation
        ac = np.correlate(mean_curve, mean_curve, mode='full')
        ac = ac[ac.size // 2:]                               # non-negative lags only
        ac /= ac[0]                                          # normalize

        # Find first peak after lag 0 (skip immediate neighbors)
        for lag in range(2, LATE_WINDOW - 1):
            if ac[lag] > ac[lag - 1] and ac[lag] > ac[lag + 1] and ac[lag] > 0.2:
                period = lag
                break
        else:
            period = -1

    return {
        'oscillatory_amplitude_per_gene': amp_per_gene,
        'oscillatory_mean_amplitude':     float(np.mean(amp_per_gene)),
        'oscillatory_period':             period,
    }


def _compute_divergent_metrics(pert_traj: np.ndarray) -> Dict[str, float]:
    """Compute divergent trajectory magnitude growth.

    Returns
    -------
    dict with keys:
      - divergent_magnitude_growth : float  log2(norm_late / norm_early)
    """
    # Early window: first LATE_WINDOW steps after fork
    early = pert_traj[T_FORK + 1: T_FORK + 1 + LATE_WINDOW]   # (LATE_WINDOW, G)
    late  = pert_traj[-LATE_WINDOW:]                           # (LATE_WINDOW, G)

    norm_early = float(np.linalg.norm(early.mean(axis=0)))
    norm_late  = float(np.linalg.norm(late.mean(axis=0)))

    if norm_early < 1e-8:
        growth = float('inf') if norm_late > 1e-8 else 0.0
    else:
        growth = float(np.log2(max(norm_late, 1e-8) / norm_early))

    return {'divergent_magnitude_growth': growth}


def _compute_collapse_metrics(pert_equilibrium: np.ndarray) -> Dict[str, object]:
    """Compute collapse trajectory metrics.

    Returns
    -------
    dict with keys:
      - collapse_zero_fraction : float   fraction of genes at zero
      - collapse_type          : str     'global' (>=50%) or 'partial'
    """
    zero_frac = float((pert_equilibrium < ZERO_THRESHOLD).mean())
    ctype = 'global' if zero_frac >= 0.5 else 'partial'
    return {
        'collapse_zero_fraction': zero_frac,
        'collapse_type':          ctype,
    }


def _build_landing_results(
    outcomes: Dict[str, dict],
) -> Dict[str, list]:
    """Build cluster-level summary and per-gene statistics from outcomes.

    Returns dict with keys 'summary_rows' and 'vector_rows'.
    """
    # Group convergent trajectories by (world, mode, repeat, state_landing_id)
    clusters: Dict[Tuple[str, str, str, int], List[Tuple[str, dict]]] = {}
    for pid, oc in outcomes.items():
        lid = oc.get('state_landing_id_post')
        if isinstance(lid, int) and lid >= 0:
            ctx = oc.get('perturbation_context', ('', ''))
            mode = ctx[0] if isinstance(ctx, (list, tuple)) else ''
            rpt  = ctx[1] if isinstance(ctx, (list, tuple)) and len(ctx) > 1 else ''
            clusters.setdefault((oc['source_world_id'], mode, rpt, lid), []).append((pid, oc))

    # Count total convergent trajectories per perturbation context (for basin_fraction)
    ctx_conv_total: Dict[Tuple[str, str, str], int] = {}
    for oc in outcomes.values():
        if oc['pert_primary_regime'] == 'Convergent':
            wid = oc['source_world_id']
            ctx = oc.get('perturbation_context', ('', ''))
            mode = ctx[0] if isinstance(ctx, (list, tuple)) else ''
            rpt  = ctx[1] if isinstance(ctx, (list, tuple)) and len(ctx) > 1 else ''
            key = (wid, mode, rpt)
            ctx_conv_total[key] = ctx_conv_total.get(key, 0) + 1

    summary_rows = []
    vector_rows = []

    for (wid, mode, rpt, lid), members in sorted(clusters.items()):
        n = len(members)
        pids = [pid for pid, _ in members]
        eq_vectors = np.stack([oc['pert_equilibrium'] for _, oc in members])      # (n, G)
        base_vectors = np.stack([oc['base_equilibrium'] for _, oc in members])    # (n, G)

        # Cluster-level per-gene statistics
        centroid  = eq_vectors.mean(axis=0)                                       # (G,)
        median_v  = np.median(eq_vectors, axis=0)                                 # (G,)
        sd_v      = eq_vectors.std(axis=0, ddof=1) if n > 1 else np.zeros(G)     # (G,)

        # Medoid: member closest to the centroid
        dists_to_centroid = np.linalg.norm(eq_vectors - centroid, axis=1)
        medoid_idx = int(np.argmin(dists_to_centroid))
        medoid_pid = pids[medoid_idx]
        medoid_vec = eq_vectors[medoid_idx]                                      # (G,)

        # Mean pairwise L2 distance within the cluster
        if n > 1:
            pw = []
            for i in range(n):
                for j in range(i + 1, n):
                    pw.append(float(np.linalg.norm(eq_vectors[i] - eq_vectors[j])))
            mean_pw_dist = np.mean(pw)
        else:
            mean_pw_dist = 0.0

        # Nearest pre-state distance: min normalized-L2 from centroid to any member's base_equilibrium
        nearest_pre = float(np.min([
            np.linalg.norm(centroid - bv) / max(float(np.linalg.norm(bv)), EPSILON_NORM)
            for bv in base_vectors
        ]))

        basin_frac = n / ctx_conv_total.get((wid, mode, rpt), n)

        # Save representative vector as .pt file
        ctx_slug = f'{mode}' if not rpt else f'{mode}_r{rpt}'
        rep_path = os.path.join(TRAJ_DIR, f'{wid}_{ctx_slug}_landing{lid}_rep.pt')
        torch.save({
            'source_world_id':           wid,
            'state_landing_id':          lid,
            'medoid_perturbation_id':    medoid_pid,
            'representative_vector':     torch.from_numpy(medoid_vec),
            'cluster_mean':              torch.from_numpy(centroid),
            'cluster_median':            torch.from_numpy(median_v),
            'cluster_sd':                torch.from_numpy(sd_v),
            'n_members':                 n,
        }, rep_path)

        # Active sets from the medoid trajectory
        medoid_oc: dict = next((oc for pid, oc in members if pid == medoid_pid), {})

        summary_rows.append({
            'source_world_id':                    wid,
            'perturbation_id':                    medoid_pid,
            'perturbation_mode':                  mode,
            'repeat_index':                       rpt,
            'state_landing_id':                   lid,
            'post_primary_regime':                'Convergent',
            'n_initial_conditions':               n,
            'basin_fraction':                     round(basin_frac, 6),
            'representative_state_vector_path':   rep_path,
            'representative_state_norm':          round(float(np.linalg.norm(medoid_vec)), 6),
            'mean_pairwise_distance_within_cluster': round(mean_pw_dist, 6),
            'nearest_pre_state_distance':         round(nearest_pre, 6),
            'active_zero_set_id':                 medoid_oc.get('post_active_zero_set_id', ''),
            'active_clipping_set_id':             medoid_oc.get('post_active_clipping_set_id', ''),
            'interpretation_note':                '',
        })

        # Per-gene rows for representative_state_vectors.csv
        for g in range(G):
            vector_rows.append({
                'source_world_id':         wid,
                'perturbation_id':         medoid_pid,
                'perturbation_mode':       mode,
                'repeat_index':            rpt,
                'state_landing_id':        lid,
                'gene_id':                 g,
                'representative_abundance': round(float(medoid_vec[g]), 8),
                'cluster_mean_abundance':   round(float(centroid[g]), 8),
                'cluster_median_abundance': round(float(median_v[g]), 8),
                'cluster_sd_abundance':     round(float(sd_v[g]), 8),
            })

    return {'summary_rows': summary_rows, 'vector_rows': vector_rows}


def s7_state_landing_clustering(outcomes, fork_results, pert_plan):
    """Cluster post-perturb convergent states within each perturbation context.

    Clustering is done independently within each (world, perturbation_mode, repeat_index)
    group, so state_landing_id is only meaningful within that context (per zyx doc).

    Modifies outcomes in-place, adding:
      - perturbation_context         (mode, repeat) tuple for downstream use
      - state_landing_id_post        (int for convergent, str for non-convergent)
      - base_active_zero_set_id      (str, baseline late-window)
      - base_active_clipping_set_id  (str, baseline late-window)
      - post_active_zero_set_id      (str, perturbed post-fork late-window)
      - post_active_clipping_set_id  (str, perturbed post-fork late-window)

    Returns
    -------
    state_landing_results : dict
        'summary_rows' → list of dicts for state_landing_summary.csv
        'vector_rows'  → list of dicts for representative_state_vectors.csv
    """
    print(f'\n{"="*70}')
    print(f'§7  State Landing Analysis (within-perturbation-context clustering)')
    print(f'{"="*70}')
    print(f'  CLUSTER_THRESHOLD = {CLUSTER_THRESHOLD}')
    print(f'  Active set window = {LATE_WINDOW} steps')
    print(f'    zero_threshold={ZERO_THRESHOLD}, zero_freq={ACTIVE_ZERO_FREQ_THRESHOLD}, '
          f'clip_freq={ACTIVE_CLIPPING_FREQ_THRESHOLD}')
    print(f'  Computing active sets for BOTH baseline and perturbed branches')

    # Build pid → (mode, repeat) mapping from pert_plan; store in outcomes
    pid_to_ctx: Dict[str, Tuple[str, str]] = {}
    for plan in pert_plan:
        pid = plan['perturbation_id']
        pid_to_ctx[pid] = (plan.get('perturbation_mode', ''), plan.get('repeat_index', ''))
    for pid, oc in outcomes.items():
        oc['perturbation_context'] = pid_to_ctx.get(pid, ('', ''))

    # ── Group outcomes by (world, perturbation_mode, repeat_index) ──
    by_ctx: Dict[Tuple[str, str, str], List[Tuple[str, dict]]] = {}
    for pid, oc in outcomes.items():
        mode, rpt = pid_to_ctx.get(pid, ('', ''))
        key = (oc['source_world_id'], mode, rpt)
        by_ctx.setdefault(key, []).append((pid, oc))

    n_clustered = 0

    for (wid, mode, rpt) in sorted(by_ctx):
        entries = by_ctx[(wid, mode, rpt)]
        conv = [(pid, oc) for pid, oc in entries
                if oc['pert_primary_regime'] == 'Convergent']
        nonconv = [(pid, oc) for pid, oc in entries
                   if oc['pert_primary_regime'] != 'Convergent']

        # ── Convergent: cluster equilibrium vectors ──
        if conv:
            eq_vectors = [oc['pert_equilibrium'] for _, oc in conv]
            labels, _, _ = util.hierarchical_cluster_relative_l2(
                eq_vectors, threshold=CLUSTER_THRESHOLD)
            for idx, (pid, oc) in enumerate(conv):
                outcomes[pid]['state_landing_id_post'] = int(labels[idx])
                n_clustered += 1
        # (non-convergent trajectories handled below with string labels)

        # ── Non-convergent: string labels ──
        _REGIME_LABEL: Dict[str, str] = {
            'Convergent':             'convergent',
            'Sustained oscillatory':  'oscillatory',
            'Divergent':              'divergent',
            'Collapse':               'collapse',
            'Ambiguous':              'ambiguous',
        }
        for pid, oc in nonconv:
            label = _REGIME_LABEL.get(oc['pert_primary_regime'],
                                       oc['pert_primary_regime'].lower().replace(' ', '_'))
            outcomes[pid]['state_landing_id_post'] = label

        # ── Per-trajectory active sets + supplementary metrics (all regimes) ──
        for pid, oc in entries:
            fr = fork_results.get(pid)
            if fr is not None:
                # Baseline branch active sets
                bzid, bcid = _compute_active_sets(
                    fr['baseline_traj'], fr['baseline_clip_mask'])
                outcomes[pid]['base_active_zero_set_id'] = bzid
                outcomes[pid]['base_active_clipping_set_id'] = bcid

                # Perturbed branch active sets
                zid, cid = _compute_active_sets(
                    fr['perturbed_traj'], fr['perturbed_clip_mask'])
                outcomes[pid]['post_active_zero_set_id'] = zid
                outcomes[pid]['post_active_clipping_set_id'] = cid

                # ── Supplementary metrics for non-convergent regimes (§12.2-12.4) ──
                regime = oc['pert_primary_regime']
                if regime == 'Sustained oscillatory':
                    osc = _compute_oscillatory_metrics(fr['perturbed_traj'])
                    outcomes[pid].update(osc)
                elif regime == 'Divergent':
                    div = _compute_divergent_metrics(fr['perturbed_traj'])
                    outcomes[pid].update(div)
                elif regime == 'Collapse':
                    col = _compute_collapse_metrics(oc['pert_equilibrium'])
                    outcomes[pid].update(col)
            else:
                outcomes[pid]['base_active_zero_set_id'] = ''
                outcomes[pid]['base_active_clipping_set_id'] = ''
                outcomes[pid]['post_active_zero_set_id'] = ''
                outcomes[pid]['post_active_clipping_set_id'] = ''

        # Per-context summary
        n_conv = len(conv)
        n_labels = len(set(int(labels[i]) for i in range(len(conv)))) if conv else 0
        ctx_label = f'{mode}' if not rpt else f'{mode}_r{rpt}'
        print(f'  [{wid} / {ctx_label}]  {len(entries)} trajectories: '
              f'{n_conv} convergent → {n_labels} state landings')

    # ── Build cluster-level statistics ──
    landing = _build_landing_results(outcomes)

    print(f'\n  DONE: {n_clustered} convergent trajectories clustered')
    print(f'        {len(landing["summary_rows"])} state landing clusters '
          f'({len(landing["vector_rows"])} per-gene rows)')
    return landing


# ═══════════════════════════════════════════════════════════════════
# §8  Table generation
# ═══════════════════════════════════════════════════════════════════

def _write_csv(path: str, columns: list, rows: list):
    """Write a list of dicts to CSV with specified column order."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = _csv.DictWriter(f, fieldnames=columns, extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)


def s8_build_tables(screening, pert_plan, outcomes, state_landing_results,
                    selected_worlds):
    """Generate all output CSVs: trajectory_outcomes, perturbation_summary,
    state_landing_summary, representative_state_vectors, edge_property_outcome.

    Parameters
    ----------
    screening : dict           world_id → screening info (init_files, cluster_info, ...)
    pert_plan : list[dict]     perturbation plan list from §3
    outcomes  : dict           perturbation_id → outcome dict from §6+§7
    state_landing_results : dict  'summary_rows' and 'vector_rows' from §7
    selected_worlds : dict     world_id → {W, meta, rho_W, scc, ...} from §1

    Returns
    -------
    None (writes files to data/)
    """
    print(f'\n{"="*70}')
    print(f'§8  Table Generation')
    print(f'{"="*70}')

    DATA = os.path.join(PROJ_ROOT, 'data')

    # ── Precompute pert_plan lookup (pid → plan) ──
    plan_by_pid: Dict[str, dict] = {p['perturbation_id']: p for p in pert_plan}
    # ── Build screening lookup: (world_id, init_id) → init_type ──
    init_type_lookup: Dict[Tuple[str, str], str] = {}
    for wid, scr in screening.items():
        for fdict in scr.get('init_files', []):
            init_type_lookup[(wid, fdict['init_id'])] = fdict['init_type']

    # ── world-level SCC & out/in degree caches ──
    world_scc: Dict[str, set] = {}
    world_outdeg: Dict[str, dict] = {}
    world_indeg: Dict[str, dict] = {}
    for wid, sw in selected_worlds.items():
        world_scc[wid] = set(sw.get('largest_scc', []))
        world_outdeg[wid] = {g: int((np.abs(sw['W'][g, :]) > 1e-12).sum())
                             for g in range(G)}
        world_indeg[wid] = {g: int((np.abs(sw['W'][:, g]) > 1e-12).sum())
                            for g in range(G)}

    # ═══════════════════════════════════════════════════════════════
    # 4.4  trajectory_outcomes.csv
    # ═══════════════════════════════════════════════════════════════
    print('  Building trajectory_outcomes.csv ...')
    traj_cols = [
        'trajectory_id',
        'source_world_id',
        'perturbation_id',
        'initial_condition_id',
        'init_type',
        'pre_perturb_state_class_id',
        'is_representative_trajectory',
        'representative_initial_condition_id',
        't_perturb',
        'T_post',
        'T_total',
        'baseline_trajectory_id',
        'perturbed_trajectory_id',
        'baseline_primary_regime',
        'baseline_secondary_label',
        'baseline_late_window_start',
        'baseline_late_window_end',
        'pre_reference_state_source',
        'pre_primary_regime',
        'pre_secondary_label',
        'post_primary_regime',
        'post_secondary_label',
        'regime_transition',
        'perturbation_outcome',
        'converged_pre',
        'converged_post',
        'convergence_time_pre',
        'convergence_time_post',
        'divergence_time_post',
        'pre_reference_state_norm',
        'post_final_state_norm',
        'pre_post_normalized_L2',
        'pre_post_cosine_distance',
        'pre_post_correlation_distance',
        'post_clipping_frequency',
        'post_active_zero_set_id',
        'post_active_clipping_set_id',
        'state_landing_id_post',
        'notes',
    ]

    traj_rows = []
    for pid, oc in sorted(outcomes.items()):
        plan = plan_by_pid.get(pid, {})
        wid = oc['source_world_id']
        init_id = plan.get('representative_init_id', '')
        itype = init_type_lookup.get((wid, init_id), '')

        traj_rows.append({
            'trajectory_id':                      pid,
            'source_world_id':                    wid,
            'perturbation_id':                    pid,
            'initial_condition_id':              init_id,
            'init_type':                          itype,
            'pre_perturb_state_class_id':        plan.get('representative_class_id', ''),
            'is_representative_trajectory':       True,
            'representative_initial_condition_id': init_id,
            't_perturb':                          T_FORK,
            'T_post':                             T_POST,
            'T_total':                            T_TOTAL,
            'baseline_trajectory_id':             pid,
            'perturbed_trajectory_id':            pid,
            'baseline_primary_regime':            oc.get('base_primary_regime', ''),
            'baseline_secondary_label':           oc.get('base_secondary_label', ''),
            'baseline_late_window_start':         T_TOTAL - LATE_WINDOW,
            'baseline_late_window_end':           T_TOTAL,
            'pre_reference_state_source':         'baseline_late_window',
            'pre_primary_regime':                 oc.get('base_primary_regime', ''),
            'pre_secondary_label':                oc.get('base_secondary_label', ''),
            'post_primary_regime':                oc.get('pert_primary_regime', ''),
            'post_secondary_label':               oc.get('pert_secondary_label', ''),
            'regime_transition':                  oc.get('regime_transition', ''),
            'perturbation_outcome':              oc.get('perturbation_outcome', ''),
            'converged_pre':                      oc.get('converged_pre', ''),
            'converged_post':                     oc.get('converged_post', ''),
            'convergence_time_pre':               oc.get('convergence_time_pre', ''),
            'convergence_time_post':              oc.get('convergence_time_post', ''),
            'divergence_time_post':               oc.get('divergence_time_post', ''),
            'pre_reference_state_norm':           oc.get('pre_reference_state_norm', ''),
            'post_final_state_norm':              oc.get('post_final_state_norm', ''),
            'pre_post_normalized_L2':            oc.get('pre_post_normalized_L2', ''),
            'pre_post_cosine_distance':          oc.get('pre_post_cosine_distance', ''),
            'pre_post_correlation_distance':     oc.get('pre_post_correlation_distance', ''),
            'post_clipping_frequency':            oc.get('post_clipping_frequency', ''),
            'post_active_zero_set_id':            oc.get('post_active_zero_set_id', ''),
            'post_active_clipping_set_id':        oc.get('post_active_clipping_set_id', ''),
            'state_landing_id_post':              oc.get('state_landing_id_post', ''),
            'notes':                              '',
        })

    _write_csv(os.path.join(DATA, 'trajectory_outcomes.csv'), traj_cols, traj_rows)
    print(f'    → {len(traj_rows)} rows written')

    # ═══════════════════════════════════════════════════════════════
    # 4.5  perturbation_summary.csv
    # ═══════════════════════════════════════════════════════════════
    print('  Building perturbation_summary.csv ...')
    pert_sum_cols = [
        'perturbation_id',
        'source_world_id',
        'perturbation_mode',
        'delete_k',
        'repeat_index',
        'N_init_screened',
        'N_pre_perturb_state_classes',
        'N_perturbed_representatives',
        'fraction_robust',
        'fraction_landing_shift',
        'fraction_regime_transition',
        'fraction_collapse',
        'fraction_ambiguous',
        'mean_state_landing_distance',
        'median_state_landing_distance',
        'delta_rho_W',
        'delta_mean_clipping_frequency',
        'mean_convergence_time_change',
        'dominant_outcome',
    ]

    # Aggregate by (source_world, perturbation_mode, repeat)
    groups = defaultdict(list)  # key=(wid, mode, repeat_idx) → list of dicts

    for plan in pert_plan:
        pid = plan['perturbation_id']
        wid = plan['source_world_id']
        mode = plan.get('perturbation_mode', '')
        rpt = plan.get('repeat_index', '')
        key = (wid, mode, rpt)
        groups[key].append({'pid': pid, 'plan': plan, 'oc': outcomes.get(pid, {})})

    pert_sum_rows = []
    for (wid, mode, rpt), members in sorted(groups.items()):
        n = len(members)
        ocs = [m['oc'] for m in members]

        # Outcome distribution
        outcomes_list = [oc.get('perturbation_outcome', '') for oc in ocs]
        def _frac(label):
            return round(sum(1 for o in outcomes_list if o == label) / n, 6) if n > 0 else 0.0

        # pre_post_normalized_L2 aggregation
        l2_vals = [oc['pre_post_normalized_L2'] for oc in ocs
                   if oc.get('pre_post_normalized_L2', '') != '']
        mean_l2 = round(float(np.mean(l2_vals)), 8) if l2_vals else ''
        median_l2 = round(float(np.median(l2_vals)), 8) if l2_vals else ''

        # delta clipping frequency aggregation
        delta_clips = []
        for oc in ocs:
            post_cf = oc.get('post_clipping_frequency')
            base_cf = oc.get('base_clipping_frequency')
            if post_cf is not None and base_cf is not None and post_cf != '' and base_cf != '':
                delta_clips.append(float(post_cf) - float(base_cf))
        mean_delta_clip = round(float(np.mean(delta_clips)), 8) if delta_clips else ''

        # convergence time change aggregation
        conv_deltas = []
        for oc in ocs:
            ct_pre = oc.get('convergence_time_pre')
            ct_post = oc.get('convergence_time_post')
            if ct_pre is not None and ct_post is not None:
                conv_deltas.append(ct_post - ct_pre)
        mean_conv_delta = round(float(np.mean(conv_deltas)), 4) if conv_deltas else ''

        # Dominant outcome (most frequent)
        outcome_counts = Counter(outcomes_list)
        dominant = outcome_counts.most_common(1)[0][0] if outcome_counts else ''

        # Representative perturbation_id: first member's pid
        rep_pid = members[0]['pid']

        pert_sum_rows.append({
            'perturbation_id':              rep_pid,
            'source_world_id':              wid,
            'perturbation_mode':            mode,
            'delete_k':                     members[0]['plan'].get('delete_k', K_DELETE),
            'repeat_index':                 rpt,
            'N_init_screened':              N_INIT_SCREEN,
            'N_pre_perturb_state_classes':  len(screening.get(wid, {}).get('cluster_info', {}).get('class_dict', {})),
            'N_perturbed_representatives':  n,
            'fraction_robust':              _frac('robust'),
            'fraction_landing_shift':       _frac('landing_shift'),
            'fraction_regime_transition':   _frac('regime_transition'),
            'fraction_collapse':            _frac('collapse'),
            'fraction_ambiguous':           _frac('ambiguous'),
            'mean_state_landing_distance':   mean_l2,
            'median_state_landing_distance': median_l2,
            'delta_rho_W':                   members[0]['plan'].get('delta_rho_W', ''),
            'delta_mean_clipping_frequency': mean_delta_clip,
            'mean_convergence_time_change':  mean_conv_delta,
            'dominant_outcome':              dominant,
        })

    _write_csv(os.path.join(DATA, 'perturbation_summary.csv'), pert_sum_cols, pert_sum_rows)
    print(f'    → {len(pert_sum_rows)} rows written')

    # ═══════════════════════════════════════════════════════════════
    # 4.6  state_landing_summary.csv
    # ═══════════════════════════════════════════════════════════════
    print('  Building state_landing_summary.csv ...')
    landing_cols = [
        'source_world_id',
        'perturbation_id',
        'perturbation_mode',
        'repeat_index',
        'state_landing_id',
        'post_primary_regime',
        'n_initial_conditions',
        'basin_fraction',
        'representative_state_vector_path',
        'representative_state_norm',
        'mean_pairwise_distance_within_cluster',
        'nearest_pre_state_distance',
        'active_zero_set_id',
        'active_clipping_set_id',
        'interpretation_note',
    ]
    summary_rows = state_landing_results.get('summary_rows', [])
    _write_csv(os.path.join(DATA, 'state_landing_summary.csv'), landing_cols, summary_rows)
    print(f'    → {len(summary_rows)} rows written')

    # ═══════════════════════════════════════════════════════════════
    # 4.7  representative_state_vectors.csv
    # ═══════════════════════════════════════════════════════════════
    print('  Building representative_state_vectors.csv ...')
    vector_cols = [
        'source_world_id',
        'perturbation_id',
        'perturbation_mode',
        'repeat_index',
        'state_landing_id',
        'gene_id',
        'representative_abundance',
        'cluster_mean_abundance',
        'cluster_median_abundance',
        'cluster_sd_abundance',
    ]
    vector_rows = state_landing_results.get('vector_rows', [])
    _write_csv(os.path.join(DATA, 'representative_state_vectors.csv'), vector_cols, vector_rows)
    print(f'    → {len(vector_rows)} rows written')

    # ═══════════════════════════════════════════════════════════════
    # 4.8  edge_property_outcome.csv
    # ═══════════════════════════════════════════════════════════════
    print('  Building edge_property_outcome.csv ...')
    edge_cols = [
        'source_world_id',
        'perturbation_id',
        'edge_id',
        'source_gene',
        'target_gene',
        'sign',
        'strength',
        'abs_strength',
        'source_outdegree',
        'target_indegree',
        'is_strong',
        'is_hub_outgoing',
        'is_in_largest_scc',
        'is_cycle_associated',
        'associated_fraction_robust',
        'associated_fraction_landing_shift',
        'associated_fraction_regime_transition',
        'associated_mean_state_distance',
    ]

    edge_rows = []
    for plan in pert_plan:
        pid = plan['perturbation_id']
        wid = plan['source_world_id']
        oc = outcomes.get(pid, {})

        src_str = plan.get('deleted_edge_sources', '')
        tgt_str = plan.get('deleted_edge_targets', '')
        sgn_str = plan.get('deleted_edge_signs', '')
        str_str = plan.get('deleted_edge_strengths', '')

        if not src_str or not tgt_str:
            continue

        sources = [int(x) for x in src_str.split(',')]
        targets = [int(x) for x in tgt_str.split(',')]
        signs   = sgn_str.split(',') if sgn_str else [''] * len(sources)
        strengths = [float(x) for x in str_str.split(',')] if str_str else [0.0] * len(sources)

        scc_set = world_scc.get(wid, set())
        outdeg = world_outdeg.get(wid, {})
        indeg  = world_indeg.get(wid, {})

        # Determine is_strong (top 5 strongest edges in this perturbation)
        abs_strs = [abs(s) for s in strengths]
        strong_threshold = sorted(abs_strs, reverse=True)[min(K_DELETE - 1, len(abs_strs) - 1)] if abs_strs else 0

        # Hub outgoing: edges from the highest-outdegree gene in this perturbation
        hub_src = max(set(sources), key=lambda g: outdeg.get(g, 0)) if sources else -1

        outcome = oc.get('perturbation_outcome', '')
        l2_dist = oc.get('pre_post_normalized_L2', '')

        for i in range(len(sources)):
            src = sources[i]
            tgt = targets[i]
            eid = f'{src}-{tgt}'

            edge_rows.append({
                'source_world_id':                   wid,
                'perturbation_id':                   pid,
                'edge_id':                           eid,
                'source_gene':                       src,
                'target_gene':                       tgt,
                'sign':                              signs[i],
                'strength':                          round(strengths[i], 8),
                'abs_strength':                      round(abs(strengths[i]), 8),
                'source_outdegree':                  outdeg.get(src, 0),
                'target_indegree':                   indeg.get(tgt, 0),
                'is_strong':                         abs(strengths[i]) >= strong_threshold,
                'is_hub_outgoing':                   src == hub_src,
                'is_in_largest_scc':                 src in scc_set and tgt in scc_set,
                'is_cycle_associated':               'NA',
                'associated_fraction_robust':         1.0 if outcome == 'robust' else 0.0,
                'associated_fraction_landing_shift':  1.0 if outcome == 'landing_shift' else 0.0,
                'associated_fraction_regime_transition': 1.0 if outcome == 'regime_transition' else 0.0,
                'associated_mean_state_distance':     l2_dist,
            })

    _write_csv(os.path.join(DATA, 'edge_property_outcome.csv'), edge_cols, edge_rows)
    print(f'    → {len(edge_rows)} rows written')

    print('  DONE: 5 CSV tables generated')


# ═══════════════════════════════════════════════════════════════════
# §9  Figure generation
# ═══════════════════════════════════════════════════════════════════

# ── Plotting constants ─────────────────────────────────────────────
REGIME_ORDER = ['Convergent', 'Sustained oscillatory', 'Divergent',
                'Collapse', 'Ambiguous']
REGIME_COLORS = {
    'Convergent':           '#2ca02c',  # green
    'Sustained oscillatory':'#ff7f0e',  # orange
    'Divergent':            '#d62728',  # red
    'Collapse':             '#9467bd',  # purple
    'Ambiguous':            '#7f7f7f',  # gray
}
OUTCOME_ORDER = ['robust', 'landing_shift', 'regime_transition',
                 'collapse', 'ambiguous']
OUTCOME_COLORS = {
    'robust':              '#2ca02c',
    'landing_shift':       '#1f77b4',
    'regime_transition':   '#d62728',
    'collapse':            '#9467bd',
    'ambiguous':           '#7f7f7f',
}
MODE_MARKERS = {
    'baseline_resim':               'o',
    'random_delete_5':              's',
    'target_strongest_5':           'D',
    'target_hub_outgoing_5':        '^',
    'target_largest_scc_internal_5':'v',
}

# Representative worlds for per-world figures (one per tier)
_REP_WORLDS = [
    'baseline_balanced_r010_t000',       # H1 stable
    'stress_balanced_r010_t003',         # H2 boundary convergent
    'stress_balanced_r010_t007',         # H2 boundary oscillatory
    'stress_balanced_r010_t005',         # H2 boundary divergent
    'stress_balanced_r010_t001',         # H2 boundary mixed
    'chen_repression_biased_r025_t007',  # H3 multi-eq
]


# ── 6.1  Regime Transition Matrix ──────────────────────────────────

# Strategy labels for display
_STRATEGY_LABELS = {
    'random_delete_5':               'A: random_delete_5',
    'target_strongest_5':            'B: target_strongest_5',
    'target_hub_outgoing_5':         'C: target_hub_outgoing_5',
    'target_largest_scc_internal_5': 'D: target_largest_scc_internal_5',
}

# World tier labels for display
_TIER_LABELS = {
    'H1_stable':   'H1: Clip-free stable anchors',
    'H2_boundary': 'H2: Boundary worlds (stress gradient)',
    'H2_control':  'H2: Other boundary controls',
    'H3_rare':     'H3: Rare specimens & outliers',
}


def _build_transition_matrix(outcomes, pids):
    """Build a pre→post regime count matrix for a set of perturbation_ids."""
    mat = pd.DataFrame(0, index=REGIME_ORDER, columns=REGIME_ORDER, dtype=int)
    for pid in pids:
        pre = outcomes[pid].get('base_primary_regime', '')
        post = outcomes[pid].get('pert_primary_regime', '')
        if pre in REGIME_ORDER and post in REGIME_ORDER:
            mat.loc[pre, post] += 1
    return mat


def _draw_one_matrix(ax, mat, title):
    """Draw a single regime transition matrix on an axis."""
    mat = mat.loc[(mat.sum(axis=1) > 0), (mat.sum(axis=0) > 0)]
    if mat.empty:
        ax.text(0.5, 0.5, 'No data', ha='center', va='center',
                transform=ax.transAxes, fontsize=10, color='gray')
        ax.set_title(title, fontsize=10)
        return
    im = ax.imshow(mat.values, cmap='Blues', aspect='auto', vmin=0)
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels(mat.columns, rotation=45, ha='right', fontsize=10)
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=10)
    ax.set_title(title, fontsize=12)
    for i in range(len(mat.index)):
        for j in range(len(mat.columns)):
            v = mat.values[i, j]
            ax.text(j, i, str(v), ha='center', va='center',
                    fontsize=10, color='white' if v > mat.values.max() * 0.5 else 'black')
    plt.colorbar(im, ax=ax, shrink=0.8)


def _fig_regime_transition_matrix(outcomes, pert_plan):
    """Two figures: (a) by deletion strategy A-D, (b) by world class tier."""
    print('  [6.1] Regime transition matrix ...')

    plan_by_pid = {p['perturbation_id']: p for p in pert_plan}

    # Build world_id → tier lookup
    wid_to_tier = {sw['world_id']: sw['tier'] for sw in SELECTED_WORLDS}

    # Partition pids by strategy and by tier
    strategy_pids: Dict[str, Set[str]] = {s: set() for s in _STRATEGY_LABELS}
    tier_pids: Dict[str, Set[str]] = {t: set() for t in _TIER_LABELS}
    for pid, oc in outcomes.items():
        mode = plan_by_pid.get(pid, {}).get('perturbation_mode', '')
        if mode in strategy_pids:
            strategy_pids[mode].add(pid)
        tier = wid_to_tier.get(oc['source_world_id'], '')
        if tier in tier_pids:
            tier_pids[tier].add(pid)

    # ── Figure (a): by strategy ──
    strategy_order = list(_STRATEGY_LABELS.keys())
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.5))
    for ax, mode in zip(axes, strategy_order):
        mat = _build_transition_matrix(outcomes, strategy_pids[mode])
        n = mat.values.sum()
        _draw_one_matrix(ax, mat,
                         f'{_STRATEGY_LABELS[mode]}\n(n={n} trajectories)')
    fig.suptitle('6.1a  Regime Transition by Deletion Strategy\n(rows = pre-perturb, columns = post-perturb)',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, 'regime_transition_by_strategy.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    → {path}')

    # ── Figure (b): by world class ──
    tier_order = list(_TIER_LABELS.keys())
    fig, axes = plt.subplots(1, 4, figsize=(22, 5.5))
    for ax, tier in zip(axes, tier_order):
        mat = _build_transition_matrix(outcomes, tier_pids[tier])
        n = mat.values.sum()
        # Show world count in title
        n_worlds = sum(1 for sw in SELECTED_WORLDS if sw['tier'] == tier)
        _draw_one_matrix(ax, mat,
                         f'{_TIER_LABELS[tier]}\n(n={n} trajectories, {n_worlds} worlds)')
    fig.suptitle('6.1b  Regime Transition by Selected World Class\n(rows = pre-perturb, columns = post-perturb)',
                 fontsize=14, y=1.02)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, 'regime_transition_by_world_class.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    → {path}')


# ── 6.2  Per-World Perturbation Outcome Heatmap ────────────────────

def _fig_per_world_heatmap(outcomes, pert_plan):
    """Heatmap: rows=worlds, cols=mode×repeat, color=dominant_outcome."""
    print('  [6.2] Per-world perturbation outcome heatmap ...')

    plan_by_pid = {p['perturbation_id']: p for p in pert_plan}

    # Group: (world_id, mode, repeat) → list of outcome labels
    groups = defaultdict(list)
    for pid, oc in outcomes.items():
        plan = plan_by_pid.get(pid, {})
        wid = oc['source_world_id']
        mode = plan.get('perturbation_mode', '')
        rpt = plan.get('repeat_index', '')
        groups[(wid, mode, rpt)].append(oc.get('perturbation_outcome', ''))

    # Compute dominant outcome per group
    _TIER_ORDER = ['H1_stable', 'H2_boundary', 'H2_control', 'H3_rare']
    _TIER_SORT = {t: i for i, t in enumerate(_TIER_ORDER)}
    _wid_to_tier_full = {sw['world_id']: sw['tier'] for sw in SELECTED_WORLDS}
    all_worlds = sorted(set(k[0] for k in groups),
                        key=lambda w: (_TIER_SORT.get(_wid_to_tier_full.get(w, ''), 99), w))
    all_modes = sorted(set(k[1] for k in groups), reverse=True)  # group by mode
    # Build column labels: mode_rpt
    col_keys = sorted(set((m, r) for _, m, r in groups),
                      key=lambda x: (all_modes.index(x[0]) if x[0] in all_modes else 99, str(x[1])))

    # Build matrix: rows=worlds, cols=mode_rpt combos
    data_raw = np.full((len(all_worlds), len(col_keys)), '', dtype=object)
    for wi, wid in enumerate(all_worlds):
        for cj, (mode, rpt) in enumerate(col_keys):
            key = (wid, mode, rpt)
            olist = groups.get(key, [])
            if olist:
                data_raw[wi, cj] = Counter(olist).most_common(1)[0][0]

    # Determine which outcomes actually appear; build compact color mapping
    _COLORBAR_ORDER = ['regime_transition', 'landing_shift', 'robust', 'collapse', 'ambiguous']
    present_set = set(o for o in data_raw.flatten() if o)
    present = [o for o in _COLORBAR_ORDER if o in present_set]
    outcome_to_idx = {o: i for i, o in enumerate(present)}
    cmap = ListedColormap([OUTCOME_COLORS[o] for o in present])
    if not present:
        present = ['robust']
        cmap = ListedColormap([OUTCOME_COLORS['robust']])
        outcome_to_idx = {'robust': 0}

    data = np.full_like(data_raw, np.nan, dtype=float)
    cell_labels = np.full_like(data_raw, '', dtype=object)
    for wi in range(len(all_worlds)):
        for cj in range(len(col_keys)):
            o = data_raw[wi, cj]
            if o:
                data[wi, cj] = outcome_to_idx[o]
                cell_labels[wi, cj] = o[:8]

    # Build col labels
    def _clean_mode(m):
        m = re.sub(r'^target_', '', m)
        m = re.sub(r'_delete', '', m)
        m = re.sub(r'_internal', '', m)
        return m
    def _format_col(m, r):
        if m.startswith('random'):
            return f'random_5\nrepeat{r}' if r else 'random_5\nrepeat0'
        return m
    col_labels = [_format_col(_clean_mode(m), r) for m, r in col_keys]
    # World labels with tier
    _TIER_SHORT = {'H1_stable': 'H1', 'H2_boundary': 'H2', 'H2_control': 'H2b', 'H3_rare': 'H3'}
    row_labels = [f'{w} [{_TIER_SHORT.get(_wid_to_tier_full.get(w, ""), "?")}]' for w in all_worlds]

    fig, ax = plt.subplots(figsize=(max(14, len(col_keys) * 1.2),
                                    max(5, len(all_worlds) * 0.45)))
    im = ax.imshow(data, cmap=cmap, aspect='auto', vmin=-0.5,
                   vmax=len(present) - 0.5)

    for wi in range(len(all_worlds)):
        for cj in range(len(col_keys)):
            if cell_labels[wi, cj]:
                ax.text(cj, wi, cell_labels[wi, cj], ha='center', va='center',
                        fontsize=5.5, color='white', fontweight='bold')

    ax.set_xticks(range(len(col_keys)))
    ax.set_xticklabels(col_labels, rotation=0, ha='center', fontsize=7)
    ax.set_yticks(range(len(all_worlds)))
    ax.set_yticklabels(row_labels, fontsize=7)
    ax.set_title('6.2  Per-World Perturbation Outcome Heatmap\n'
                 'H1: Clip-free stable anchors | H2: Boundary worlds (stress gradient) | '
                 'H2b: Other boundary controls | H3: Rare specimens & outliers',
                 fontsize=10)

    # Colorbar: only show outcomes that actually appear
    cbar = plt.colorbar(im, ax=ax, ticks=range(len(present)), shrink=0.8, pad=0.02)
    cbar.ax.set_yticklabels(present, fontsize=8)

    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, 'per_world_perturbation_outcome_heatmap.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    → {path}')


# ── 6.3  State Landing Plots (2×2 grid: H1/H2/H2b/H3) ────────────

# Marker / short-name for each perturbation mode
_MODE_MARKER = {
    'baseline':                    'o',
    'random_delete_5':             's',
    'target_strongest_5':          '^',
    'target_hub_outgoing_5':       'D',
    'target_largest_scc_internal_5': 'v',
}
_MODE_SHORT = {
    'baseline':                       'baseline',
    'random_delete_5':                'random',
    'target_strongest_5':             'strongest',
    'target_hub_outgoing_5':          'hub',
    'target_largest_scc_internal_5':  'SCC',
}
_TIER_PICK = {
    'H1': 'baseline_balanced_r010_t000',
    'H2': 'stress_balanced_r010_t003',
    'H2b': 'stress_balanced_r010_t007',
    'H3': 'chen_repression_biased_r025_t007',
}


def _fig_state_landing(outcomes, pert_plan):
    """Single 2×2 figure: one panel per representative world (per Output Spec §6.3).
    All post-perturb final states are plotted (convergent + non-convergent).
    Color = state_landing_id, marker = perturbation_mode.
    Each panel runs its own within-world PCA (cross-world PCA forbidden by spec).
    """
    print('  [6.3] State landing plots (2×2: H1/H2/H2b/H3, color = lid, marker = mode) ...')

    plan_by_pid = {p['perturbation_id']: p for p in pert_plan}
    strategies = ['random_delete_5', 'target_strongest_5',
                  'target_hub_outgoing_5', 'target_largest_scc_internal_5']
    strat_short = {
        'random_delete_5':              'random',
        'target_strongest_5':           'strongest',
        'target_hub_outgoing_5':        'hub',
        'target_largest_scc_internal_5': 'largest_scc',
    }

    def _base_mode(m):
        for s in strategies:
            if m == s or m.startswith(s + '_r'):
                return s
        return m

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes_flat = axes.flatten()
    panel_titles = [
        'H1: clip-free stable',
        'H2: boundary (stress)',
        'H2b: boundary (other)',
        'H3: rare specimens',
    ]

    for ax, (tier, wid), ptitle, legend_loc in zip(
            axes_flat, _TIER_PICK.items(), panel_titles,
            ['upper right', 'upper right', 'lower right', 'lower right']):
        pts, lids, modes = [], [], []
        for pid, oc in outcomes.items():
            if oc['source_world_id'] != wid:
                continue
            mode = plan_by_pid.get(pid, {}).get('perturbation_mode', '')
            base = _base_mode(mode)
            if base not in strategies:
                continue
            pts.append(oc['pert_equilibrium'])
            lids.append(oc.get('state_landing_id_post'))
            modes.append(base)

        if len(pts) < 3:
            ax.text(0.5, 0.5, f'{wid}\n(no trajectories)', ha='center', va='center',
                    transform=ax.transAxes, fontsize=10, color='gray')
            ax.set_xticks([]); ax.set_yticks([])
            ax.set_title(f'{tier} — {ptitle}', fontsize=11)
            continue

        X = np.stack(pts)
        pca = PCA(n_components=2)
        X2 = pca.fit_transform(X)
        # Two-pass drawing: no-landing first (bottom layer, gray), then landing points (top, colored)
        idxs_no = [i for i, lid in enumerate(lids) if not (isinstance(lid, int) and lid >= 0)]
        idxs_ok = [i for i, lid in enumerate(lids) if isinstance(lid, int) and lid >= 0]
        cmap = plt.cm.tab10
        # Pass 1: no-landing, gray, behind everything
        for i in idxs_no:
            m = modes[i]
            marker = _MODE_MARKER.get(m, 'o')
            ax.scatter(X2[i, 0], X2[i, 1], c=['gray'], marker=marker,
                       s=25, alpha=0.5, edgecolors='none')
        # Pass 2: landing points on top
        seen_pairs = set()
        combo_counts: Dict[tuple, int] = {}
        for i in idxs_ok:
            lid = lids[i]
            l_text = f'L{lid}'
            key = (strat_short[modes[i]], l_text)
            combo_counts[key] = combo_counts.get(key, 0) + 1
        for i in idxs_ok:
            lid = lids[i]
            color = cmap(lid % 10)
            marker = _MODE_MARKER.get(modes[i], 'o')
            l_text = f'L{lid}'
            pair = (strat_short[modes[i]], l_text)
            cnt = combo_counts[pair]
            if pair not in seen_pairs:
                label = f'{strat_short[modes[i]]} · {l_text} (n={cnt})'
                seen_pairs.add(pair)
            else:
                label = None
            ax.scatter(X2[i, 0], X2[i, 1], c=[color], marker=marker,
                       s=25, alpha=0.7, edgecolors='k', linewidths=0.3,
                       label=label)
        # Add no-landing legend entry once
        if idxs_no:
            ax.scatter([], [], c=['gray'], marker='s', s=25, alpha=0.5,
                       label=f'no landing (n={len(idxs_no)})')

        ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%})')
        ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%})')
        ax.set_title(f'{tier} — {ptitle}\n{wid} (n={X.shape[0]})', fontsize=10)
        ax.legend(fontsize=7, loc=legend_loc, framealpha=0.85, ncol=2)
        ax.grid(True, alpha=0.3)

    fig.suptitle('6.3  State Landing — one world per tier '
                 '(color = state_landing_id, marker = perturbation_mode)',
                 fontsize=12, y=1.00)
    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, 'state_landing_2x2.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    → {path}')


# ── 6.4  Pre/Post Trajectory Plots (per representative world) ──────

def _fig_trajectory_pre_post(outcomes, fork_results, pert_plan):
    """For each tier-representative world (H1/H2/H2b/H3), pick the first rep,
    then choose 4 perturbations (one per deletion strategy A/B/C/D).
    Each world -> one figure: rows = deletion strategies, cols = (full / zoom / late).
    Baseline is shown as a faint reference background on every panel.
    """
    print('  [6.4] Pre/post trajectory plots (4 worlds × 4 strategies × 3 views) ...')

    plan_by_pid = {p['perturbation_id']: p for p in pert_plan}
    strategies = ['random_delete_5', 'target_strongest_5',
                  'target_hub_outgoing_5', 'target_largest_scc_internal_5']
    strat_short = {
        'random_delete_5':              'A: random',
        'target_strongest_5':           'B: strongest',
        'target_hub_outgoing_5':        'C: hub',
        'target_largest_scc_internal_5': 'D: largest_scc',
    }

    for tier, wid in _TIER_PICK.items():
        # Pick the first rep's pid for each of the 4 strategies
        # A random -> first "random_delete_5" entry (not _rN variants)
        # B/C/D -> first pid in each mode
        mode_pids: Dict[str, str] = {}
        for pid, oc in outcomes.items():
            if oc['source_world_id'] != wid:
                continue
            mode = plan_by_pid.get(pid, {}).get('perturbation_mode', '')
            # for random we want only the base "random_delete_5", not _r1..r9
            if mode == 'random_delete_5':
                if 'random_delete_5' not in mode_pids:
                    mode_pids['random_delete_5'] = pid
            elif mode in ('target_strongest_5', 'target_hub_outgoing_5',
                          'target_largest_scc_internal_5'):
                if mode not in mode_pids:
                    mode_pids[mode] = pid

        # Need the baseline trajectory of the same rep used for the chosen pids.
        # All chosen pids should share the same X0 (rep), so any non-baseline pid's
        # baseline_traj in fork_results is the corresponding reference.
        # We pick whichever pid exists first and use ITS baseline_traj.
        ref_pid = next((pid for pid in mode_pids.values()), None)
        if ref_pid is None or ref_pid not in fork_results:
            print(f'    [{wid}]  no trajectories found, skipping')
            continue

        # Build 4-row x 3-col figure
        fig, axes = plt.subplots(4, 3, figsize=(15, 12),
                                  gridspec_kw={'hspace': 0.85, 'wspace': 0.3})

        for row, mode in enumerate(strategies):
            pid = mode_pids.get(mode)
            if pid is None or pid not in fork_results:
                for col in range(3):
                    axes[row, col].text(0.5, 0.5, 'no data',
                                          ha='center', va='center',
                                          transform=axes[row, col].transAxes,
                                          fontsize=8, color='gray')
                continue
            fr = fork_results[pid]
            base_traj = fr['baseline_traj']   # (T+1, G)
            pert_traj = fr['perturbed_traj']  # (T+1, G)
            n_genes = base_traj.shape[1]
            t_full = np.arange(base_traj.shape[0])

            oc = outcomes.get(pid, {})
            outcome = oc.get('perturbation_outcome', 'unknown')
            pre_regime  = oc.get('base_primary_regime', '?')
            post_regime = oc.get('pert_primary_regime', '?')
            regime_str  = f'\n({pre_regime} → {post_regime})'
            title_prefix = f"{strat_short[mode]} → {outcome}{regime_str}"

            def draw_panel(ax, t_slice, label_extra=''):
                # baseline gene trajectories as faint dashed blue (only post-perturbation)
                t_in_slice_local = t_slice
                b_in = (t_in_slice_local >= T_FORK)
                if b_in.any():
                    for g in range(n_genes):
                        ax.plot(t_in_slice_local[b_in], base_traj[t_in_slice_local[b_in], g],
                                color='#1f77b4', alpha=0.3, linewidth=0.5,
                                linestyle='--')
                # perturbed gene trajectories (full t range covered by t_slice)
                for g in range(n_genes):
                    ax.plot(t_in_slice_local, pert_traj[t_in_slice_local, g],
                            color='#d62728', alpha=0.55, linewidth=0.5)
                ax.axvline(T_FORK, color='gray', linestyle='--',
                           linewidth=1.0, alpha=0.7)
                ax.set_xlim(t_slice[0], t_slice[-1])

            # 1) Full trajectory
            ax = axes[row, 0]
            draw_panel(ax, t_full)
            ax.set_title(f'{title_prefix}\nfull trajectory ({n_genes} genes)',
                         fontsize=11)
            ax.set_xlabel('step')
            ax.set_ylabel('X (abundance)')

            # 2) Zoom around t_perturb
            ax = axes[row, 1]
            z0, z1 = max(0, T_FORK - 50), min(T_TOTAL, T_FORK + 50)
            tz = t_full[z0:z1 + 1]
            draw_panel(ax, tz)
            ax.set_title(f'{title_prefix}\nzoom ±50 steps around t_perturb={T_FORK}', fontsize=11)
            ax.set_xlabel('step')

            # 3) Late-time window
            ax = axes[row, 2]
            lz = np.arange(T_TOTAL - LATE_WINDOW, T_TOTAL + 1)
            draw_panel(ax, lz)
            ax.set_title(f'{title_prefix}\nlate-time window (last {LATE_WINDOW} steps)', fontsize=11)
            ax.set_xlabel('step')

        # Legend only once in upper-left of the first panel (row=0, col=0)
        # (use fixed G=20, the canonical gene count)
        axes[0, 0].plot([], [], color='#1f77b4', linestyle='--',
                         linewidth=1.0, label='baseline (post-perturb only, dashed)')
        axes[0, 0].plot([], [], color='#d62728', linewidth=1.0,
                         label='perturbed (all 20 genes)')
        axes[0, 0].legend(fontsize=6, loc='best', framealpha=0.85)

        fig.suptitle(f'6.4  Pre/Post Trajectory — {wid}  [{tier}]\n'
                     f'(rows: deletion strategy A/B/C/D;\n'
                     f'cols: full / zoom / late; faint blue = baseline reference)',
                     fontsize=14, y=1.02)
        plt.tight_layout()
        path = os.path.join(FIGURES_DIR,
                            f'trajectory_pre_post_world_{wid}.png')
        fig.savefig(path, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f'    → {path}')


# ── 6.5  Delta rho(W) vs Outcome ───────────────────────────────────

def _fig_delta_rho_vs_outcome(outcomes, pert_plan):
    """Scatter: delta_rho_W vs pre_post_normalized_L2.

    Design rationale (v2):
      - log y-axis to handle 10+ orders of magnitude naturally (no broken axis).
      - Color = perturbation_outcome (the run_006 primary classification).
      - Shape = perturbation_mode (A/B/C/D).
      - Alpha transparency + slight x-jitter to relieve overplotting.
      - Reference lines at LANDING_SAME_THRESHOLD (robust landing_shift boundary)
        and y=1.0 (divergent saturation artifact).
      - Marginal count annotation per outcome category.
    """
    print('  [6.5] Delta rho(W) vs outcome ...')

    plan_by_pid = {p['perturbation_id']: p for p in pert_plan}

    # ── gather data ──
    Y_FLOOR_LOG = 1e-4  # floor for log-scale visibility
    pts = []
    n_divergent_excluded = 0
    for pid, oc in outcomes.items():
        plan = plan_by_pid.get(pid, {})
        dr = plan.get('delta_rho_W')
        l2 = oc.get('pre_post_normalized_L2')
        mode = plan.get('perturbation_mode', '')
        outcome = oc.get('perturbation_outcome', '')
        base_regime = oc.get('base_primary_regime', '')
        post_regime = oc.get('pert_primary_regime', '')
        if dr is None or l2 is None or l2 == '' or mode not in MODE_MARKERS:
            continue
        # Skip divergent trajectories: normalized L2 is not meaningful for them
        # (doc §13.1: "do not interpret final point as equilibrium; compare
        # boundedness/growth-rate/magnitude summary instead")
        if base_regime == 'Divergent' or post_regime == 'Divergent':
            n_divergent_excluded += 1
            continue
        val = max(float(l2), Y_FLOOR_LOG) if float(l2) > 0 else Y_FLOOR_LOG
        pts.append({
            'dr': float(dr), 'l2': val,
            'mode': mode, 'outcome': outcome,
        })

    # ── figure layout ──
    fig, ax = plt.subplots(figsize=(11, 7))

    MODE_SHORT = {
        'random_delete_5':               'A: random',
        'target_strongest_5':            'B: strongest',
        'target_hub_outgoing_5':         'C: hub',
        'target_largest_scc_internal_5': 'D: SCC',
    }

    rng = np.random.default_rng(42)
    JITTER_SCALE = 0.003  # x-axis jitter

    for mode, marker in MODE_MARKERS.items():
        if mode == 'baseline_resim':
            continue
        sub = [p for p in pts if p['mode'] == mode]
        if not sub:
            continue
        xs = [p['dr'] + rng.uniform(-JITTER_SCALE, JITTER_SCALE) for p in sub]
        ys = [p['l2'] for p in sub]
        cs = [OUTCOME_COLORS.get(p['outcome'], '#7f7f7f') for p in sub]
        ax.scatter(xs, ys, c=cs, marker=marker, s=40, alpha=0.55,
                   edgecolors='none',
                   label=MODE_SHORT.get(mode, mode))

    # ── reference lines ──
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.4)

    ax.set_xlabel(r'$\Delta \rho(W)$  (negative → spectral radius decreased)',
                  fontsize=11)
    ax.set_ylabel('normalized L2 distance  (linear)', fontsize=11)
    ax.set_title('Delta ρ(W) vs Perturbation Outcome' + '\n(excluding trajectories with Divergent baseline or perturbed regime)',
                 fontsize=13, pad=12)

    ax.grid(True, alpha=0.25)

    # ── xlim from data ──
    all_x = [p['dr'] for p in pts]
    x_margin = 0.02
    ax.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)

    # ── marginal count annotation (top-right inset) ──
    outcome_counts = Counter(p['outcome'] for p in pts)
    count_lines = ['Outcome  |  N']
    count_lines.append('─' * 18)
    for label in OUTCOME_ORDER:
        n = outcome_counts.get(label, 0)
        pct = n / len(pts) * 100 if pts else 0
        count_lines.append(f'{label:18s} {n:4d}  ({pct:4.1f}%)')
    count_text = '\n'.join(count_lines)
    ax.text(1.01, 0.35, count_text, transform=ax.transAxes,
            fontsize=7.5, fontfamily='monospace', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='#cccccc', alpha=0.85),
            clip_on=False)

    # ── single combined legend (only show categories present in data) ──
    from matplotlib.lines import Line2D
    present_outcomes = set(p['outcome'] for p in pts)
    present_modes = set(p['mode'] for p in pts)

    combined_handles = []
    if present_outcomes:
        combined_handles.append(Line2D([], [], linestyle='none', label='═ Color: outcome'))
        for label in OUTCOME_ORDER:
            if label in present_outcomes:
                c = OUTCOME_COLORS[label]
                combined_handles.append(
                    Line2D([0], [0], marker='s', color='w', markerfacecolor=c,
                           markersize=9, label=label))
    if present_modes:
        combined_handles.append(Line2D([], [], linestyle='none', label=''))
        combined_handles.append(Line2D([], [], linestyle='none', label='═ Shape: perturbation'))
        for mode_key, mode_label in MODE_SHORT.items():
            if mode_key in present_modes:
                combined_handles.append(
                    Line2D([0], [0], marker=MODE_MARKERS[mode_key], color='gray',
                           markerfacecolor='gray', markersize=8, label=mode_label))

    ax.legend(handles=combined_handles, fontsize=7.5,
              loc='upper left', bbox_to_anchor=(1.01, 1),
              title='Legend', framealpha=0.85, borderaxespad=0)

    path = os.path.join(FIGURES_DIR, 'delta_rhoW_vs_outcome.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    → {path}')


# ── 6.5b  Delta rho(W) vs Divergent Magnitude Ratio ───────────────

def _fig_divergent_magnitude_vs_outcome(outcomes, pert_plan):
    """Scatter: delta_rho_W vs magnitude_ratio for divergent-involved trajectories.

    For trajectories where baseline or perturbed regime is Divergent,
    state-landing L2 distance is not meaningful. Instead we compute:
        magnitude_ratio = ||post_final|| / max(||pre_reference||, eps)

    ratio < 1 → deletion suppressed divergence
    ratio > 1 → deletion accelerated divergence
    ratio ≈ 1 → no effect on divergence magnitude
    """
    print('  [6.5b] Delta rho(W) vs divergent magnitude ratio ...')

    plan_by_pid = {p['perturbation_id']: p for p in pert_plan}

    # ── gather data ──
    pts = []
    n_total = 0
    transition_labels = []
    for pid, oc in outcomes.items():
        plan = plan_by_pid.get(pid, {})
        dr = plan.get('delta_rho_W')
        mode = plan.get('perturbation_mode', '')
        if dr is None or mode not in MODE_MARKERS:
            continue

        base_regime = oc.get('base_primary_regime', '')
        post_regime = oc.get('pert_primary_regime', '')

        # Only include divergent-involved trajectories
        if base_regime != 'Divergent' and post_regime != 'Divergent':
            continue

        n_total += 1
        pre_norm = oc.get('pre_reference_state_norm', 0)
        post_norm = oc.get('post_final_state_norm', 0)
        if pre_norm is None or post_norm is None:
            continue

        eps_norm = 1e-8
        ratio = post_norm / max(pre_norm, eps_norm)

        # Classify transition type
        if base_regime == 'Divergent' and post_regime == 'Divergent':
            if ratio < 0.9:
                trans_type = 'Div→Div (suppressed)'
            elif ratio <= 1.1:
                trans_type = 'Div→Div (unchanged)'
            else:
                trans_type = 'Div→Div (accelerated)'
        elif base_regime != 'Divergent' and post_regime == 'Divergent':
            trans_type = f'{base_regime[:4]}→Div'
        else:
            trans_type = f'{base_regime}→{post_regime}'

        pts.append({
            'dr': float(dr), 'ratio': ratio,
            'mode': mode, 'trans_type': trans_type,
        })
        transition_labels.append(trans_type)

    # unique transition types — define colour palette
    present_types = sorted(set(transition_labels), reverse=True)
    type_colors = {
        'Div→Div (accelerated)': '#d62728',
        'Conv→Div':               '#ff7f0e',
        'Osc→Div':                '#e377c2',
        'Div→Div (unchanged)':    '#7f7f7f',
        'Div→Div (suppressed)':   '#2ca02c',
    }
    # only include those present
    type_colors = {k: v for k, v in type_colors.items() if k in present_types}

    MODE_SHORT = {
        'random_delete_5':               'A: random',
        'target_strongest_5':            'B: strongest',
        'target_hub_outgoing_5':         'C: hub',
        'target_largest_scc_internal_5': 'D: SCC',
    }

    # ── plot ──
    fig, ax = plt.subplots(figsize=(11, 7))

    rng = np.random.default_rng(42)
    JITTER_SCALE = 0.003

    for mode, marker in MODE_MARKERS.items():
        if mode == 'baseline_resim':
            continue
        sub = [p for p in pts if p['mode'] == mode]
        if not sub:
            continue
        xs = [p['dr'] + rng.uniform(-JITTER_SCALE, JITTER_SCALE) for p in sub]
        ys = [p['ratio'] for p in sub]
        cs = [type_colors.get(p['trans_type'], '#7f7f7f') for p in sub]
        ax.scatter(xs, ys, c=cs, marker=marker, s=40, alpha=0.55,
                   edgecolors='none',
                   label=MODE_SHORT.get(mode, mode))

    ax.set_yscale('log')
    ax.axhline(1.0, color='gray', linewidth=0.8, linestyle='--', alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.5, alpha=0.4)

    ax.set_xlabel(r'$\Delta \rho(W)$  (negative → spectral radius decreased)',
                  fontsize=11)
    ax.set_ylabel('magnitude ratio  (post / pre, log scale)', fontsize=11)
    ax.set_title(r'$\Delta \rho(W)$ vs Divergent Magnitude Ratio', fontsize=13)

    ax.grid(True, alpha=0.25, which='both')

    # xlim from data
    all_x = [p['dr'] for p in pts]
    x_margin = 0.02
    ax.set_xlim(min(all_x) - x_margin, max(all_x) + x_margin)

    # Annotate counts per transition type
    from collections import Counter
    tcounts = Counter(p['trans_type'] for p in pts)
    count_lines = ['Transition type        |  N']
    count_lines.append('─' * 30)
    for ttype in ['Div→Div (accelerated)', 'Conv→Div', 'Osc→Div',
                  'Div→Div (unchanged)', 'Div→Div (suppressed)']:
        n = tcounts.get(ttype, 0)
        if n > 0:
            pct = n / len(pts) * 100 if pts else 0
            count_lines.append(f'{ttype:30s} {n:4d}  ({pct:4.1f}%)')
    count_text = '\n'.join(count_lines)
    ax.text(1.01, 0.98, count_text, transform=ax.transAxes,
            fontsize=7.5, fontfamily='monospace', va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white',
                      edgecolor='#cccccc', alpha=0.85),
            clip_on=False)

    # ── combined legend (color + shape) ──
    from matplotlib.lines import Line2D
    combined_handles = []
    if present_types:
        combined_handles.append(Line2D([], [], linestyle='none', label='═ Color: transition type'))
        for ttype in ['Div→Div (accelerated)', 'Conv→Div', 'Osc→Div',
                      'Div→Div (unchanged)', 'Div→Div (suppressed)']:
            if ttype in present_types:
                c = type_colors[ttype]
                combined_handles.append(
                    Line2D([0], [0], marker='s', color='w', markerfacecolor=c,
                           markersize=9, label=ttype))
    combined_handles.append(Line2D([], [], linestyle='none', label=''))
    combined_handles.append(Line2D([], [], linestyle='none', label='═ Shape: perturbation'))
    for mode_key in ['random_delete_5', 'target_strongest_5',
                     'target_hub_outgoing_5', 'target_largest_scc_internal_5']:
        combined_handles.append(
            Line2D([0], [0], marker=MODE_MARKERS[mode_key], color='gray',
                   markerfacecolor='gray', markersize=8, label=MODE_SHORT[mode_key]))
    ax.legend(handles=combined_handles, fontsize=7.5,
              loc='upper left', bbox_to_anchor=(1.01, 0.28),
              title='Legend', framealpha=0.85, borderaxespad=0)

    plt.tight_layout()
    path = os.path.join(FIGURES_DIR, 'delta_rhoW_vs_divergent_magnitude.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'    → {path}  (n_total={n_total}, plotted={len(pts)})')


# ═══════════════════════════════════════════════════════════════════
# §9 Dispatcher
# ═══════════════════════════════════════════════════════════════════

def s9_make_figures(outcomes, state_landing_results, fork_results, pert_plan):
    """Generate all required figures (§6.1–6.5) plus §6.5b (divergent magnitude)."""
    print(f'\n{"="*70}')
    print(f'§9  Figure Generation')
    print(f'{"="*70}')

    _fig_regime_transition_matrix(outcomes, pert_plan)
    _fig_per_world_heatmap(outcomes, pert_plan)
    _fig_state_landing(outcomes, pert_plan)
    _fig_trajectory_pre_post(outcomes, fork_results, pert_plan)
    _fig_delta_rho_vs_outcome(outcomes, pert_plan)
    _fig_divergent_magnitude_vs_outcome(outcomes, pert_plan)

    print(f'  DONE: 6 figures saved to {FIGURES_DIR}/')


# ═══════════════════════════════════════════════════════════════════
# Main
# ═══════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()
    print('=' * 70)
    print('Run 006 — GRN Edge Perturbation Sensitivity Experiment')
    print('  Execution grid: core')
    print(f'  {len(SELECTED_IDS)} selected worlds from run_005')
    print(f'  Fork design: t_fork={T_FORK}, T_post={T_POST}')
    print(f'  K_delete={K_DELETE}, Strategies: A(random×{R_RANDOM}) B(strongest) C(hub) D(scc)')
    print('=' * 70)

    # ── §1  Data loading & preparation ──
    selected_worlds = s1_load_and_prepare()

    # ── §2  Initial-condition screening ──
    screening = s2_screen_init_conditions(selected_worlds)

    # ── §3  Edge deletion strategies ──
    pert_plan = s3_generate_perturbations(selected_worlds, screening)

    # ── §4  Fork simulation ──
    fork_results = s4_run_fork_simulations(selected_worlds, screening, pert_plan)

    # ── §5  Regime classification ──
    traj_results = s5_classify_regimes(fork_results)

    # ── §6  Outcome assignment ──
    outcomes = s6_assign_outcomes(traj_results, fork_results)

    # ── §7  State landing clustering ──
    state_landing_results = s7_state_landing_clustering(outcomes, fork_results, pert_plan)

    # ── §8  Tables ──
    s8_build_tables(screening, pert_plan, outcomes, state_landing_results,
                    selected_worlds)

    # ── §9  Figures ──
    s9_make_figures(outcomes, state_landing_results, fork_results, pert_plan)

    elapsed = time.time() - t0
    print(f'\n{"="*70}')
    print(f'Run 006 complete — elapsed {elapsed/60:.1f} min ({elapsed:.1f} s)')
    print(f'{"="*70}')


if __name__ == '__main__':
    main()
