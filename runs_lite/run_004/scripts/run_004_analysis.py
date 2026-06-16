"""
Run 004 — Model A Decay / Memory Boundary
===========================================

Per: docs/01_DDC_Lite_Curriculum/run_004/

Full Grid:
  2 topology_types × 3 δ × 4 strength_regimes × 2 sign_ratios = 48 combos
  N_world = 10, N_init = 5, T = 1000
  → 480 worlds, 2400 trajectories

Steps:
  1.  Full Grid config generation
  2.  World sampling + trajectory simulation
  3.  Per-trajectory analysis (equilibrium / stability / oscillation / classify)
  4.  Aggregation (world-level → combo-level → cross-dimension)
  5.  Tables (world_summary / regime_summary / boundary_summary)
  6.  Figures (heatmaps / topology×δ / clipping / ρ(A) / trajectory exemplars)
  7.  Perturbation (canonical KD at t=500)
  8.  Summary (Q1–Q7 answers)

Author: zhanghl
"""
import json
import os
import sys
import csv
import copy
import torch
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Any

# Ensure src/ is on the path for ddc_model_a / ddc_model_a_tf imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'src'))


# ═══════════════════════════════════════════════════════════════════
# Step 1:  Configuration & Constants
# ═══════════════════════════════════════════════════════════════════

# ── Run identity ──
RUN_ID: str = 'run_004'
EXECUTION_GRID: str = 'full'
RUNTIME_VERSION: str = 'v0.3-delta'

# ── Model constants (matched to ddc_model_a / ddc_model_a_tf) ──
G: int = 20
N_TF: int = 5
TF_GENES: List[int] = [0, 1, 2, 3, 4]
N_EDGES: int = 38
EDGE_DENSITY_R: float = 0.1
T_SIM: int = 1000
DEFAULT_B: float = 0.1

# ── Topology ──
TOPOLOGY_TYPES: List[str] = [
    'unrestricted_sparse',
    'tf_restricted',
]

# ── δ scan ──
DELTA_VALUES: List[float] = [
    0.1,   # self-retention = 0.9  (high memory / low decay)
    0.2,   # self-retention = 0.8  (reference: run_002 / run_003)
    0.4,   # self-retention = 0.6  (low memory / high decay)
]

# ──  delta label helper (0.1 → d0p1, 0.2 → d0p2, 0.4 → d0p4) ──
def delta_label(d: float) -> str:
    """Convert delta float to label: 0.1 → d0p1, 0.2 → d0p2, 0.4 → d0p4."""
    #  Avoid floating-point rounding: use fixed mapping.
    mapping = {0.1: 'd0p1', 0.2: 'd0p2', 0.4: 'd0p4'}
    return mapping[d]


# ── Interaction strength regimes (per Analysis Plan §4.3) ──
STRENGTH_REGIMES: Dict[str, Tuple[float, float, str]] = {
    'baseline':        (0.02, 0.10, 'conservative baseline'),
    'chen_moderate':   (0.10, 0.20, 'Chen-moderate'),
    'stress':          (0.10, 0.30, 'stress-test'),
    'chen_stress':     (0.20, 0.40, 'Chen-stress'),
}
STRENGTH_ORDER: List[str] = ['baseline', 'chen_moderate', 'stress', 'chen_stress']
STRENGTH_INDEX: Dict[str, int] = {k: i for i, k in enumerate(STRENGTH_ORDER)}

# ── Sign ratios (per Analysis Plan §4.4) ──
SIGN_RATIOS: Dict[str, float] = {
    'balanced':          0.5,     # activation : repression = 1 : 1
    'repression_biased': 0.333,   # activation : repression = 1 : 2
}
SIGN_RATIO_LABELS: Dict[str, str] = {
    'balanced':          'Balanced',
    'repression_biased': 'Repression-biased',
}
SIGN_RATIO_RATIOS: Dict[str, str] = {
    'balanced':          '1:1',
    'repression_biased': '1:2',
}

# ── Experimental scale ──
N_WORLD: int = 10       # worlds per combo
N_INIT:  int = 5        # initial states per world

# ── Seed allocation ──
SIGN_BASE: Dict[str, int] = {
    'balanced':          2000,
    'repression_biased': 6000,
}
DELTA_OFFSET: Dict[float, int] = {
    0.1:      0,
    0.2:  10000,
    0.4:  20000,
}
TOPOLOGY_OFFSET: Dict[str, int] = {
    'unrestricted_sparse':     0,
    'tf_restricted':        5000,
}
# topo_seed  = SIGN_BASE[sign] + DELTA_OFFSET[δ] + TOPOLOGY_OFFSET[topo] + ti × 100
# world_seed = topo_seed + 50 + STRENGTH_INDEX[strength] × 10
# cell_seed  = world_seed + 1 + ci                                           (ci ∈ 0..4)

# ── Thresholds (matched to run_002 / run_003) ──
EPSILON:          float = 1e-4    # convergence: ‖X(t+1)-X(t)‖ < ε
CONVERGENCE_WINDOW: int = 50      # consecutive steps for convergence
COLLAPSE_THRESHOLD: float = 1e-3  # final mean expression < this → Collapse
DIVERGENCE_THRESHOLD: float = 1e3 # max expression > this → Divergent
CLIPPING_FRAC_THRESHOLD: float = 0.1  # clipping_freq > 0.1 → clipping_dominated
SLOW_CONVERGENCE_THRESHOLD: int = 200  # fast ≤ 200 < slow
SPARSITY_THRESHOLD: float = 1e-3      # x_eq_i < this considered near-zero

# ── Perturbation ──
PERTURBATION_TIME: int = 500   # intervention step

# ── Output directories ──
RUN_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS_DIR = os.path.join(RUN_DIR, 'scripts')
FIGURES_DIR = os.path.join(RUN_DIR, 'figures')
TABLES_DIR = os.path.join(RUN_DIR, 'tables')
TRAJECTORIES_DIR = os.path.join(RUN_DIR, 'trajectories')
METADATA_DIR = os.path.join(RUN_DIR, 'world_metadata')
ANALYSIS_DIR = os.path.join(RUN_DIR, 'analysis')
SUMMARY_DIR = os.path.join(RUN_DIR, 'summary')
PERTURBATIONS_DIR = os.path.join(RUN_DIR, 'perturbations')

# ── Primary regime labels (per Problem Definition §7) ──
PRIMARY_REGIME_ORDER = [
    'Convergent',
    'Sustained oscillatory',
    'Divergent',
    'Collapse',
    'Ambiguous',
]
PRIMARY_REGIME_COLORS = {
    'Convergent':           '#2ecc71',
    'Sustained oscillatory': '#e74c3c',
    'Divergent':            '#9b59b6',
    'Collapse':             '#95a5a6',
    'Ambiguous':            '#bdc3c7',
}
PRIMARY_REGIME_DESCRIPTIONS = {
    'Convergent':           'converges to equilibrium',
    'Sustained oscillatory': 'sustained bounded oscillation',
    'Divergent':            'runaway divergence',
    'Collapse':             'numerical collapse',
    'Ambiguous':            'ambiguous or numerical failure',
}

# ── Secondary labels ──
SECONDARY_LABELS = [
    'fast_convergence',
    'slow_convergence',
    'damped_oscillatory_transient',
    'clipping_dominated',
    'low_clipping_linear_like',
    'runaway_divergence',
    'numerical_collapse',
]

# ── Legacy Type → Primary / Secondary mapping ──
LEGACY_TO_PRIMARY: Dict[str, Tuple[str, List[str]]] = {
    'Type A': ('Convergent', ['fast_convergence']),
    'Type B': ('Convergent', ['slow_convergence']),
    'Type C': ('Convergent', ['damped_oscillatory_transient']),
    'Type D': ('Sustained oscillatory', []),
    'Type E': ('Divergent', ['runaway_divergence']),
    'Type F': ('Collapse', ['numerical_collapse']),
    'Type G': ('Ambiguous', []),
}

# ── Attractor names (for legacy compatibility in plots) ──
ATTRACTOR_NAMES = {
    'Type A': 'Stable equilibrium',
    'Type B': 'Slow convergence',
    'Type C': 'Damped oscillation',
    'Type D': 'Sustained oscillation',
    'Type E': 'Runaway divergence',
    'Type F': 'Numerical collapse',
    'Type G': 'Others',
}

# ── Topology colors (for TF / non-TF distinction) ──
TF_COLOR = 'red'
NON_TF_COLOR = 'grey'

# ── Sign ratio colors ──
SIGN_COLORS = {
    'balanced':          '#3498db',
    'repression_biased': '#9b59b6',
}

# ── δ colors ──
DELTA_COLORS = {
    0.1: '#e74c3c',   # red   – high memory
    0.2: '#f39c12',   # orange – reference
    0.4: '#3498db',   # blue  – low memory
}

# ═══════════════════════════════════════════════════════════════════
# Step 1:  Full Grid Config Generation
# ═══════════════════════════════════════════════════════════════════

def generate_all_combos() -> List[Dict[str, Any]]:
    """Generate the Full Grid experimental config: 48 combos.

    Per Analysis Plan §5.2:
      2 topology_types × 3 δ × 4 strength_regimes × 2 sign_ratios = 48 combos

    Each combo dict contains all parameters needed to drive
    world sampling and simulation.
    """
    combos: List[Dict[str, Any]] = []

    for topology_type in TOPOLOGY_TYPES:
        for delta_val in DELTA_VALUES:
            dlabel = delta_label(delta_val)
            for strength_regime, (a_min, a_max, regime_label) in STRENGTH_REGIMES.items():
                for sign_ratio_key, sign_ratio_val in SIGN_RATIOS.items():
                    combo_key = (
                        f"{topology_type}_{dlabel}"
                        f"_{strength_regime}_{sign_ratio_key}"
                    )
                    combo = {
                        'combo_key':        combo_key,
                        'topology_type':    topology_type,
                        'delta':            delta_val,
                        'delta_label':      dlabel,
                        'self_retention':   1.0 - delta_val,
                        'strength_regime':  strength_regime,
                        'a_min':            a_min,
                        'a_max':            a_max,
                        'regime_label':     regime_label,
                        'sign_ratio':       sign_ratio_key,
                        'sign_ratio_value': sign_ratio_val,
                        'strength_idx':     STRENGTH_INDEX[strength_regime],
                        'sign_base':        SIGN_BASE[sign_ratio_key],
                        'delta_offset':     DELTA_OFFSET[delta_val],
                        'topology_offset':  TOPOLOGY_OFFSET[topology_type],
                    }
                    combos.append(combo)

    return combos


def make_world_id(topology_type: str, delta_label_str: str,
                  strength_regime: str, sign_ratio: str,
                  topo_idx: int) -> str:
    """Construct world_id per Analysis Plan §6 naming convention.

    E.g.: unrestricted_d0p1_baseline_balanced_t000
          tf_d0p2_chen_stress_repression_biased_t009
    """
    return f"{topology_type}_{delta_label_str}_{strength_regime}_{sign_ratio}_t{topo_idx:03d}"


def compute_topology_seed(combo: Dict[str, Any], topo_idx: int) -> int:
    """Compute deterministic topology seed."""
    return (combo['sign_base'] + combo['delta_offset']
            + combo['topology_offset'] + topo_idx * 100)


def compute_world_seed(topo_seed: int, strength_idx: int) -> int:
    """Compute deterministic world seed."""
    return topo_seed + 50 + strength_idx * 10


def compute_cell_seed(world_seed: int, cell_idx: int) -> int:
    """Compute deterministic cell seed."""
    return world_seed + 1 + cell_idx


# ═══════════════════════════════════════════════════════════════════
# Step 1:  Validation
# ═══════════════════════════════════════════════════════════════════

def validate_seed_ranges(combos: List[Dict[str, Any]]) -> bool:
    """Verify that all seed intervals across combos are non-overlapping.

    Prints a diagnostic table and returns True if no conflicts are found.
    """
    intervals: List[Tuple[int, int, str]] = []

    for combo in combos:
        for ti in range(N_WORLD):
            topo_seed = compute_topology_seed(combo, ti)
            world_seed = compute_world_seed(topo_seed, combo['strength_idx'])
            cell_min = compute_cell_seed(world_seed, 0)
            cell_max = compute_cell_seed(world_seed, N_INIT - 1)
            intervals.append((cell_min, cell_max, combo['combo_key']))

    intervals.sort()
    ok = True
    print("\n── Seed Range Validation ──")

    for lo, hi, ck in intervals:
        print(f"  [{lo:>6}, {hi:>6}]  {ck}")

    for i in range(len(intervals) - 1):
        _, hi1, ck1 = intervals[i]
        lo2, _, ck2 = intervals[i + 1]
        if hi1 >= lo2:
            print(f"\n  *** OVERLAP *** {ck1} [{hi1}] vs {ck2} [{lo2}]")
            ok = False

    if ok:
        print(f"\n  All {len(intervals)} cell-seed intervals are non-overlapping. ✓")
    return ok


def print_grid_summary(combos: List[Dict[str, Any]]) -> None:
    """Print a summary of the experimental grid."""
    print(f"\n── Full Grid Summary ──")
    print(f"  topology_types:    {TOPOLOGY_TYPES}")
    print(f"  δ values:          {DELTA_VALUES}")
    print(f"  strength regimes:  {list(STRENGTH_REGIMES.keys())}")
    print(f"  sign ratios:       {list(SIGN_RATIOS.keys())}")
    print(f"  Total combos:      {len(combos)}")
    print(f"  Worlds per combo:  {N_WORLD}")
    print(f"  Initial states:    {N_INIT}")
    print(f"  Total worlds:      {len(combos) * N_WORLD}")
    print(f"  Total trajectories: {len(combos) * N_WORLD * N_INIT}")
    print(f"  T per trajectory:  {T_SIM}")
    print(f"  execution_grid:    {EXECUTION_GRID}")

    # per-topology breakdown
    for topo in TOPOLOGY_TYPES:
        n = sum(1 for c in combos if c['topology_type'] == topo)
        print(f"    {topo}: {n} combos → {n * N_WORLD} worlds")


# ═══════════════════════════════════════════════════════════════════
# Step 2:  World Generation & Trajectory Simulation  (stubs)
# ═══════════════════════════════════════════════════════════════════

def get_model_module(topology_type: str):
    """Return the correct model module based on topology type.

    Args:
        topology_type: 'unrestricted_sparse' → ddc_model_a
                       'tf_restricted'       → ddc_model_a_tf
    """
    if topology_type == 'unrestricted_sparse':
        import ddc_model_a as dma
        return dma
    elif topology_type == 'tf_restricted':
        import ddc_model_a_tf as dma
        return dma
    else:
        raise ValueError(f"Unknown topology_type: {topology_type}")


def build_world_metadata(world, world_id: str, combo: Dict[str, Any],
                         topo_idx: int, topo_seed: int,
                         world_seed: int) -> Dict[str, Any]:
    """Construct the metadata dict for a single world.

    Per Output Spec §4, must include: world_id, execution_grid,
    topology_type, δ, self_retention, b, strength_regime, sign_ratio,
    combo_key, replicate_index, edge_count, edge sign map,
    edge strength distribution, seeds, T, N_init, runtime_version,
    and topology metadata (candidate count, edge densities,
    in/out-degree distributions, etc.).
    """
    dma = get_model_module(combo['topology_type'])
    wdict = world.to_dict()

    wdict.update({
        # ── core identity ──
        'world_id':         world_id,
        'execution_grid':   EXECUTION_GRID,
        'runtime_version':  RUNTIME_VERSION,

        # ── parameter axes ──
        'topology_type':    combo['topology_type'],
        'delta':            combo['delta'],
        'self_retention':   combo['self_retention'],
        'b':                DEFAULT_B,
        'strength_regime':  combo['strength_regime'],
        'sign_ratio':       combo['sign_ratio'],
        'sign_ratio_value': combo['sign_ratio_value'],

        # ── identification ──
        'combo_key':        combo['combo_key'],
        'replicate_index':  topo_idx,
        'topo_idx':         topo_idx,
        'delta_label':      combo['delta_label'],

        # ── topology summary ──
        'edge_count':       dma.N_EDGES,
        'edge_density_r':   dma.EDGE_DENSITY_R,
        'G':                dma.G,
        # candidate_edge_count per Output Spec §4
        'candidate_edge_count': (dma.G * dma.G - dma.G
                                 if combo['topology_type'] == 'unrestricted_sparse'
                                 else dma.N_TF * dma.G - dma.N_TF),

        # ── runtime config ──
        't_sim':            T_SIM,
        'n_init':           N_INIT,
        'topo_seed':        topo_seed,
        'world_seed':       world_seed,
    })

    # ── world.to_dict() already wrote 'seed' = world_seed; remove duplicate ──
    wdict.pop('seed', None)

    # ── strength regime parameters ──
    wdict['a_min'] = combo['a_min']
    wdict['a_max'] = combo['a_max']

    # ── hub-likeness (per Output Spec §4, matched to run_002) ──
    # tf_restricted is inherently hub-like (edges restricted to TF pool)
    # unrestricted_sparse: out_degree_quartile assigned post-hoc via assign_hub_likeness()
    if combo['topology_type'] == 'tf_restricted':
        wdict['hub_likeness'] = True
        wdict['out_degree_quartile'] = None
    else:
        out_deg_vals = np.array(list(world.out_degrees().values()), dtype=float)
        wdict['_out_degree_std'] = float(np.std(out_deg_vals, ddof=0))  # temporary
        wdict['hub_likeness'] = False
        wdict['out_degree_quartile'] = None

    # ── TF-specific metadata (Output Spec §4) ──
    if combo['topology_type'] == 'tf_restricted':
        import ddc_model_a_tf
        tf_meta = ddc_model_a_tf.tf_metadata(world)
        wdict.update(tf_meta)
        wdict['tf_genes'] = ddc_model_a_tf.TF_GENES

    return wdict


def simulate_one_combo(combo: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Generate worlds and simulate trajectories for a single combo.

    For world i (0..9):
      1. sample topology (seed = sign_base + δ_offset + i*100)
      2. sample world with edge strengths (seed = topo_seed + 50 + strength_idx*10)
      3. override world.delta ← combo['delta']
      4. for cell j (0..4):
           sample initial state, run simulation, save trajectory

    Returns: list of per-world result dicts (will be extended by step 3 analysis).
    """
    dma = get_model_module(combo['topology_type'])
    results: List[Dict[str, Any]] = []

    print(f"\n  combo: {combo['combo_key']}  "
          f"(topo={combo['topology_type']}, δ={combo['delta']}, "
          f"strength=[{combo['a_min']},{combo['a_max']}], "
          f"sign={combo['sign_ratio']})")

    for ti in range(N_WORLD):
        # ── topology ──
        topo_seed = compute_topology_seed(combo, ti)
        w_topo = dma.sample_topology(topo_seed, sign_ratio=combo['sign_ratio_value'])

        # ── world ──
        world_seed = compute_world_seed(topo_seed, combo['strength_idx'])
        world = dma.sample_world(world_seed, world=w_topo,
                                 a_min=combo['a_min'], a_max=combo['a_max'])

        # ── override δ ──
        world.delta = torch.full((dma.G,), combo['delta'], dtype=dma.DTYPE)
        world._A = None  # force rebuild of transition matrix

        # ── metadata ──
        world_id = make_world_id(
            combo['topology_type'], combo['delta_label'],
            combo['strength_regime'], combo['sign_ratio'], ti,
        )
        wdict = build_world_metadata(world, world_id, combo, ti, topo_seed, world_seed)
        os.makedirs(METADATA_DIR, exist_ok=True)
        with open(os.path.join(METADATA_DIR, f'{world_id}.json'), 'w') as f:
            json.dump(wdict, f, indent=2)

        # ── trajectories ──
        per_world_analyses = []
        for ci in range(N_INIT):
            cell_seed = compute_cell_seed(world_seed, ci)
            X0 = dma.sample_initial_state(cell_seed)
            traj = dma.simulate_single_cell(world, X0, t_steps=T_SIM)

            # TODO: placeholder — analysis will be populated in step 3
            analysis_placeholder = {
                'world_id':   world_id,
                'cell_seed':  cell_seed,
                'cell_idx':   ci,
            }
            per_world_analyses.append(analysis_placeholder)

            # ── save trajectory ──
            traj_save = {
                'X_traj':        traj['X_traj'],
                'clip_count':    traj['clip_count'],
                'total_clips':   traj['total_clips'],
                'world_id':      world_id,
                'world_seed':    world_seed,
                'cell_seed':     cell_seed,
                'delta':         combo['delta'],
                'topology_type': combo['topology_type'],
            }
            os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
            torch.save(traj_save,
                       os.path.join(TRAJECTORIES_DIR, f'{world_id}_cell{ci:02d}.pt'))

        # ── per-world result stub ──
        result_stub = {
            'world_id':       world_id,
            'world_seed':     world_seed,
            'combo_key':      combo['combo_key'],
            'topo_idx':       ti,
            'topology_type':  combo['topology_type'],
            'delta':          combo['delta'],
            'strength_regime': combo['strength_regime'],
            'a_min':          combo['a_min'],
            'a_max':          combo['a_max'],
            'sign_ratio':     combo['sign_ratio'],
            'sign_ratio_value': combo['sign_ratio_value'],
            'per_cell':       per_world_analyses,
        }
        results.append(result_stub)

    return results


def simulate_all_combos(combos: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Run world generation + trajectory simulation for all combos.

    Returns: combined list of per-world result dicts for all combos.
    """
    all_results: List[Dict[str, Any]] = []

    for i, combo in enumerate(combos):
        print(f"\n=== [{i+1}/{len(combos)}] {combo['combo_key']} ===")
        results = simulate_one_combo(combo)
        all_results.extend(results)

    return all_results


def assign_hub_likeness() -> None:
    """Post-process unrestricted_sparse world metadata: compute out-degree std
    quartiles across all unrestricted worlds, then assign out_degree_quartile
    and hub_likeness flag.

    tf_restricted worlds are already marked hub_likeness=True (inherently hub-like).

    Q3 and Q4 → hub_likeness = True  (topologies with concentrated regulatory influence)
    Q1 and Q2 → hub_likeness = False

    After assignment, the temporary _out_degree_std field is removed.
    Per run_002 convention: quartile bins based on cross-world ranking.
    """
    # 1. Collect _out_degree_std only from unrestricted_sparse metadata files
    world_stds: List[Tuple[str, float]] = []
    for fname in os.listdir(METADATA_DIR):
        if not fname.endswith('.json'):
            continue
        fpath = os.path.join(METADATA_DIR, fname)
        with open(fpath) as f:
            meta = json.load(f)
        std_val = meta.get('_out_degree_std', None)
        if std_val is not None:
            world_stds.append((fpath, std_val))

    if not world_stds:
        print("  Warning: no unrestricted_sparse world metadata found for hub-likeness assignment")
        return

    # 2. Compute quartile boundaries
    stds = np.array([ws[1] for ws in world_stds])
    q25, q50, q75 = np.percentile(stds, [25, 50, 75])

    def _quartile(v: float) -> str:
        if v < q25:
            return 'Q1 (low)'
        if v < q50:
            return 'Q2'
        if v < q75:
            return 'Q3'
        return 'Q4 (high)'

    # 3. Assign out_degree_quartile and hub_likeness, write back
    for fpath, std_val in world_stds:
        with open(fpath) as f:
            meta = json.load(f)
        q = _quartile(std_val)
        meta['out_degree_quartile'] = q
        meta['hub_likeness'] = q.startswith('Q3') or q.startswith('Q4')
        del meta['_out_degree_std']
        with open(fpath, 'w') as f:
            json.dump(meta, f, indent=2)

    print(f"\n  Hub-likeness assigned for {len(world_stds)} unrestricted_sparse worlds "
          f"(Q1≤{q25:.3f}, Q2≤{q50:.3f}, Q3≤{q75:.3f})")


# ═══════════════════════════════════════════════════════════════════
# Step 3:  Per-Trajectory Analysis  (stubs)
# ═══════════════════════════════════════════════════════════════════

def detect_equilibrium(X_traj: torch.Tensor) -> Dict[str, Any]:
    """Detect convergence & compute equilibrium statistics.

    Per Output Spec §8: window = 50, ε = 1e-4.
    Scans forward from t=0 to find the first step where the
    following 50 consecutive steps all have ‖ΔX‖ < ε.
    """
    T = X_traj.shape[0] - 1  # t_steps
    converged = False
    convergence_time = None

    # Scan forward: find the earliest step from which the next
    # 50 steps are all within ε.
    if T >= CONVERGENCE_WINDOW:
        for t in range(T - CONVERGENCE_WINDOW + 1):
            all_ok = True
            for w in range(CONVERGENCE_WINDOW):
                diff = torch.norm(X_traj[t + w + 1] - X_traj[t + w]).item()
                if diff >= EPSILON:
                    all_ok = False
                    break
            if all_ok:
                converged = True
                convergence_time = t
                break

    # Equilibrium statistics (final window, even if not formally converged)
    x_final_window = X_traj[-CONVERGENCE_WINDOW:, :]  # (50, G)
    x_eq = x_final_window.mean(dim=0)  # (G,)
    equilibrium_magnitude = float(x_eq.mean().item())
    equilibrium_sparsity = float((x_eq < SPARSITY_THRESHOLD).float().mean().item())
    final_window_variance = float(x_final_window.var(dim=0).mean().item())

    return {
        'converged': converged,
        'convergence_time': convergence_time,
        'equilibrium_magnitude': equilibrium_magnitude,
        'equilibrium_sparsity': equilibrium_sparsity,
        'final_window_variance': final_window_variance,
    }


def analyze_stability(X_traj: torch.Tensor, clip_count: torch.Tensor,
                      total_clips: int) -> Dict[str, Any]:
    """Analyze boundedness, divergence, collapse, and clipping.

    Per Output Spec §5 & §7.
    """
    T = X_traj.shape[0] - 1
    G = X_traj.shape[1]

    # Global abundance-scale metrics (Output Spec §5)
    max_amplitude = float(X_traj.max().item())
    final_amplitude = float(X_traj[-1].max().item())

    # Boundedness
    is_bounded = max_amplitude < DIVERGENCE_THRESHOLD

    # Divergence: first step where any gene exceeds DIVERGENCE_THRESHOLD
    divergence_time = None
    for t in range(X_traj.shape[0]):
        if X_traj[t].max().item() >= DIVERGENCE_THRESHOLD:
            divergence_time = t
            break

    # Collapse: final mean expression < threshold
    final_mean = float(X_traj[-1].mean().item())
    is_collapse = final_mean < COLLAPSE_THRESHOLD

    # Clipping
    # total_clips counts events across T update steps (t=1..T)
    # per Analysis Plan §7: clipping_frequency = clipping_events / (T × G)
    clipping_frequency = total_clips / (T * G) if T > 0 and G > 0 else 0.0
    clipping_dominated = clipping_frequency > CLIPPING_FRAC_THRESHOLD

    return {
        'max_amplitude': max_amplitude,
        'final_amplitude': final_amplitude,
        'is_bounded': is_bounded,
        'is_collapse': is_collapse,
        'divergence_time': divergence_time,
        'clipping_frequency': clipping_frequency,
        'clipping_dominated': clipping_dominated,
    }


def analyze_oscillation(X_traj: torch.Tensor,
                        converged: bool, conv_time: int) -> Dict[str, Any]:
    """Detect oscillation, classify sustained vs damped.

    Per Output Spec §9; detection logic matched to run_002.
    Uses local extrema detection (skip BURN_IN=200), per-gene
    peak-pair amplitude, and early/late damping ratio.
    """
    BURN_IN = 200
    T_total = X_traj.shape[0]

    # Slice: skip burn-in; if converged, stop at conv_time
    if converged and conv_time is not None and conv_time > BURN_IN + 100:
        X = X_traj[BURN_IN:conv_time].numpy()
    else:
        X = X_traj[BURN_IN:].numpy()

    T, G_dim = X.shape
    if T < 50:
        return {
            'oscillation_exists': False,
            'is_sustained': False,
            'is_damped_transient': False,
            'oscillation_amplitude': None,
            'oscillatory_genes': None,
            'estimated_period': None,
        }

    oscillatory_genes = []
    amplitudes = []
    frequencies = []
    damping_rates = []

    MIN_REL_AMP = 0.01

    for g in range(G_dim):
        x = X[:, g]
        signal_range = float(x.max() - x.min())
        if signal_range < EPSILON:
            continue

        # Find local extrema
        extrema = np.zeros(T, dtype=np.int8)
        for t in range(1, T - 1):
            if x[t] > x[t - 1] and x[t] > x[t + 1]:
                extrema[t] = 1   # peak
            elif x[t] < x[t - 1] and x[t] < x[t + 1]:
                extrema[t] = -1  # trough

        ext_idx = np.where(extrema != 0)[0]
        if len(ext_idx) < 3:
            continue

        # Consecutive peak-trough / trough-peak amplitudes
        peak_pairs = []
        for i in range(len(ext_idx) - 1):
            s, e = ext_idx[i], ext_idx[i + 1]
            peak_pairs.append(abs(float(x[e] - x[s])))

        if len(peak_pairs) < 2:
            continue

        gene_amplitude = float(np.median(peak_pairs))
        relative_amplitude = gene_amplitude / signal_range if signal_range > EPSILON else 0.0

        if relative_amplitude < MIN_REL_AMP:
            continue

        gene_freq = len(ext_idx) / (2.0 * T)

        # Damping: compare early vs late half of peak pairs
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
            'is_sustained': False,
            'is_damped_transient': False,
            'oscillation_amplitude': None,
            'oscillatory_genes': None,
            'estimated_period': None,
        }

    avg_damping = float(np.median(damping_rates))
    is_sustained = avg_damping <= 0.05
    is_damped_transient = not is_sustained

    mean_amplitude = float(np.mean(amplitudes))
    mean_frequency = float(np.mean(frequencies))
    estimated_period = 1.0 / mean_frequency if mean_frequency > 0 else None

    # oscillation_amplitude per Output Spec §9:
    # peak-to-trough amplitude in the final analysis window
    if len(oscillatory_genes) > 0:
        final_win = X_traj[-CONVERGENCE_WINDOW:, oscillatory_genes]
        oscillation_amplitude = float((final_win.max() - final_win.min()).item())
    else:
        oscillation_amplitude = None

    return {
        'oscillation_exists': True,
        'is_sustained': is_sustained,
        'is_damped_transient': is_damped_transient,
        'oscillation_amplitude': oscillation_amplitude,
        'oscillatory_genes': oscillatory_genes,
        'estimated_period': estimated_period,
    }


def classify_primary_regime(eq: Dict, st: Dict, osc: Dict) -> Tuple[str, List[str], str]:
    """Return (primary_regime, secondary_labels, legacy_type).

    Priority order (per Problem Definition §7, Analysis Plan §8):
      Collapse > Divergent > Sustained oscillatory > Convergent > Ambiguous

    Secondary labels as per Output Spec §6.
    """
    secondary: List[str] = []

    # Add clipping secondary labels
    if st['clipping_dominated']:
        secondary.append('clipping_dominated')
    if st['clipping_frequency'] < 0.01 and eq['converged']:
        secondary.append('low_clipping_linear_like')

    # Classification by priority
    if st['is_collapse']:
        primary = 'Collapse'
        secondary.append('numerical_collapse')
        legacy = 'Type F'

    elif not st['is_bounded']:
        primary = 'Divergent'
        secondary.append('runaway_divergence')
        legacy = 'Type E'

    elif osc['is_sustained'] and not eq['converged']:
        primary = 'Sustained oscillatory'
        legacy = 'Type D'

    elif eq['converged']:
        primary = 'Convergent'
        conv_t = eq['convergence_time']
        if conv_t is not None and conv_t <= SLOW_CONVERGENCE_THRESHOLD:
            secondary.append('fast_convergence')
            legacy = 'Type A'
        else:
            secondary.append('slow_convergence')
            legacy = 'Type B'
        if osc['is_damped_transient']:
            secondary.append('damped_oscillatory_transient')
            legacy = 'Type C'

    else:
        primary = 'Ambiguous'
        legacy = 'Type G'

    return primary, secondary, legacy


def compute_spectral_info(world_dict: Dict[str, Any]) -> Tuple[float, Dict]:
    """Compute spectral radius ρ(A) and eigenvalues of the linear-core A.

    Per Output Spec §7: A is the clipping-free transition matrix.
    A = diag(1 - δ_i) + signed_strength_matrix
    """
    G_val = world_dict.get('G', 20)
    delta_val = world_dict.get('delta', 0.2)

    # Build A: start with diagonal (1 - δ)
    A = np.diag(np.full(G_val, 1.0 - delta_val))

    # Add signed edge strengths
    edge_signs = world_dict.get('edge_signs', {})
    edge_strengths = world_dict.get('edge_strengths', {})

    for i_str, targets in edge_signs.items():
        i = int(i_str)
        for j_str, sign in targets.items():
            j = int(j_str)
            strength_key = str(j_str)
            i_edge = edge_strengths.get(i_str, {})
            strength = i_edge.get(strength_key, 0.0)
            A[i, j] += float(sign) * float(strength)

    # Eigen decomposition
    eigvals = np.linalg.eigvals(A)
    spectral_radius = float(np.max(np.abs(eigvals)))

    return spectral_radius, {
        'eigenvalues_real': [float(v.real) for v in eigvals],
        'eigenvalues_imag': [float(v.imag) for v in eigvals],
        'spectral_radius': spectral_radius,
    }


def analyze_trajectory(X_traj: torch.Tensor, clip_count: torch.Tensor,
                       total_clips: int, world_dict: Dict[str, Any]) -> Dict[str, Any]:
    """Full per-trajectory analysis: eq + stability + oscillation + classify + ρ(A).

    Returns a flat dict with all analysis fields merged.
    """
    eq  = detect_equilibrium(X_traj)
    st  = analyze_stability(X_traj, clip_count, total_clips)
    osc = analyze_oscillation(X_traj, eq['converged'], eq['convergence_time'])
    primary_regime, secondary_labels, legacy_type = classify_primary_regime(eq, st, osc)
    spectral_radius, spectral_dict = compute_spectral_info(world_dict)

    return {
        # ── equilibrium ──
        'converged': eq['converged'],
        'convergence_time': eq['convergence_time'],
        'equilibrium_magnitude': eq['equilibrium_magnitude'],
        'equilibrium_sparsity': eq['equilibrium_sparsity'],
        'final_window_variance': eq['final_window_variance'],
        # ── stability ──
        'max_amplitude': st['max_amplitude'],
        'final_amplitude': st['final_amplitude'],
        'is_bounded': st['is_bounded'],
        'is_collapse': st['is_collapse'],
        'divergence_time': st['divergence_time'],
        'clipping_frequency': st['clipping_frequency'],
        'clipping_dominated': st['clipping_dominated'],
        # ── oscillation ──
        'oscillation_exists': osc['oscillation_exists'],
        'is_sustained': osc['is_sustained'],
        'is_damped_transient': osc['is_damped_transient'],
        'oscillation_amplitude': osc['oscillation_amplitude'],
        'oscillatory_genes': osc['oscillatory_genes'],
        'estimated_period': osc['estimated_period'],
        # ── classification ──
        'primary_regime': primary_regime,
        'secondary_labels': secondary_labels,
        'legacy_type': legacy_type,
        # ── spectral ──
        'spectral_radius': spectral_radius,
        'eigenvalues_real': spectral_dict['eigenvalues_real'],
        'eigenvalues_imag': spectral_dict['eigenvalues_imag'],
    }


def run_per_trajectory_analysis(all_results: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Load all trajectories from disk and run per-trajectory analysis.

    Updates all_results in-place: fills each per_cell entry with full analysis fields.
    Also saves per-trajectory analysis JSONs to analysis/.
    """
    print("\n── Step 3: Per-Trajectory Analysis ──")

    total_trajs = 0

    for world_result in all_results:
        world_id = world_result['world_id']

        # Load world metadata
        meta_path = os.path.join(METADATA_DIR, f'{world_id}.json')
        with open(meta_path) as f:
            world_dict = json.load(f)

        for cell_entry in world_result['per_cell']:
            ci = cell_entry['cell_idx']

            # Load trajectory from disk
            traj_path = os.path.join(TRAJECTORIES_DIR, f'{world_id}_cell{ci:02d}.pt')
            traj_data = torch.load(traj_path, map_location='cpu', weights_only=False)

            X_traj = traj_data['X_traj']
            clip_count = traj_data['clip_count']
            total_clips = traj_data['total_clips']

            # Run analysis
            analysis = analyze_trajectory(X_traj, clip_count, total_clips, world_dict)

            # Fill per_cell entry (replace placeholder)
            cell_entry.update(analysis)

            # Save per-trajectory analysis
            os.makedirs(ANALYSIS_DIR, exist_ok=True)
            analysis_path = os.path.join(ANALYSIS_DIR, f'{world_id}_cell{ci:02d}.json')
            with open(analysis_path, 'w') as f:
                json.dump({
                    'world_id': world_id,
                    'cell_idx': ci,
                    'cell_seed': cell_entry['cell_seed'],
                    **analysis,
                }, f, indent=2)

            total_trajs += 1

    print(f"  {total_trajs} trajectories analyzed.")
    return all_results


# ═══════════════════════════════════════════════════════════════════
# Step 4:  Aggregation
# ═══════════════════════════════════════════════════════════════════

def aggregate_world(per_cell_results: List[Dict]) -> Dict[str, Any]:
    """World-level majority vote on primary regime, with per-cell means.

    Input: list of 5 per_cell dicts (already filled with analysis fields).
    Returns: dict with world-level primary_regime and mean metrics.
    """
    # Majority vote on primary regime
    regime_counts = Counter(c['primary_regime'] for c in per_cell_results)
    primary_regime = regime_counts.most_common(1)[0][0]

    # Majority vote on is_sustained
    sustained_votes = sum(1 for c in per_cell_results if c['is_sustained'])
    sustained_oscillation = sustained_votes >= 3

    # Majority vote on divergence
    divergent_votes = sum(1 for c in per_cell_results if c['primary_regime'] == 'Divergent')
    divergence = divergent_votes >= 3

    # Majority vote on collapse
    collapse_votes = sum(1 for c in per_cell_results if c['primary_regime'] == 'Collapse')
    collapse = collapse_votes >= 3

    # Collect secondary labels across cells
    all_secondary = []
    for c in per_cell_results:
        all_secondary.extend(c['secondary_labels'])
    secondary_labels = sorted(set(all_secondary))

    # Mean numerical fields
    all_converged = all(c['converged'] for c in per_cell_results)
    if all_converged:
        mean_convergence_time = float(np.mean(
            [c['convergence_time'] for c in per_cell_results]))
    else:
        mean_convergence_time = None

    mean_clipping_frequency = float(np.mean([c['clipping_frequency'] for c in per_cell_results]))
    # Matched to run_002/003: world is clipping_dominated if any cell is
    clipping_dominated = any(c['clipping_dominated'] for c in per_cell_results)

    # spectral_radius is a world-level property (same A matrix for all cells)
    # Matched to run_002/003: take from first cell
    spectral_radius = per_cell_results[0]['spectral_radius']
    mean_equilibrium_magnitude = float(np.mean([c['equilibrium_magnitude'] for c in per_cell_results]))
    mean_equilibrium_sparsity = float(np.mean([c['equilibrium_sparsity'] for c in per_cell_results]))

    max_amplitude = float(max(c['max_amplitude'] for c in per_cell_results))
    final_amplitude = float(max(c['final_amplitude'] for c in per_cell_results))

    return {
        'primary_regime': primary_regime,
        'secondary_labels': secondary_labels,
        'convergence_time': mean_convergence_time,
        'sustained_oscillation': sustained_oscillation,
        'divergence': divergence,
        'collapse': collapse,
        'clipping_frequency': mean_clipping_frequency,
        'clipping_dominated': clipping_dominated,
        'max_amplitude': max_amplitude,
        'final_amplitude': final_amplitude,
        'spectral_radius': spectral_radius,
        'equilibrium_magnitude': mean_equilibrium_magnitude,
        'equilibrium_sparsity': mean_equilibrium_sparsity,
        'regime_tally': dict(regime_counts),
    }


def aggregate_across_combos(all_results: List[Dict]) -> Dict[str, Any]:
    """Combo-level statistics.

    Per Output Spec §11.2:
    For each combo_key, aggregate across its N_WORLD worlds.
    Returns a dict keyed by combo_key.
    """
    print("\n── Step 4: Aggregation ──")

    # Step A: Aggregate each world
    world_aggregates = {}
    for world_result in all_results:
        wid = world_result['world_id']
        wa = aggregate_world(world_result['per_cell'])
        world_aggregates[wid] = {
            **wa,
            'combo_key': world_result['combo_key'],
            'topology_type': world_result['topology_type'],
            'delta': world_result['delta'],
            'strength_regime': world_result['strength_regime'],
            'sign_ratio': world_result['sign_ratio'],
        }

    # Step B: Group by combo_key
    combo_groups: Dict[str, List[Dict]] = defaultdict(list)
    for wid, wa in world_aggregates.items():
        combo_groups[wa['combo_key']].append(wa)

    combo_stats: Dict[str, Dict] = {}
    for ck, worlds in combo_groups.items():
        n = len(worlds)
        regimes = [w['primary_regime'] for w in worlds]
        reg_counter = Counter(regimes)

        conv_times = [w['convergence_time'] for w in worlds
                      if w['convergence_time'] is not None]
        cfs = [w['clipping_frequency'] for w in worlds]
        srs = [w['spectral_radius'] for w in worlds]

        # Representative metadata from first world
        rep = worlds[0]

        combo_stats[ck] = {
            'combo_key': ck,
            'topology_type': rep['topology_type'],
            'delta': rep['delta'],
            'strength_regime': rep['strength_regime'],
            'sign_ratio': rep['sign_ratio'],
            'N_world': n,
            'fraction_convergent': reg_counter.get('Convergent', 0) / n,
            'fraction_sustained_oscillatory': reg_counter.get('Sustained oscillatory', 0) / n,
            'fraction_divergent': reg_counter.get('Divergent', 0) / n,
            'fraction_collapse': reg_counter.get('Collapse', 0) / n,
            'fraction_ambiguous': reg_counter.get('Ambiguous', 0) / n,
            'mean_convergence_time': float(np.mean(conv_times)) if conv_times else None,
            'mean_clipping_frequency': float(np.mean(cfs)),
            'mean_spectral_radius': float(np.mean(srs)),
            'regime_tally': dict(reg_counter),
        }

    print(f"  {len(combo_stats)} combos aggregated.")
    return combo_stats


def analyze_cross_dimension(all_results: List[Dict]) -> Dict[str, Any]:
    """Cross-dimensional interaction statistics.

    Per Analysis Plan §9.2-9.4:
      topology_type × δ
      δ × strength_regime
      sign_ratio × δ
    """
    # Build per-world aggregates first
    world_agg = {}
    for world_result in all_results:
        wid = world_result['world_id']
        wa = aggregate_world(world_result['per_cell'])
        world_agg[wid] = {
            **wa,
            'topology_type': world_result['topology_type'],
            'delta': world_result['delta'],
            'strength_regime': world_result['strength_regime'],
            'sign_ratio': world_result['sign_ratio'],
        }

    worlds = list(world_agg.values())

    def _cross_stats(group_key_fn, dimensions) -> Dict:
        """Helper: group worlds by key_fn, compute per-group stats."""
        groups: Dict[Tuple, List[Dict]] = defaultdict(list)
        for w in worlds:
            key = group_key_fn(w)
            groups[key].append(w)

        result = {}
        for key, gw in sorted(groups.items()):
            n = len(gw)
            regimes = [w['primary_regime'] for w in gw]
            reg_counter = Counter(regimes)
            cfs = [w['clipping_frequency'] for w in gw]
            srs = [w['spectral_radius'] for w in gw]
            conv_times = [w['convergence_time'] for w in gw
                          if w['convergence_time'] is not None]

            label = '_x_'.join(str(k) for k in key) if isinstance(key, tuple) else str(key)
            result[label] = {
                'dimensions': {dim: k for dim, k in zip(dimensions, key)} if isinstance(key, tuple) else {dimensions[0]: key},
                'N_world': n,
                'fraction_convergent': reg_counter.get('Convergent', 0) / n,
                'fraction_sustained_oscillatory': reg_counter.get('Sustained oscillatory', 0) / n,
                'fraction_divergent': reg_counter.get('Divergent', 0) / n,
                'fraction_collapse': reg_counter.get('Collapse', 0) / n,
                'fraction_ambiguous': reg_counter.get('Ambiguous', 0) / n,
                'mean_clipping_frequency': float(np.mean(cfs)),
                'mean_spectral_radius': float(np.mean(srs)),
                'mean_convergence_time': float(np.mean(conv_times)) if conv_times else None,
                'regime_tally': dict(reg_counter),
            }
        return result

    cross_dim = {
        'topology_x_delta': _cross_stats(
            lambda w: (w['topology_type'], w['delta']),
            ['topology_type', 'delta']
        ),
        'delta_x_strength_regime': _cross_stats(
            lambda w: (w['delta'], w['strength_regime']),
            ['delta', 'strength_regime']
        ),
        'sign_ratio_x_delta': _cross_stats(
            lambda w: (w['sign_ratio'], w['delta']),
            ['sign_ratio', 'delta']
        ),
    }

    print("  Cross-dimension analysis complete.")
    return cross_dim


# ═══════════════════════════════════════════════════════════════════
# Step 5:  Tables
# ═══════════════════════════════════════════════════════════════════

def build_world_summary_table(all_results: List[Dict]) -> List[Dict]:
    """Per Output Spec §11.1: one row per world with 23 fields."""
    rows: List[Dict] = []

    for world_result in all_results:
        world_id = world_result['world_id']

        # Load metadata
        meta_path = os.path.join(METADATA_DIR, f'{world_id}.json')
        with open(meta_path) as f:
            meta = json.load(f)

        # World-level aggregate
        wa = aggregate_world(world_result['per_cell'])

        rows.append({
            'world_id':             world_id,
            'execution_grid':       meta.get('execution_grid', ''),
            'topology_type':        world_result['topology_type'],
            'delta':                world_result['delta'],
            'self_retention':       meta.get('self_retention', 1.0 - world_result['delta']),
            'b':                    meta.get('b', DEFAULT_B),
            'strength_regime':      world_result['strength_regime'],
            'sign_ratio':           world_result['sign_ratio'],
            'replicate_index':      world_result['topo_idx'],
            'edge_count':           meta.get('edge_count', N_EDGES),
            'primary_regime':       wa['primary_regime'],
            'secondary_labels':     ','.join(wa['secondary_labels']),
            'convergence_time':     wa['convergence_time'],
            'sustained_oscillation': wa['sustained_oscillation'],
            'divergence':           wa['divergence'],
            'collapse':             wa['collapse'],
            'clipping_frequency':   wa['clipping_frequency'],
            'clipping_dominated':   wa['clipping_dominated'],
            'max_amplitude':        wa['max_amplitude'],
            'final_amplitude':      wa['final_amplitude'],
            'spectral_radius':      wa['spectral_radius'],
            'equilibrium_magnitude': wa['equilibrium_magnitude'],
            'equilibrium_sparsity': wa['equilibrium_sparsity'],
        })

    return rows


def build_regime_summary_table(aggregation: Dict, all_results: List[Dict]) -> List[Dict]:
    """Per Output Spec §11.2: combo-level regime summary."""
    # aggregation is combo_stats dict keyed by combo_key
    rows: List[Dict] = []
    for ck in sorted(aggregation.keys()):
        cs = aggregation[ck]
        rows.append({
            'combo_key':                    cs['combo_key'],
            'topology_type':                cs['topology_type'],
            'delta':                        cs['delta'],
            'strength_regime':              cs['strength_regime'],
            'sign_ratio':                   cs['sign_ratio'],
            'N_world':                      cs['N_world'],
            'fraction_convergent':          cs['fraction_convergent'],
            'fraction_sustained_oscillatory': cs['fraction_sustained_oscillatory'],
            'fraction_divergent':           cs['fraction_divergent'],
            'fraction_collapse':            cs['fraction_collapse'],
            'fraction_ambiguous':           cs['fraction_ambiguous'],
            'mean_convergence_time':        cs['mean_convergence_time'],
            'mean_clipping_frequency':      cs['mean_clipping_frequency'],
            'mean_spectral_radius':         cs['mean_spectral_radius'],
        })
    return rows


def build_boundary_summary_table(all_results: List[Dict]) -> List[Dict]:
    """Per Output Spec §11.3: boundary summary.

    Aggregated by topology_type × δ, and by topology_type × δ × strength_regime.
    """
    # Assemble world-level data
    world_agg = {}
    for world_result in all_results:
        wid = world_result['world_id']
        wa = aggregate_world(world_result['per_cell'])
        world_agg[wid] = {
            **wa,
            'topology_type': world_result['topology_type'],
            'delta': world_result['delta'],
            'strength_regime': world_result['strength_regime'],
        }

    worlds = list(world_agg.values())

    rows: List[Dict] = []

    # topology_type × δ level
    for topo in TOPOLOGY_TYPES:
        topo_worlds = [w for w in worlds if w['topology_type'] == topo]
        for d in DELTA_VALUES:
            td_worlds = [w for w in topo_worlds if w['delta'] == d]
            n = len(td_worlds)
            if n == 0:
                continue
            regimes = [w['primary_regime'] for w in td_worlds]
            reg_counter = Counter(regimes)
            clipping_frac = sum(1 for w in td_worlds if w['clipping_dominated']) / n
            sustained_count = sum(1 for w in td_worlds if w['sustained_oscillation'])

            # First strength regime with substantial divergence (> 30%)
            first_div = None
            for sr in STRENGTH_ORDER:
                sr_worlds = [w for w in td_worlds if w['strength_regime'] == sr]
                if sr_worlds:
                    div_frac = sum(1 for w in sr_worlds if w['primary_regime'] == 'Divergent') / len(sr_worlds)
                    if div_frac > 0.3:
                        first_div = sr
                        break

            rows.append({
                'level': 'topology_x_delta',
                'topology_type': topo,
                'delta': d,
                'strength_regime': '',
                'N_world': n,
                'fraction_convergent': reg_counter.get('Convergent', 0) / n,
                'fraction_divergent': reg_counter.get('Divergent', 0) / n,
                'fraction_collapse': reg_counter.get('Collapse', 0) / n,
                'fraction_sustained_oscillatory': reg_counter.get('Sustained oscillatory', 0) / n,
                'clipping_dominated_fraction': clipping_frac,
                'sustained_oscillation_count': sustained_count,
                'first_divergent_regime': first_div,
            })

    # topology_type × δ × strength_regime level
    for topo in TOPOLOGY_TYPES:
        for d in DELTA_VALUES:
            for sr in STRENGTH_ORDER:
                td_worlds = [w for w in worlds
                             if w['topology_type'] == topo
                             and w['delta'] == d
                             and w['strength_regime'] == sr]
                n = len(td_worlds)
                if n == 0:
                    continue
                regimes = [w['primary_regime'] for w in td_worlds]
                reg_counter = Counter(regimes)
                clipping_frac = sum(1 for w in td_worlds if w['clipping_dominated']) / n
                sustained_count = sum(1 for w in td_worlds if w['sustained_oscillation'])

                rows.append({
                    'level': 'topology_x_delta_x_strength',
                    'topology_type': topo,
                    'delta': d,
                    'strength_regime': sr,
                    'N_world': n,
                    'fraction_convergent': reg_counter.get('Convergent', 0) / n,
                    'fraction_divergent': reg_counter.get('Divergent', 0) / n,
                    'fraction_collapse': reg_counter.get('Collapse', 0) / n,
                    'fraction_sustained_oscillatory': reg_counter.get('Sustained oscillatory', 0) / n,
                    'clipping_dominated_fraction': clipping_frac,
                    'sustained_oscillation_count': sustained_count,
                    'first_divergent_regime': '',
                })

    return rows


def write_tables(all_results: List[Dict], aggregation: Dict) -> None:
    """Write all TSV tables and analysis JSONs.

    Per Output Spec §11: world_summary, regime_summary, boundary_summary.
    Also saves the full aggregation dictionary as JSON.
    """
    print("\n── Step 5: Writing Tables ──")
    os.makedirs(TABLES_DIR, exist_ok=True)

    # 11.1 World Summary Table
    world_rows = build_world_summary_table(all_results)
    _write_tsv(world_rows, os.path.join(TABLES_DIR, 'world_summary.tsv'))
    print(f"  world_summary.tsv: {len(world_rows)} rows")

    # 11.2 Regime Summary Table
    regime_rows = build_regime_summary_table(aggregation, all_results)
    _write_tsv(regime_rows, os.path.join(TABLES_DIR, 'regime_summary.tsv'))
    print(f"  regime_summary.tsv: {len(regime_rows)} rows")

    # 11.3 Boundary Summary Table
    boundary_rows = build_boundary_summary_table(all_results)
    _write_tsv(boundary_rows, os.path.join(TABLES_DIR, 'boundary_summary.tsv'))
    print(f"  boundary_summary.tsv: {len(boundary_rows)} rows")

    # Full aggregation JSON
    agg_save = {
        'combo_stats': {},
        **{k: v for k, v in aggregation.items()},
    }
    # Make JSON-serializable
    def _convert(obj):
        if isinstance(obj, dict):
            return {str(k): _convert(v) for k, v in obj.items()}
        if isinstance(obj, list):
            return [_convert(v) for v in obj]
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return obj

    agg_json = _convert(agg_save)
    with open(os.path.join(TABLES_DIR, 'aggregation.json'), 'w') as f:
        json.dump(agg_json, f, indent=2)
    print(f"  aggregation.json saved.")


def _write_tsv(rows: List[Dict], path: str) -> None:
    """Write a list of dicts as a TSV file."""
    if not rows:
        return
    fieldnames = list(rows[0].keys())
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in rows:
            # Convert non-serializable values
            clean = {}
            for k, v in row.items():
                if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
                    clean[k] = None
                else:
                    clean[k] = v
            writer.writerow(clean)


# ═══════════════════════════════════════════════════════════════════
# Step 6:  Figures
# ═══════════════════════════════════════════════════════════════════

def _assemble_world_summaries(all_results: List[Dict]) -> List[Dict]:
    """Build flat world-level summary rows for figure generation.

    Each row = aggregate_world + combo metadata.
    Cached for reuse across all plot functions.
    """
    rows = []
    for world_result in all_results:
        wa = aggregate_world(world_result['per_cell'])
        rows.append({
            'world_id':         world_result['world_id'],
            'combo_key':        world_result['combo_key'],
            'topology_type':    world_result['topology_type'],
            'delta':            world_result['delta'],
            'strength_regime':  world_result['strength_regime'],
            'sign_ratio':       world_result['sign_ratio'],
            'primary_regime':   wa['primary_regime'],
            'secondary_labels': wa['secondary_labels'],
            'convergence_time': wa['convergence_time'],
            'sustained_oscillation': wa['sustained_oscillation'],
            'divergence':       wa['divergence'],
            'collapse':         wa['collapse'],
            'clipping_frequency': wa['clipping_frequency'],
            'clipping_dominated': wa['clipping_dominated'],
            'max_amplitude':    wa['max_amplitude'],
            'final_amplitude':  wa['final_amplitude'],
            'spectral_radius':  wa['spectral_radius'],
            'equilibrium_magnitude': wa['equilibrium_magnitude'],
            'equilibrium_sparsity': wa['equilibrium_sparsity'],
        })
    return rows


# ── §12.1: Primary Regime Heatmaps ──────────────────────────────

def plot_primary_regime_heatmaps(all_results: List[Dict]) -> None:
    """Per Output Spec §12.1: δ × strength heatmap, 2 topology × 2 sign = 4 subplots."""
    ws = _assemble_world_summaries(all_results)

    regime_colors = {
        'Convergent':           '#2ecc71',
        'Sustained oscillatory':'#e74c3c',
        'Divergent':            '#9b59b6',
        'Collapse':             '#95a5a6',
        'Ambiguous':            '#bdc3c7',
    }

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    x_labels = [f'{sr}\na ∈ [{STRENGTH_REGIMES[sr][0]}, {STRENGTH_REGIMES[sr][1]}]'
                for sr in STRENGTH_ORDER]
    y_labels = [f'{d}' for d in DELTA_VALUES]

    n_rows = len(DELTA_VALUES)
    n_cols = len(STRENGTH_ORDER)

    for row, topo in enumerate(TOPOLOGY_TYPES):
        for col, sign_ratio in enumerate(SIGN_RATIOS):
            ax = axes[row, col]
            img = np.full((n_rows, n_cols, 3), 0.8)

            for ri, d in enumerate(DELTA_VALUES):
                for ci, sr in enumerate(STRENGTH_ORDER):
                    subset = [w for w in ws
                              if w['topology_type'] == topo
                              and w['delta'] == d
                              and w['strength_regime'] == sr
                              and w['sign_ratio'] == sign_ratio]
                    if not subset:
                        continue
                    regimes = [w['primary_regime'] for w in subset]
                    regime_counts = Counter(regimes)
                    total = len(subset)
                    dominant = regime_counts.most_common(1)[0][0]
                    frac = regime_counts[dominant] / total
                    c = regime_colors.get(dominant, '#bdc3c7')
                    hex_rgb = np.array([int(c[i+1:i+3], 16) / 255.0 for i in (0, 2, 4)])
                    # Blend with white by fraction → continuous color intensity
                    white = np.array([1.0, 1.0, 1.0])
                    blended = white + (hex_rgb - white) * frac
                    img[ri, ci] = np.clip(blended, 0, 1)

                    # Show all regimes present in this cell, one per line
                    # e.g. "Conv 80%\nDive 20%"; smaller font when 2+ regimes
                    items = regime_counts.most_common()
                    if len(items) == 1:
                        label = f'{items[0][0][:4]} {items[0][1]/total:.0%}'
                        fs = 10
                    else:
                        label = '\n'.join(
                            f'{r[:4]} {n/total:.0%}' for r, n in items
                        )
                        fs = 8
                    ax.text(ci, ri, label,
                            ha='center', va='center', fontsize=fs, fontweight='bold')

            ax.imshow(img, aspect='equal', origin='upper')
            ax.set_xticks(range(n_cols))
            # Only bottom row shows x-tick labels
            if row == 1:
                ax.set_xticklabels(x_labels, fontsize=10)
            else:
                ax.set_xticklabels([''] * n_cols)
            ax.set_yticks(range(n_rows))
            # Only left column shows y-tick labels
            if col == 0:
                ax.set_yticklabels(y_labels, fontsize=11)
            else:
                ax.set_yticklabels([''] * n_rows)
            ax.set_title(f'{sign_ratio} ({SIGN_RATIO_RATIOS.get(sign_ratio, "")})',
                         fontsize=12)
            if col == 0:
                ax.set_ylabel('δ', fontsize=13)
            if row == 1:
                ax.set_xlabel('strength regime', fontsize=11)

    fig.suptitle('Primary Regime Maps',
                 fontsize=16, fontweight='bold')

    # Layout: generous hspace so x-tick labels, row titles, and sub-titles don't overlap
    fig.subplots_adjust(left=0.08, right=0.98, top=0.90, bottom=0.06,
                        hspace=0.36, wspace=0.06)

    # Row titles placed dynamically from subplot top edges (must be after subplots_adjust)
    row_labels = ['Unrestricted Sparse GRN', 'TF Restricted GRN']
    for row_idx, label in enumerate(row_labels):
        top_y = axes[row_idx, 0].get_position().y1
        fig.text(0.5, top_y + 0.035, label,
                 fontsize=13, fontweight='bold', ha='center', va='bottom')

    fig.savefig(os.path.join(FIGURES_DIR, '01_primary_regime_heatmaps.png'),
                dpi=150, bbox_inches='tight')
    plt.close(fig)


# ── §12.2: Topology × δ Comparison ───────────────────────────────

def plot_topology_delta_comparison(all_results: List[Dict]) -> None:
    """Per Output Spec §12.2: unrestricted vs tf_restricted across δ.

    Four panels: Convergent fraction, Divergent fraction, mean clipping_frequency, mean ρ(A).
    """
    ws = _assemble_world_summaries(all_results)
    topo_colors = {'unrestricted_sparse': '#e74c3c', 'tf_restricted': '#3498db'}

    metrics = [
        ('fraction_convergent', 'Convergent Fraction'),
        ('fraction_divergent',  'Divergent Fraction'),
        ('mean_clipping_frequency', 'Mean Clipping Frequency'),
        ('spectral_radius',     'Mean Spectral Radius ρ(A)'),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))
    axes_flat = axes.flatten()

    for ai, (metric_key, metric_label) in enumerate(metrics):
        ax = axes_flat[ai]
        x = np.arange(len(DELTA_VALUES))
        width = 0.35

        # Compute all values first (to determine y-axis headroom)
        all_vals: Dict[str, List[float]] = {}
        for topo in TOPOLOGY_TYPES:
            vals = []
            for d in DELTA_VALUES:
                subset = [w for w in ws
                          if w['topology_type'] == topo and w['delta'] == d]
                if metric_key == 'fraction_convergent':
                    v = sum(1 for w in subset if w['primary_regime'] == 'Convergent') / len(subset)
                elif metric_key == 'fraction_divergent':
                    v = sum(1 for w in subset if w['primary_regime'] == 'Divergent') / len(subset)
                elif metric_key == 'mean_clipping_frequency':
                    v = np.mean([w['clipping_frequency'] for w in subset])
                else:
                    v = np.mean([w['spectral_radius'] for w in subset])
                vals.append(v)
            all_vals[topo] = vals

        # Y-axis limits with headroom for text labels
        global_max = max(max(v) for v in all_vals.values())
        if metric_key in ('fraction_convergent', 'fraction_divergent'):
            ax.set_ylim(0, 1.1)
        else:
            ax.set_ylim(0, global_max * 1.25)

        # Draw bars + text
        text_offs = global_max * 0.02
        for ti, topo in enumerate(TOPOLOGY_TYPES):
            vals = all_vals[topo]
            pos = x + (ti - 0.5) * width
            ax.bar(pos, vals, width, color=topo_colors[topo],
                   label=topo, edgecolor='white', linewidth=0.5, alpha=0.7)
            for px, vy in zip(pos, vals):
                if metric_key in ('fraction_convergent', 'fraction_divergent', 'mean_clipping_frequency'):
                    ax.text(px, vy + text_offs, f'{vy:.1%}', ha='center', fontsize=7.5)
                else:
                    ax.text(px, vy + text_offs, f'{vy:.2f}', ha='center', fontsize=7.5)

        ax.set_xticks(x)
        ax.set_xticklabels([f'δ={d}' for d in DELTA_VALUES], fontsize=9)
        ax.set_title(metric_label, fontsize=11)
        ax.grid(axis='y', linestyle='--', alpha=0.3)

    fig.suptitle('Topology × δ Comparison — unrestricted_sparse vs tf_restricted',
                 fontsize=13, fontweight='bold', y=0.98)
    fig.legend(
        [plt.Rectangle((0,0),1,1,facecolor=topo_colors[t],alpha=0.7) for t in TOPOLOGY_TYPES],
        TOPOLOGY_TYPES, fontsize=10, loc='upper center', ncol=2,
        bbox_to_anchor=(0.5, 0.95))
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(FIGURES_DIR, '02_topology_delta_comparison.png'), dpi=150)
    plt.close(fig)


# ── §12.3: Clipping Boundary ─────────────────────────────────────

def plot_clipping_boundary(all_results: List[Dict]) -> None:
    """Per Output Spec §12.3: δ vs clipping_frequency, 3 subplots.

    Subplot 1: grouped by topology_type
    Subplot 2: grouped by strength_regime
    Subplot 3: grouped by sign_ratio
    All subplots: x = δ, y = clipping_frequency, boxplot per cross-cell.
    """
    ws = _assemble_world_summaries(all_results)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    ax1, ax2, ax3 = axes

    # ── shared setup ──
    n_deltas = len(DELTA_VALUES)
    base_colors = {
        'unrestricted_sparse':   '#e74c3c',
        'tf_restricted':         '#3498db',
        'baseline':              '#2ecc71',
        'chen_moderate':         '#3498db',
        'stress':                '#e74c3c',
        'chen_stress':           '#9b59b6',
        'balanced':              '#2ecc71',
        'repression_biased':     '#e74c3c',
    }

    # ── Subplot 1: topology_type ──
    topo_categories = TOPOLOGY_TYPES
    n_cats = len(topo_categories)
    group_w1 = n_cats + 0.8
    for di, d in enumerate(DELTA_VALUES):
        x0 = di * group_w1
        for ci, cat in enumerate(topo_categories):
            vals = [w['clipping_frequency'] for w in ws
                    if w['delta'] == d and w['topology_type'] == cat]
            bp = ax1.boxplot(vals, positions=[x0 + ci], widths=0.7,
                             patch_artist=True)
            for box in bp['boxes']:
                box.set_facecolor(base_colors[cat])
                box.set_alpha(0.7)
    ticks1 = [di * group_w1 + (n_cats - 1) / 2 for di in range(n_deltas)]
    ax1.set_xticks(ticks1)
    ax1.set_xticklabels([f'δ={d}' for d in DELTA_VALUES], fontsize=9)
    ax1.set_ylabel('clipping frequency', fontsize=9)
    ax1.set_title('By topology_type', fontsize=13)
    ax1.axhline(y=CLIPPING_FRAC_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax1.grid(axis='y', linestyle='--', alpha=0.3)
    ax1.legend(
        [plt.Rectangle((0,0),1,1,facecolor=base_colors[c],alpha=0.7) for c in topo_categories],
        topo_categories, fontsize=8, loc='upper right')

    # ── Subplot 2: strength_regime ──
    sr_categories = STRENGTH_ORDER
    n_cats = len(sr_categories)
    group_w2 = n_cats + 0.8
    for di, d in enumerate(DELTA_VALUES):
        x0 = di * group_w2
        for ci, cat in enumerate(sr_categories):
            vals = [w['clipping_frequency'] for w in ws
                    if w['delta'] == d and w['strength_regime'] == cat]
            bp = ax2.boxplot(vals, positions=[x0 + ci], widths=0.7,
                             patch_artist=True)
            for box in bp['boxes']:
                box.set_facecolor(base_colors[cat])
                box.set_alpha(0.7)
    ticks2 = [di * group_w2 + (n_cats - 1) / 2 for di in range(n_deltas)]
    ax2.set_xticks(ticks2)
    ax2.set_xticklabels([f'δ={d}' for d in DELTA_VALUES], fontsize=9)
    ax2.set_ylabel('clipping frequency', fontsize=9)
    ax2.set_title('By strength_regime', fontsize=13)
    ax2.axhline(y=CLIPPING_FRAC_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax2.grid(axis='y', linestyle='--', alpha=0.3)
    ax2.legend(
        [plt.Rectangle((0,0),1,1,facecolor=base_colors[c],alpha=0.7) for c in sr_categories],
        sr_categories, fontsize=8, loc='upper right')

    # ── Subplot 3: sign_ratio ──
    sign_categories = list(SIGN_RATIOS.keys())
    n_cats = len(sign_categories)
    group_w3 = n_cats + 0.8
    for di, d in enumerate(DELTA_VALUES):
        x0 = di * group_w3
        for ci, cat in enumerate(sign_categories):
            vals = [w['clipping_frequency'] for w in ws
                    if w['delta'] == d and w['sign_ratio'] == cat]
            bp = ax3.boxplot(vals, positions=[x0 + ci], widths=0.7,
                             patch_artist=True)
            for box in bp['boxes']:
                box.set_facecolor(base_colors[cat])
                box.set_alpha(0.7)
    ticks3 = [di * group_w3 + (n_cats - 1) / 2 for di in range(n_deltas)]
    ax3.set_xticks(ticks3)
    ax3.set_xticklabels([f'δ={d}' for d in DELTA_VALUES], fontsize=10)
    ax3.set_ylabel('clipping frequency', fontsize=10)
    ax3.set_title('By sign_ratio', fontsize=13)
    ax3.axhline(y=CLIPPING_FRAC_THRESHOLD, color='red', linestyle='--', alpha=0.5)
    ax3.grid(axis='y', linestyle='--', alpha=0.3)
    ax3.legend(
        [plt.Rectangle((0,0),1,1,facecolor=base_colors[c],alpha=0.7) for c in sign_categories],
        sign_categories, fontsize=8, loc='upper right')

    fig.suptitle('Clipping Boundary — δ vs Clipping Frequency',
                 fontsize=13, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(os.path.join(FIGURES_DIR, '03_clipping_boundary.png'), dpi=150)
    plt.close(fig)


# ── §12.4: Spectral Radius ───────────────────────────────────────

def plot_spectral_radius(all_results: List[Dict]) -> None:
    """Per Output Spec §12.4: ρ(A) vs primary_regime (boxplot), matched to run_002/003 style.

    Also includes: ρ(A) vs clipping_frequency scatter, ρ(A) across δ.
    """
    ws = _assemble_world_summaries(all_results)

    # ── 12.4a: ρ(A) by primary_regime (boxplot, run_002/003 style) ──
    regime_order = PRIMARY_REGIME_ORDER
    by_regime = defaultdict(list)
    by_regime_no_clip = defaultdict(list)
    by_regime_clip_count = defaultdict(int)

    for w in ws:
        pr = w['primary_regime']
        by_regime[pr].append(w['spectral_radius'])
        if not w['clipping_dominated']:
            by_regime_no_clip[pr].append(w['spectral_radius'])
        if w['clipping_dominated']:
            by_regime_clip_count[pr] += 1

    labels, data, label_text = [], [], []
    for pr in regime_order:
        if by_regime[pr]:
            n_total = len(by_regime[pr])
            n_clip = by_regime_clip_count[pr]
            labels.append(pr)
            data.append(by_regime[pr])
            label_text.append(f"{pr}\n(n={n_total}, clip={n_clip}/{n_total})")

    if data:
        fig, ax = plt.subplots(figsize=(10, 6))
        pos = list(range(1, len(labels) + 1))
        bp = ax.boxplot(data, positions=pos, tick_labels=label_text,
                        patch_artist=True, widths=0.6)
        for patch in bp['boxes']:
            patch.set_facecolor('#a8c5e8')
        ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.7, label='ρ=1')
        ax.set_xticklabels(label_text, fontsize=8)
        ax.set_ylabel('Spectral Radius ρ(A)', fontsize=11)
        ax.set_title('ρ(A) by Primary Regime\n'
                     '(n = world count, clip = clipping-dominated worlds)')
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR, '04_spectral_radius_by_regime.png'), dpi=150)
        plt.close(fig)

    # ── 12.4b: ρ(A) vs clipping_frequency ──
    topo_markers = {'unrestricted_sparse': 'o', 'tf_restricted': '^'}
    fig, ax = plt.subplots(figsize=(8, 5))
    for topo in TOPOLOGY_TYPES:
        topo_ws = [w for w in ws if w['topology_type'] == topo]
        if not topo_ws:
            continue
        cfs = [w['clipping_frequency'] for w in topo_ws]
        srs = [w['spectral_radius'] for w in topo_ws]
        ax.scatter(cfs, srs, marker=topo_markers.get(topo, 'o'), alpha=0.5,
                   label=topo, s=20)

    ax.set_xlabel('clipping frequency', fontsize=11)
    ax.set_ylabel('spectral radius ρ(A)', fontsize=11)
    ax.set_title('ρ(A) vs Clipping Frequency')
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='ρ=1')
    ax.legend(fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, '04_spectral_radius_vs_clipping.png'), dpi=150)
    plt.close(fig)

    # ── 12.4c: ρ(A) across δ ── (boxplot, δ as categorical)
    fig, ax = plt.subplots(figsize=(9, 5))
    topo_colors = {'unrestricted_sparse': '#3498db', 'tf_restricted': '#e67e22'}
    positions = np.arange(len(DELTA_VALUES))
    width = 0.35
    for i, topo in enumerate(TOPOLOGY_TYPES):
        data_per_d = []
        for d in DELTA_VALUES:
            vals = [w['spectral_radius'] for w in ws
                    if w['topology_type'] == topo and w['delta'] == d]
            data_per_d.append(vals if vals else [np.nan])
        offset = (i - 0.5) * width
        bp = ax.boxplot(data_per_d, positions=positions + offset, widths=width * 0.9,
                        patch_artist=True, showfliers=False,
                        boxprops=dict(facecolor=topo_colors[topo], alpha=0.7,
                                      edgecolor='black'),
                        medianprops=dict(color='black', linewidth=1.5),
                        whiskerprops=dict(color='black'),
                        capprops=dict(color='black'))
        ax.plot([], [], 's', color=topo_colors[topo], label=topo, markersize=8)

    ax.set_xticks(positions)
    ax.set_xticklabels([f'{d}' for d in DELTA_VALUES], fontsize=11)
    ax.set_xlabel('δ', fontsize=11)
    ax.set_ylabel('Spectral Radius ρ(A)', fontsize=11)
    ax.set_title('ρ(A) across δ', fontsize=13)
    ax.axhline(y=1.0, color='red', linestyle='--', alpha=0.5, label='ρ=1')
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, linestyle='--', alpha=0.3, axis='y')
    plt.tight_layout()
    fig.savefig(os.path.join(FIGURES_DIR, '04_spectral_radius_by_delta.png'), dpi=150)
    plt.close(fig)


# ── §12.5: Trajectory Exemplars ──────────────────────────────────

def _get_max_step(regime: str, X_len: int) -> int:
    """Return the max step index to plot based on primary regime.

    Matched to run_002/003 convention for zoom behavior.
    """
    if regime == 'Convergent':
        return min(51, X_len)
    else:
        return X_len


def plot_trajectory_exemplars(all_results: List[Dict]) -> None:
    """Per Output Spec §12.5: exemplar trajectories per primary regime.

    All regimes in a single 2-column figure: full horizon (left) + early zoom (right).
    For tf_restricted worlds, TF genes (0-4) shown in red, non-TF in grey.
    """

    # Find one exemplar world per primary regime (first occurrence)
    exemplars: Dict[str, str] = {}
    for world_result in all_results:
        wa = aggregate_world(world_result['per_cell'])
        pr = wa['primary_regime']
        if pr not in exemplars:
            exemplars[pr] = world_result['world_id']

    if not exemplars:
        return

    # Build world_id → topology_type lookup
    topo_map: Dict[str, str] = {r['world_id']: r['topology_type'] for r in all_results}

    # Order by PRIMARY_REGIME_ORDER, keep only those with exemplars
    regime_order = [r for r in PRIMARY_REGIME_ORDER if r in exemplars]
    n_regimes = len(regime_order)
    if n_regimes == 0:
        return

    fig, axes = plt.subplots(n_regimes, 2, figsize=(14, 4 * n_regimes))
    if n_regimes == 1:
        axes = np.array([axes])

    for row, regime in enumerate(regime_order):
        wid = exemplars[regime]
        traj_files = sorted(
            [f for f in os.listdir(TRAJECTORIES_DIR)
             if f.startswith(wid) and f.endswith('.pt')]
        )
        if not traj_files:
            continue

        data = torch.load(os.path.join(TRAJECTORIES_DIR, traj_files[0]),
                          map_location='cpu', weights_only=False)
        X = data['X_traj'].numpy()
        T_total = X.shape[0]
        is_tf_restricted = (topo_map.get(wid, '') == 'tf_restricted')

        ax_full = axes[row, 0]
        ax_zoom = axes[row, 1]

        # Full horizon
        for g in range(G):
            if is_tf_restricted and g in TF_GENES:
                ax_full.plot(X[:, g], color='red', alpha=0.8, linewidth=1.0)
            elif is_tf_restricted:
                ax_full.plot(X[:, g], color='grey', alpha=0.35, linewidth=0.8)
            else:
                ax_full.plot(X[:, g], alpha=0.9, linewidth=0.8)
        ax_full.set_title(f'{regime}\n{wid}', fontsize=12, linespacing=1.3)
        ax_full.set_xlabel('t', fontsize=11)
        ax_full.set_ylabel('gene activity', fontsize=11)
        ax_full.spines['top'].set_visible(False)
        ax_full.spines['right'].set_visible(False)
        ax_full.tick_params(labelsize=9)

        # Early zoom [0, 200]
        zoom_end = min(200, T_total)
        for g in range(G):
            if is_tf_restricted and g in TF_GENES:
                ax_zoom.plot(X[:zoom_end, g], color='red', alpha=0.8, linewidth=1.0)
            elif is_tf_restricted:
                ax_zoom.plot(X[:zoom_end, g], color='grey', alpha=0.35, linewidth=0.8)
            else:
                ax_zoom.plot(X[:zoom_end, g], alpha=0.9, linewidth=0.8)
        ax_zoom.set_title(f'Early Dynamics Zoom: t = 0 ... {zoom_end}', fontsize=12)
        ax_zoom.set_xlabel('t', fontsize=11)
        ax_zoom.set_ylabel('gene activity', fontsize=11)
        ax_zoom.spines['top'].set_visible(False)
        ax_zoom.spines['right'].set_visible(False)
        ax_zoom.tick_params(labelsize=9)

        # Legend inside subplot for tf_restricted worlds
        if is_tf_restricted:
            legend_handles = [
                plt.Line2D([0], [0], color='red', lw=1.5, label='TF'),
                plt.Line2D([0], [0], color='grey', lw=1.5, label='non-TF'),
            ]
            ax_full.legend(handles=legend_handles, fontsize=10,
                           loc='upper right', framealpha=0.7)
            ax_zoom.legend(handles=legend_handles, fontsize=10,
                           loc='upper right', framealpha=0.7)

    fig.suptitle('Trajectory Exemplars by Primary Regime', fontsize=16, fontweight='bold', y=0.99)
    plt.tight_layout(rect=[0, 0, 1, 0.99])
    fig.savefig(os.path.join(FIGURES_DIR, '05_trajectory_exemplars.png'), dpi=150)
    plt.close(fig)


# ── Sign Ratio Comparison ────────────────────────────────────────

def plot_sign_ratio_comparison(all_results: List[Dict]) -> None:
    """Stacked bar: regime distribution per (topology_type, δ, strength_regime, sign_ratio).

    Balanced vs repression_biased across δ × strength, faceted by topology_type.
    """
    ws = _assemble_world_summaries(all_results)

    regime_colors = {
        'Convergent':           '#2ecc71',
        'Sustained oscillatory':'#e74c3c',
        'Divergent':            '#9b59b6',
        'Collapse':             '#95a5a6',
        'Ambiguous':            '#bdc3c7',
    }

    for topo in TOPOLOGY_TYPES:
        n_combos = len(DELTA_VALUES) * len(STRENGTH_ORDER) * 2  # 2 sign ratios
        if n_combos < 1:
            continue

        fig, axes = plt.subplots(len(DELTA_VALUES), 1,
                                 figsize=(12, 4 * len(DELTA_VALUES)),
                                 squeeze=False)

        for di, d in enumerate(DELTA_VALUES):
            ax = axes[di][0]
            # Build per-(strength, sign) regime counts
            groups = []
            labels = []
            for sr in STRENGTH_ORDER:
                for sign_ratio in SIGN_RATIOS:
                    subset = [w for w in ws
                              if w['topology_type'] == topo
                              and w['delta'] == d
                              and w['strength_regime'] == sr
                              and w['sign_ratio'] == sign_ratio]
                    counts = Counter(w['primary_regime'] for w in subset)
                    groups.append(counts)
                    labels.append(f'{sr}\n{sign_ratio[:5]}')

            # Stacked bar
            x = np.arange(len(groups))
            bottoms = np.zeros(len(groups))
            for regime in PRIMARY_REGIME_ORDER:
                vals = np.array([g.get(regime, 0) for g in groups], dtype=float)
                ax.bar(x, vals, bottom=bottoms,
                       color=regime_colors.get(regime, '#bdc3c7'),
                       label=regime, edgecolor='white', linewidth=0.5)
                bottoms += vals

            ax.set_xticks(x)
            ax.set_xticklabels(labels, fontsize=8)
            ax.set_ylabel('World Count', fontsize=10)
            ax.set_title(f'δ = {d}', fontsize=11)
            ax.set_ylim(0, N_WORLD * 1.15)
            if di == 0:
                ax.legend(fontsize=7, loc='upper right', ncol=3)

        fig.suptitle(f'Sign Ratio Comparison — {topo}\n'
                     f'(regime distribution per δ × strength × sign_ratio)',
                     fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig.savefig(os.path.join(FIGURES_DIR,
                    f'06_sign_ratio_comparison_{topo}.png'), dpi=150,
                    bbox_inches='tight')
        plt.close(fig)


# ── All Worlds Trajectories (run_002/003 style grid) ─────────────

def plot_all_worlds_trajectories(all_results: List[Dict]) -> None:
    """Grid of all worlds' trajectories (like run_002 / run_003).

    For each (topology_type, delta, sign_ratio):
      rows = strength_regimes (4)
      cols = world replicates (MAX_COLS per batch)
    Two batches per combo to cover all 10 worlds (batch1: ti 0-4, batch2: ti 5-9).
    """
    traj_fig_dir = os.path.join(FIGURES_DIR, 'traj')
    os.makedirs(traj_fig_dir, exist_ok=True)

    MAX_COLS = 5  # worlds per batch

    for topo in TOPOLOGY_TYPES:
        for d in DELTA_VALUES:
            for sign_ratio in SIGN_RATIOS:
                # Find all worlds for this combo slice
                slice_worlds = [r for r in all_results
                                if r['topology_type'] == topo
                                and r['delta'] == d
                                and r['sign_ratio'] == sign_ratio]
                if not slice_worlds:
                    continue

                # Group by strength_regime → list of world_ids (sorted)
                by_strength: Dict[str, List[str]] = defaultdict(list)
                for r in slice_worlds:
                    by_strength[r['strength_regime']].append(r['world_id'])
                for sr in by_strength:
                    by_strength[sr].sort()

                strength_order = [sr for sr in STRENGTH_ORDER if sr in by_strength]
                n_rows = len(strength_order)
                if n_rows == 0:
                    continue

                # Number of batches needed (ceil(max_worlds / MAX_COLS))
                max_worlds = max(len(by_strength[sr]) for sr in strength_order)
                n_batches = (max_worlds + MAX_COLS - 1) // MAX_COLS

                dl = delta_label(d)

                for batch in range(n_batches):
                    start = batch * MAX_COLS
                    end = start + MAX_COLS

                    # Worlds in this batch per strength regime
                    batch_worlds: Dict[str, List[str]] = {
                        sr: by_strength[sr][start:end] for sr in strength_order
                    }
                    batch_n_cols = MAX_COLS

                    fig, axes = plt.subplots(n_rows, batch_n_cols,
                                             figsize=(batch_n_cols * 4.5, n_rows * 2.8),
                                             squeeze=False)

                    for ri, sr in enumerate(strength_order):
                        world_ids = batch_worlds[sr]
                        for ci in range(batch_n_cols):
                            ax = axes[ri, ci]
                            if ci >= len(world_ids):
                                ax.set_visible(False)
                                continue

                            wid = world_ids[ci]
                            traj_files = sorted(
                                [f for f in os.listdir(TRAJECTORIES_DIR)
                                 if f.startswith(wid) and f.endswith('.pt')]
                            )
                            if not traj_files:
                                ax.set_title(f'{wid} — no data', fontsize=9)
                                continue

                            data = torch.load(os.path.join(TRAJECTORIES_DIR,
                                              traj_files[0]),
                                              map_location='cpu', weights_only=False)
                            X = data['X_traj'].numpy()
                            wa = aggregate_world(next(
                                r['per_cell'] for r in all_results
                                if r['world_id'] == wid))
                            pr = wa['primary_regime']

                            for g in range(G):
                                if topo == 'tf_restricted' and g in TF_GENES:
                                    ax.plot(X[:, g], color='red', alpha=0.8, linewidth=1.0)
                                elif topo == 'tf_restricted':
                                    ax.plot(X[:, g], color='grey', alpha=0.35, linewidth=0.8)
                                else:
                                    ax.plot(X[:, g], alpha=0.9, linewidth=1.2)
                            ax.set_title(f'[{pr}] {wid}', fontsize=9)
                            ax.tick_params(labelsize=7)
                            ax.spines['top'].set_visible(False)
                            ax.spines['right'].set_visible(False)

                    for ri, sr in enumerate(strength_order):
                        a_min, a_max, _ = STRENGTH_REGIMES[sr]
                        axes[ri, 0].set_ylabel(
                            f'{sr}\na [{a_min}, {a_max}]',
                            fontsize=9, rotation=0, ha='right', va='center',
                            labelpad=10)

                    batch_label = f'_batch{batch + 1}' if n_batches > 1 else ''
                    fig.suptitle(
                        f'World Trajectories (Cell 0) — {topo} / ={d} / {sign_ratio}'
                        f'{batch_label}',
                        fontsize=12, fontweight='bold')
                    if topo == 'tf_restricted':
                        legend_handles = [
                            plt.Line2D([0], [0], color='red', lw=2, label='TF'),
                            plt.Line2D([0], [0], color='grey', lw=2, label='non-TF'),
                        ]
                        fig.legend(handles=legend_handles, fontsize=9,
                                   loc='lower center', bbox_to_anchor=(0.5, -0.01), ncol=2)
                    plt.tight_layout()
                    fig.savefig(os.path.join(traj_fig_dir,
                                f'all_trajectories_{topo}_{dl}_{sign_ratio}{batch_label}.png'),
                                dpi=150)
                    plt.close(fig)

    print("  all_worlds_trajectories complete.")


# ── Master Figure Generator ──────────────────────────────────────

def generate_all_figures(all_results: List[Dict]) -> None:
    """Master figure generation: calls all plot_* functions."""
    print("\n── Step 6: Generating Figures ──")

    print("  §12.1 Primary regime heatmaps ...")
    plot_primary_regime_heatmaps(all_results)

    print("  §12.2 Topology × δ comparison ...")
    plot_topology_delta_comparison(all_results)

    print("  §12.3 Clipping boundary ...")
    plot_clipping_boundary(all_results)

    print("  §12.4 Spectral radius ...")
    plot_spectral_radius(all_results)

    print("  §12.5 Trajectory exemplars ...")
    plot_trajectory_exemplars(all_results)

    print("  Sign ratio comparison ...")
    plot_sign_ratio_comparison(all_results)

    print("  All worlds trajectories ...")
    plot_all_worlds_trajectories(all_results)

    print(f"\n── Figures complete. ──")


# ═══════════════════════════════════════════════════════════════════
# Step 7:  Perturbation  (stubs)
# ═══════════════════════════════════════════════════════════════════

def select_perturbation_candidates(all_results: List[Dict]) -> List[Dict]:
    """Select representative worlds for canonical KD perturbation.

    Strategy: cover each existing (primary_regime, delta) combination with one world,
    balancing topology_type across selections (aligned with run_002 / run_003 principle
    of per-regime representative selection).

    Returns list of dicts with keys: world_id, cell_idx, topology_type, primary_regime,
    delta, strength_regime, sign_ratio, combo_key.
    """
    # Group worlds by (primary_regime, delta)
    pool: Dict[Tuple[str, float], List[Dict]] = defaultdict(list)
    for r in all_results:
        wa = aggregate_world(r['per_cell'])
        key = (wa['primary_regime'], r['delta'])
        pool[key].append({
            'world_id': r['world_id'],
            'cell_idx': 0,
            'topology_type': r['topology_type'],
            'combo_key': r['combo_key'],
            'delta': r['delta'],
            'strength_regime': r['strength_regime'],
            'sign_ratio': r['sign_ratio'],
            'primary_regime': wa['primary_regime'],
        })

    # Pick one per (regime, delta), balancing topology coverage greedily
    candidates: List[Dict] = []
    topo_counts: Counter = Counter()
    for key in sorted(pool.keys()):
        worlds = pool[key]
        best = min(worlds, key=lambda w: topo_counts.get(w['topology_type'], 0))
        candidates.append(best)
        topo_counts[best['topology_type']] += 1

    n_regimes = len(set(c['primary_regime'] for c in candidates))
    print(f"  selected {len(candidates)} perturbation candidates "
          f"({n_regimes} regimes, across {len(pool)} (regime,δ) groups)")
    return candidates


def run_perturbation_analysis(all_results: List[Dict]) -> None:
    """Run canonical KD perturbation for each candidate (per Output Spec §10).

    Protocol: t=500, KD the highest-expression gene at that time, one-step to 0.
    Saves perturbed trajectories and perturbation_summary.json.
    """
    candidates = select_perturbation_candidates(all_results)
    if not candidates:
        print("  skip perturbation: no candidates")
        return

    os.makedirs(PERTURBATIONS_DIR, exist_ok=True)
    summary: List[Dict] = []

    for ci, cand in enumerate(candidates):
        wid = cand['world_id']
        topo = cand['topology_type']
        cell_idx = cand['cell_idx']
        cell_label = f'cell{cell_idx:02d}'
        original_regime = cand['primary_regime']

        print(f"\n  [{ci+1}/{len(candidates)}] perturb: {wid} / {cell_label} "
              f"({original_regime})")

        # Load world metadata and reconstruct World
        meta_path = os.path.join(METADATA_DIR, f'{wid}.json')
        if not os.path.exists(meta_path):
            print(f"    skip: metadata not found")
            continue
        with open(meta_path) as f:
            meta = json.load(f)

        dma = get_model_module(topo)
        # from_dict requires 'seed' key; metadata strips it, but world_seed holds the value
        meta['seed'] = meta['world_seed']
        w = dma.World(meta['world_seed'])
        w.from_dict(meta)

        # Load original trajectory at t=PERTURBATION_TIME
        traj_path = os.path.join(TRAJECTORIES_DIR, f'{wid}_{cell_label}.pt')
        if not os.path.exists(traj_path):
            print(f"    skip: trajectory not found")
            continue
        orig_traj = torch.load(traj_path, map_location='cpu', weights_only=False)
        X_at_pert = orig_traj['X_traj'][PERTURBATION_TIME]
        kd_gene = int(torch.argmax(X_at_pert).item())
        max_expr = float(X_at_pert[kd_gene].item())
        kd_gene_role = 'TF' if topo == 'tf_restricted' and kd_gene in TF_GENES else 'non-TF'
        print(f"    KD gene {kd_gene} ({kd_gene_role}), expr={max_expr:.4f} at t={PERTURBATION_TIME}")

        # Re-run simulation with KD at t=500
        cell_seed = compute_cell_seed(orig_traj['world_seed'], cell_idx)
        X0 = dma.sample_initial_state(cell_seed)
        result = dma.simulate_single_cell(
            w, X0, T_SIM,
            intervention_time=PERTURBATION_TIME,
            intervention_config={'knockdown': [kd_gene]},
        )

        X_traj_full = result['X_traj']
        clip_count_full = result['clip_count']

        # ── Post-perturbation analysis ──
        X_post = X_traj_full[PERTURBATION_TIME:]
        clip_post = clip_count_full[PERTURBATION_TIME:]
        total_clips_post = int(clip_post.sum().item())

        world_dict = meta.copy()
        world_dict['world_id'] = wid
        analysis = analyze_trajectory(X_post, clip_post, total_clips_post, world_dict)
        post_regime = analysis['primary_regime']
        regime_changed = post_regime != original_regime
        print(f"    regime: {original_regime} → {post_regime}"
              f"{' (CHANGED)' if regime_changed else ''}")

        # ── Recovery metrics ──
        eq_converged = analysis['converged']
        eq_time = analysis['convergence_time']
        recovery_time = None
        recovery_failure = False
        if eq_converged:
            recovery_time = PERTURBATION_TIME + eq_time  # absolute time (aligned with run_002)
        else:
            recovery_failure = True

        # max deviation from baseline (pre-perturbation baseline = 50-step window)
        pre_baseline = X_traj_full[max(0, PERTURBATION_TIME - 50):PERTURBATION_TIME]
        pre_mean = pre_baseline.mean(dim=0)
        post_devs = (X_traj_full[PERTURBATION_TIME:] - pre_mean).abs()
        max_deviation = float(post_devs.max().item())

        # Save perturbed trajectory
        torch.save({
            'X_traj': X_traj_full,
            'clip_count': clip_count_full,
            'total_clips': total_clips_post,
            'world_seed': orig_traj['world_seed'],
            'cell_seed': cell_seed,
            'world_id': wid,
            'original_regime': original_regime,
            'perturbed_regime': post_regime,
            'knockdown_gene': kd_gene,
            'knockdown_gene_role': kd_gene_role,
            'knockdown_gene_expr_at_intervention': max_expr,
            'intervention_time': PERTURBATION_TIME,
            'perturbation_type': 'knockdown',
        }, os.path.join(PERTURBATIONS_DIR, f'{wid}_{cell_label}_perturb.pt'))

        summary.append({
            'world_id': wid,
            'cell': cell_label,
            'topology_type': topo,
            'delta': cand['delta'],
            'strength_regime': cand['strength_regime'],
            'sign_ratio': cand['sign_ratio'],
            'combo_key': cand['combo_key'],
            'primary_regime': original_regime,
            'perturbed_regime': post_regime,
            'perturbed_secondary_labels': analysis['secondary_labels'],
            'regime_changed': regime_changed,
            'knockdown_gene': kd_gene,
            'knockdown_gene_role': kd_gene_role,
            'knockdown_gene_expr_at_intervention': max_expr,
            'intervention_time': PERTURBATION_TIME,
            'recovery_time': recovery_time,
            'recovery_failure': recovery_failure,
            'max_deviation_from_baseline': max_deviation,
            'equilibrium_converged': eq_converged,
            'spectral_radius': analysis['spectral_radius'],
            'max_amplitude': analysis['max_amplitude'],
            'final_amplitude': analysis['final_amplitude'],
            'clipping_frequency': analysis['clipping_frequency'],
            'clipping_dominated': analysis['clipping_dominated'],
            'is_sustained': analysis['is_sustained'],
        })

    # Save summary
    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    with open(os.path.join(ANALYSIS_DIR, 'perturbation_summary.json'), 'w') as f:
        json.dump(summary, f, indent=2)
    changed = sum(1 for s in summary if s['regime_changed'])
    print(f"\n  perturbation done: {changed}/{len(summary)} regimes changed")


def plot_perturbation_comparison() -> None:
    """Before/after perturbation trajectory plots, grouped by primary_regime.

    Per Output Spec §10: show original + perturbed side by side for each candidate,
    faceted by original primary regime.
    """
    pert_file = os.path.join(ANALYSIS_DIR, 'perturbation_summary.json')
    if not os.path.exists(pert_file):
        print("  skip perturbation plots: no summary found")
        return
    with open(pert_file) as f:
        pert_summary = json.load(f)

    if not pert_summary:
        return

    pert_dir = os.path.join(FIGURES_DIR, 'perturbation')
    os.makedirs(pert_dir, exist_ok=True)

    # Group by primary_regime
    grouped: Dict[str, List[Dict]] = defaultdict(list)
    for s in pert_summary:
        grouped[s['primary_regime']].append(s)

    for regime, entries in sorted(grouped.items()):
        n = len(entries)
        fig, axes = plt.subplots(n, 2, figsize=(14, 4 * n))
        if n == 1:
            axes = np.array([axes])

        for r, entry in enumerate(entries):
            wid = entry['world_id']
            cell_label = entry['cell']
            topo = entry['topology_type']

            # Original trajectory
            tp_o = os.path.join(TRAJECTORIES_DIR, f'{wid}_{cell_label}.pt')
            if os.path.exists(tp_o):
                orig = torch.load(tp_o, map_location='cpu', weights_only=False)
                Xo = orig['X_traj'].numpy()
                for g in range(G):
                    if topo == 'tf_restricted' and g in TF_GENES:
                        axes[r, 0].plot(Xo[:, g], color='red', alpha=0.8, linewidth=1.0)
                    elif topo == 'tf_restricted':
                        axes[r, 0].plot(Xo[:, g], color='grey', alpha=0.35, linewidth=0.8)
                    else:
                        axes[r, 0].plot(Xo[:, g], alpha=0.9, linewidth=0.8)
                axes[r, 0].axvline(x=PERTURBATION_TIME, color='gray', linestyle='--',
                                   linewidth=0.8, alpha=0.6)
                axes[r, 0].set_title(f'{wid}\n(original)', fontsize=10, linespacing=1.3)
                axes[r, 0].set_ylabel('gene activity')
                axes[r, 0].tick_params(labelsize=8)

            # Perturbed trajectory
            tp_p = os.path.join(PERTURBATIONS_DIR, f'{wid}_{cell_label}_perturb.pt')
            if os.path.exists(tp_p):
                pert = torch.load(tp_p, map_location='cpu', weights_only=False)
                Xp = pert['X_traj'].numpy()
                post_regime = pert.get('perturbed_regime', '?')
                for g in range(G):
                    if topo == 'tf_restricted' and g in TF_GENES:
                        axes[r, 1].plot(Xp[:, g], color='red', alpha=0.8, linewidth=1.0)
                    elif topo == 'tf_restricted':
                        axes[r, 1].plot(Xp[:, g], color='grey', alpha=0.35, linewidth=0.8)
                    else:
                        axes[r, 1].plot(Xp[:, g], alpha=0.9, linewidth=0.8)
                axes[r, 1].axvline(x=PERTURBATION_TIME, color='red', linestyle='--',
                                   linewidth=0.8, alpha=0.6)
                axes[r, 1].set_title(
                    f'{wid} [{post_regime}]\n(KD gene{pert.get("knockdown_gene","?")})',
                    fontsize=10, linespacing=1.3)
                axes[r, 1].set_ylabel('gene activity')
                axes[r, 1].tick_params(labelsize=8)

            for ax in axes[r]:
                ax.set_xlabel('t')
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

        safe_regime = regime.lower().replace(' ', '_')
        fig.suptitle(f'Perturbation — {regime}', fontsize=14, fontweight='bold', y=0.99)
        plt.tight_layout(rect=[0, 0, 1, 0.98])
        fig.savefig(os.path.join(pert_dir, f'perturbation_{safe_regime}.png'),
                    dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"  perturbation figure: perturbation_{safe_regime}.png ({n} entries)")


# ═══════════════════════════════════════════════════════════════════
# main
# ═══════════════════════════════════════════════════════════════════

def _load_existing_results() -> List[Dict[str, Any]]:
    """Load existing analysis JSONs and reconstruct all_results list.

    Parses world_id to extract combo metadata (topology_type, delta, strength_regime, sign_ratio).
    Groups per-cell analysis JSONs by world_id.
    """
    import re
    from collections import defaultdict

    # Collect all analysis JSONs, grouped by world_id
    world_cells: Dict[str, List[Dict]] = defaultdict(list)
    for fname in os.listdir(ANALYSIS_DIR):
        if not fname.endswith('.json'):
            continue
        if fname == 'perturbation_summary.json':
            continue
        with open(os.path.join(ANALYSIS_DIR, fname)) as f:
            data = json.load(f)
        world_cells[data['world_id']].append(data)

    # Parse world_id: e.g. unrestricted_sparse_d0p1_chen_stress_balanced_t009
    # pattern: {topo}_{delta_label}_{strength_regime}_{sign_ratio}_t{ti}
    _parse_pattern = re.compile(
        r'^(?P<topology_type>unrestricted_sparse|tf_restricted)'
        r'_(?P<delta_label>d0p[124])'
        r'_(?P<strength_regime>baseline|chen_moderate|stress|chen_stress)'
        r'_(?P<sign_ratio>balanced|repression_biased)'
        r'_t(?P<ti>\d{3})$'
    )

    delta_map = {'d0p1': 0.1, 'd0p2': 0.2, 'd0p4': 0.4}

    all_results: List[Dict[str, Any]] = []
    for world_id, cells in sorted(world_cells.items()):
        m = _parse_pattern.match(world_id)
        if not m:
            continue
        topo = m.group('topology_type')
        dlabel = m.group('delta_label')
        dval = delta_map[dlabel]
        sr = m.group('strength_regime')
        sign = m.group('sign_ratio')
        ti = int(m.group('ti'))

        # Sort cells by cell_idx
        cells.sort(key=lambda c: c['cell_idx'])

        # Build per_cell list (strip world_id top-level key,
        # keep cell_idx needed by run_per_trajectory_analysis)
        per_cell = []
        for c in cells:
            entry = {k: v for k, v in c.items() if k not in ('world_id',)}
            per_cell.append(entry)

        result_stub = {
            'world_id':        world_id,
            'combo_key':       f"{topo}_{dlabel}_{sr}_{sign}",
            'topology_type':   topo,
            'delta':           dval,
            'strength_regime': sr,
            'sign_ratio':      sign,
            'topo_idx':        ti,
            'per_cell':        per_cell,
        }
        all_results.append(result_stub)

    print(f"  Loaded {len(all_results)} worlds ({sum(len(c) for c in world_cells.values())} cells)")
    return all_results


def main():
    # ── Step 0:  Ensure directories exist ──
    for d in [SCRIPTS_DIR, FIGURES_DIR, TABLES_DIR, TRAJECTORIES_DIR,
              METADATA_DIR, ANALYSIS_DIR, SUMMARY_DIR, PERTURBATIONS_DIR]:
        os.makedirs(d, exist_ok=True)

    # ── Step 1:  Generate and validate Full Grid config ──
    print("=" * 60)
    print(f"  {RUN_ID} — Model A Decay / Memory Boundary")
    print(f"  execution_grid = {EXECUTION_GRID}")
    print("=" * 60)

    all_combos = generate_all_combos()
    print_grid_summary(all_combos)
    assert validate_seed_ranges(all_combos), "Seed overlap detected!"
    print(f"\n── Example world IDs ──")
    for combo in all_combos[:3]:
        for ti in range(2):
            wid = make_world_id(combo['topology_type'], combo['delta_label'],
                                combo['strength_regime'], combo['sign_ratio'], ti)
            print(f"  {wid}")
    print(f"\n── Step 1 complete.  {len(all_combos)} combos configured. ──")

    # ── Step 2:  World generation & trajectory simulation ──
    all_results = simulate_all_combos(all_combos)
    assign_hub_likeness()   # post-process: compute quartiles & is_hub_like

    # ── Step 3:  Per-trajectory analysis ──
    all_results = run_per_trajectory_analysis(all_results)

    # ── Step 4:  Aggregation ──
    aggregation = aggregate_across_combos(all_results)
    cross_dim = analyze_cross_dimension(all_results)

    # ── Step 5:  Tables ──
    write_tables(all_results, aggregation)

    # ── Step 6:  Figures ──
    generate_all_figures(all_results)

    # ── Step 7:  Perturbation ──
    run_perturbation_analysis(all_results)
    plot_perturbation_comparison()

    print("\nDone.")


if __name__ == '__main__':
    main()
