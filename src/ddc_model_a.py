"""
DDC Model A - Classical Linear GRN (Chen-style sparse)
========================================================

Single-layer linear Gene Regulatory Network with Chen-style sparse topology.
No protein layer. No nonlinear regulation.

Equation:
    raw_x_i(t+1) = (1 - delta_i) * x_i(t) + b_i + sum_j a_ij * x_j(t)
    x_i(t+1)     = max(0, raw_x_i(t+1))

Matrix form (linear core, before clipping):
    X(t+1) = A @ X(t) + b

where:
    A[i][i] = 1 - delta_i      (diagonal decay)
    A[i][j] = sign_ij * strength_ij  (signed linear regulatory strength)
    A[i][k] = 0                 (all other off-diagonal entries)

Topology: Chen-style sparse random graph
    G = 20, r = 0.1
    380 candidate edges (j -> i, j != i), 38 sampled without replacement
    Each gene may have 0, 1, or multiple regulators

Design reference:
    docs/01_DDC_Lite_Curriculum/run_002/   

Author: zhanghl
Version: v0.2
"""
import json
import os
import torch
from typing import Dict, List, Tuple, Any


G: int = 20

EDGE_DENSITY_R: float = 0.1

T_SIM: int = 1000

DTYPE: torch.dtype = torch.float64

ACTIVATION: int = 1
REPRESSION: int = -1

SIGN_RATIO: float = 0.5

DEFAULT_DELTA: float = 0.2
DEFAULT_B: float = 0.1

A_MIN_BASELINE: float = 0.02
A_MAX_BASELINE: float = 0.10

A_MIN_STRESS: float = 0.10
A_MAX_STRESS: float = 0.30

A_MIN_CHEN_MODERATE: float = 0.10
A_MAX_CHEN_MODERATE: float = 0.20

A_MIN_CHEN_STRESS: float = 0.20
A_MAX_CHEN_STRESS: float = 0.40

Tensor = torch.Tensor

N_EDGES: int = int(G * (G - 1) * EDGE_DENSITY_R)


class World:
    def __init__(self, seed: int):
        self.seed: int = seed
        self.P_graph: Dict[int, List[int]] = {}
        self.edge_signs: Dict[int, Dict[int, int]] = {}
        self.edge_strengths: Dict[int, Dict[int, float]] = {}
        self.delta: Tensor = torch.full((G,), DEFAULT_DELTA, dtype=DTYPE)
        self.b: Tensor = torch.full((G,), DEFAULT_B, dtype=DTYPE)
        self._A: Tensor = None

    def build_A(self) -> Tensor:
        if self._A is not None:
            return self._A
        A = torch.diag(1.0 - self.delta)
        for i in range(G):
            for j in self.P_graph.get(i, []):
                sign = self.edge_signs[i][j]
                strength = self.edge_strengths[i][j]
                A[i, j] = sign * strength
        self._A = A
        return A

    def in_degrees(self) -> Dict[int, int]:
        return {i: len(self.P_graph.get(i, [])) for i in range(G)}

    def out_degrees(self) -> Dict[int, int]:
        od = {i: 0 for i in range(G)}
        for i, regs in self.P_graph.items():
            for j in regs:
                od[j] = od.get(j, 0) + 1
        return od

    def isolated_genes(self) -> List[int]:
        return [i for i in range(G) if len(self.P_graph.get(i, [])) == 0]

    def to_dict(self) -> Dict[str, Any]:
        return {
            'seed': self.seed,
            'P_graph': {str(k): v for k, v in self.P_graph.items()},
            'edge_signs': {
                str(k): {str(sk): sv for sk, sv in v.items()}
                for k, v in self.edge_signs.items()
            },
            'edge_strengths': {
                str(k): {str(sk): sv for sk, sv in v.items()}
                for k, v in self.edge_strengths.items()
            },
            'parameters': {
                'delta': self.delta.tolist(),
                'b': self.b.tolist(),
            },
            'in_degrees': {str(k): v for k, v in self.in_degrees().items()},
            'out_degrees': {str(k): v for k, v in self.out_degrees().items()},
            'isolated_genes': self.isolated_genes(),
        }

    def from_dict(self, data: Dict[str, Any]) -> None:
        self.seed = data['seed']
        self.P_graph = {int(k): [int(r) for r in v] for k, v in data['P_graph'].items()}
        self.edge_signs = {
            int(k): {int(sk): sv for sk, sv in v.items()}
            for k, v in data['edge_signs'].items()
        }
        self.edge_strengths = {
            int(k): {int(sk): sv for sk, sv in v.items()}
            for k, v in data.get('edge_strengths', {}).items()
        }
        params = data['parameters']
        self.delta = torch.tensor(params['delta'], dtype=DTYPE)
        self.b = torch.tensor(params['b'], dtype=DTYPE)
        self._A = None


def _sample_edge_set(seed: int, sign_ratio: float) -> Tuple[Dict[int, List[int]], Dict[int, Dict[int, int]]]:
    rng = torch.Generator()
    rng.manual_seed(seed)

    candidates = [(j, i) for i in range(G) for j in range(G) if j != i]
    indices = torch.randperm(len(candidates), generator=rng)[:N_EDGES].tolist()

    n_act = round(N_EDGES * sign_ratio)
    sign_pool = [ACTIVATION] * n_act + [REPRESSION] * (N_EDGES - n_act)
    sign_order = torch.randperm(N_EDGES, generator=rng).tolist()

    P_graph: Dict[int, List[int]] = {}
    edge_signs: Dict[int, Dict[int, int]] = {}

    for k in range(N_EDGES):
        idx = indices[k]
        j, i = candidates[idx]
        if i not in P_graph:
            P_graph[i] = []
        P_graph[i].append(j)
        if i not in edge_signs:
            edge_signs[i] = {}
        edge_signs[i][j] = sign_pool[sign_order[k]]

    for i in range(G):
        if i not in P_graph:
            P_graph[i] = []

    return P_graph, edge_signs


def sample_topology(seed: int, sign_ratio: float = SIGN_RATIO) -> World:
    world = World(seed)
    world.P_graph, world.edge_signs = _sample_edge_set(seed, sign_ratio)
    return world


def sample_world(seed: int, world: World = None,
                 a_min: float = A_MIN_BASELINE, a_max: float = A_MAX_BASELINE) -> World:
    rng = torch.Generator()
    rng.manual_seed(seed)

    if world is None:
        w = sample_topology(seed)
    else:
        # 不能直接 w = world：World 是可变对象，直接赋引用会导致后续
        # 对 w.edge_strengths / w.delta / w.b 的修改污染传入的 world。
        # 这里创建新对象 + 深拷贝拓扑，保证传入 world 始终只读。
        w = World(seed)
        w.P_graph = {k: v[:] for k, v in world.P_graph.items()}
        w.edge_signs = {
            k: {sk: sv for sk, sv in v.items()}
            for k, v in world.edge_signs.items()
        }

    w.edge_strengths = {}
    for i in range(G):
        if i in w.P_graph and len(w.P_graph[i]) > 0:
            w.edge_strengths[i] = {}
            for j in w.P_graph[i]:
                w.edge_strengths[i][j] = float(
                    torch.empty(1, dtype=DTYPE).uniform_(a_min, a_max, generator=rng).item()
                )
        else:
            w.edge_strengths[i] = {}

    w.delta = torch.full((G,), DEFAULT_DELTA, dtype=DTYPE)
    w.b = torch.full((G,), DEFAULT_B, dtype=DTYPE)
    w._A = None

    return w


def sample_initial_state(cell_seed: int) -> Tensor:
    rng = torch.Generator()
    rng.manual_seed(cell_seed)
    X0 = torch.empty(G, dtype=DTYPE).uniform_(0, 1, generator=rng)
    return X0


def apply_intervention(X: Tensor, config: Dict[str, Any]) -> Tensor:
    X = X.clone()
    if 'knockdown' in config:
        for gene in config['knockdown']:
            X[gene] = 0.0
    if 'scale' in config:
        for gene, scale in config['scale']:
            X[gene] *= scale
    if 'set' in config:
        for gene, val in config['set']:
            X[gene] = float(val)
    return X


def simulate_single_cell(
    world: World,
    X0: Tensor,
    t_steps: int = T_SIM,
    intervention_time: int = None,
    intervention_config: Dict[str, Any] = None,
) -> Dict[str, Any]:
    A = world.build_A()
    b = world.b
    X = X0.clone()
    X_traj = torch.zeros((t_steps + 1, G), dtype=DTYPE)
    clip_count = torch.zeros(t_steps + 1, dtype=torch.int64)
    intervention_history = None
    X_traj[0] = X
    clip_count[0] = 0

    for t in range(t_steps):
        if intervention_time is not None and t == intervention_time:
            X = apply_intervention(X, intervention_config)
            X_traj[t] = X
            clip_count[t] = 0
            intervention_history = {
                'time': intervention_time,
                'knockdown_genes': intervention_config.get('knockdown'),
                'scale_genes': intervention_config.get('scale'),
                'set_genes': intervention_config.get('set'),
            }
            continue

        raw = A @ X + b
        X_next = torch.clamp(raw, min=0.0)
        X = X_next
        X_traj[t + 1] = X
        clip_count[t + 1] = (raw != X).sum().item()

    return {
        'X_traj': X_traj,
        'clip_count': clip_count,
        'total_clips': int(clip_count.sum().item()),
        'intervention_history': intervention_history,
        'world_seed': world.seed,
        'world': world.to_dict(),
    }


def run_simulation(
    seed: int,
    world: World = None,
    save_path: str = None,
    a_min: float = None, a_max: float = None,
    intervention_time: int = None,
    intervention_config: Dict[str, Any] = None,
) -> Dict[str, Any]:
    world_seed = seed
    cell_seed = seed + 1

    _amin = a_min if a_min is not None else A_MIN_BASELINE
    _amax = a_max if a_max is not None else A_MAX_BASELINE

    if world is None:
        w = sample_world(world_seed, a_min=_amin, a_max=_amax)
    else:
        w = sample_world(world_seed, world=world, a_min=_amin, a_max=_amax)

    X0 = sample_initial_state(cell_seed)
    traj = simulate_single_cell(w, X0, T_SIM,
        intervention_time=intervention_time,
        intervention_config=intervention_config)
    traj['cell_seed'] = cell_seed

    if save_path is not None:
        torch.save({
            'X_traj': traj['X_traj'],
            'clip_count': traj['clip_count'],
            'total_clips': traj['total_clips'],
            'world_seed': traj['world_seed'],
            'cell_seed': traj['cell_seed'],
            'world': traj['world'],
        }, save_path)

    return traj


def compute_stability_metrics(
    traj: Dict[str, Any],
    eps: float = 1e-4,
    window: int = 50,
    collapse_threshold: float = 1e-3,
    divergence_threshold: float = 1e3,
) -> Dict[str, Any]:
    X_traj = traj['X_traj']
    max_expression = float(X_traj.max().item())
    final_mean_expression = float(X_traj[-1].mean().item())
    bounded = max_expression < divergence_threshold
    collapsed = final_mean_expression < collapse_threshold
    clipping_dominated = traj.get('total_clips', 0) > X_traj.shape[0] * G * 0.1
    converged = False
    convergence_time = -1
    t_steps = X_traj.shape[0] - 1

    for t in range(t_steps - window + 1):
        consecutive = True
        for w in range(window):
            step_diff = float(torch.norm(X_traj[t + w + 1] - X_traj[t + w]).item())
            if step_diff >= eps:
                consecutive = False
                break
        if consecutive:
            converged = True
            convergence_time = t
            break

    return {
        'converged': converged,
        'convergence_time': convergence_time,
        'bounded': bounded,
        'collapsed': collapsed,
        'clipping_dominated': clipping_dominated,
        'total_clips': traj.get('total_clips', 0),
        'max_expression': max_expression,
        'final_mean_expression': final_mean_expression,
    }


def aggregate_stability_stats(metrics_list: List[Dict[str, Any]]) -> Dict[str, Any]:
    total = len(metrics_list)
    if total == 0:
        return {'total_worlds': 0}
    n_converged = sum(1 for m in metrics_list if m['converged'])
    n_bounded = sum(1 for m in metrics_list if m['bounded'])
    n_collapsed = sum(1 for m in metrics_list if m['collapsed'])
    n_clipping_dominated = sum(1 for m in metrics_list if m.get('clipping_dominated', False))
    return {
        'total_worlds': total,
        'converged': n_converged,
        'bounded': n_bounded,
        'collapsed': n_collapsed,
        'clipping_dominated': n_clipping_dominated,
        'convergence_rate': float(n_converged) / total,
        'bounded_rate': float(n_bounded) / total,
        'collapse_rate': float(n_collapsed) / total,
        'clipping_dominated_rate': float(n_clipping_dominated) / total,
    }
