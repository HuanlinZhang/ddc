"""
DDC Lite Level 0 - Minimal Learnable DDC
=========================================

A minimal synthetic GRN system for AI learnability benchmarking.
Single-input signed-regulation: one TF per gene, activation or repression.

Based on: docs/01_DDC_Lite_Curriculum/Level_0_Minimal_TF/*

Author: zhanghl
Version: v0.1 (commit 1: module skeleton)

================================================================================
Differences from Phase 0 (src/ddc.py)
================================================================================
- G = 20 (vs 50)
- 1 gene <- 1 TF (vs 1~3 TF per gene)
- Signed single-input Hill response (vs multiplicative combinatorial TFinput)
- No tilde_P normalization
- No chromatin gating (Z removed)
- No resource projection
- No fate (N removed)
- State = {X, P} only
- Edge signs: activation (+1) or repression (-1)
"""
import json
import os
import torch
import copy
from typing import Dict, Tuple, List, Any


G: int = 20

N_TF: int = 4

T: int = 200

DTYPE: torch.dtype = torch.float64

ACTIVATION: int = 1
REPRESSION: int = -1

Tensor = torch.Tensor
State = Dict[str, Any]

TF_GENES: List[int] = list(range(0, N_TF))


class World:
    def __init__(self, seed: int):
        self.seed: int = seed
        self.P_graph: Dict[int, List[int]] = {}
        self.edge_signs: Dict[int, Dict[int, int]] = {}
        self.rho: Tensor = torch.zeros(G, dtype=DTYPE)
        self.K: Tensor = torch.zeros(G, dtype=DTYPE)
        self.n: Tensor = torch.zeros(G, dtype=DTYPE)
        self.delta_x: Tensor = torch.zeros(G, dtype=DTYPE)
        self.delta_p: Tensor = torch.zeros(G, dtype=DTYPE)
        self.gamma: Tensor = torch.zeros(G, dtype=DTYPE)

    def to_dict(self) -> Dict[str, Any]:
        return {
            'seed': self.seed,
            'P_graph': {str(k): v for k, v in self.P_graph.items()},
            'edge_signs': {
                str(k): {str(sk): sv for sk, sv in v.items()}
                for k, v in self.edge_signs.items()
            },
            'parameters': {
                'rho': self.rho.tolist(),
                'K': self.K.tolist(),
                'n': self.n.tolist(),
                'delta_x': self.delta_x.tolist(),
                'delta_p': self.delta_p.tolist(),
                'gamma': self.gamma.tolist(),
            },
        }

    def from_dict(self, data: Dict[str, Any]) -> None:
        self.seed = data['seed']
        self.P_graph = {int(k): v for k, v in data['P_graph'].items()}
        self.edge_signs = {
            int(k): {int(sk): sv for sk, sv in v.items()}
            for k, v in data['edge_signs'].items()
        }
        params = data['parameters']
        self.rho = torch.tensor(params['rho'], dtype=DTYPE)
        self.K = torch.tensor(params['K'], dtype=DTYPE)
        self.n = torch.tensor(params['n'], dtype=DTYPE)
        self.delta_x = torch.tensor(params['delta_x'], dtype=DTYPE)
        self.delta_p = torch.tensor(params['delta_p'], dtype=DTYPE)
        self.gamma = torch.tensor(params['gamma'], dtype=DTYPE)


def sample_world(seed: int) -> World:
    rng: torch.Generator = torch.Generator()
    rng.manual_seed(seed)

    world: World = World(seed)

    world.rho = torch.empty(G, dtype=DTYPE).uniform_(0.5, 2.0, generator=rng)
    world.K = torch.empty(G, dtype=DTYPE).uniform_(1.0, 5.0, generator=rng)
    world.n = torch.full((G,), 2.0, dtype=DTYPE)
    world.delta_x = torch.empty(G, dtype=DTYPE).uniform_(0.1, 0.5, generator=rng)
    world.delta_p = torch.empty(G, dtype=DTYPE).uniform_(0.05, 0.3, generator=rng)
    world.gamma = torch.full((G,), 1.0, dtype=DTYPE)

    for i in range(G):
        tf_pool: List[int] = [tf for tf in TF_GENES if tf != i]
        idx: int = int(torch.randint(0, len(tf_pool), (1,), generator=rng).item())
        regulator: int = tf_pool[idx]
        world.P_graph[i] = [regulator]

        world.edge_signs[i] = {}
        sign: int = ACTIVATION if torch.rand(1, generator=rng).item() < 0.5 else REPRESSION
        world.edge_signs[i][regulator] = sign

    return world


def sample_initial_state(cell_seed: int, world: World) -> Tuple[Tensor, Tensor]:
    raise NotImplementedError


def compute_regulatory_response(P: Tensor, world: World) -> Tensor:
    R: Tensor = torch.zeros(G, dtype=DTYPE)
    for i in range(G):
        j: int = world.P_graph[i][0]
        s: int = world.edge_signs[i][j]
        pj_n: Tensor = P[j] ** world.n[i]
        Ki_n: Tensor = world.K[i] ** world.n[i]
        denom: Tensor = Ki_n + pj_n
        if s == ACTIVATION:
            R[i] = pj_n / denom
        else:
            R[i] = Ki_n / denom
    return R


def update_mRNA(X: Tensor, P: Tensor, world: World) -> Tensor:
    R: Tensor = compute_regulatory_response(P, world)
    X_next: Tensor = (1.0 - world.delta_x) * X + world.rho * R
    return X_next


def update_protein(P: Tensor, X: Tensor, world: World) -> Tensor:
    P_next: Tensor = (1.0 - world.delta_p) * P + world.gamma * X
    return P_next


def simulate_single_cell(
    world: World,
    X0: Tensor,
    P0: Tensor,
    t_steps: int = T,
    intervention_time: int = None,
    intervention_config: Dict = None,
) -> Dict[str, Tensor]:
    raise NotImplementedError


def run_simulation(
    seed: int,
    save_path: str = None,
    intervention_time: int = None,
    intervention_config: Dict = None,
) -> Dict[str, Any]:
    raise NotImplementedError


def generate_dataset(
    world_seed: int,
    M: int,
    save_path: str = None,
) -> Tuple[Tensor, World, List[int]]:
    raise NotImplementedError


def apply_perturbation(
    world: World,
    state: State,
    config: Dict[str, Any],
) -> Tuple[World, State]:
    raise NotImplementedError


def apply_intervention(
    state: Tuple[Tensor, Tensor],
    config: Dict[str, Any],
) -> Tuple[Tensor, Tensor]:
    raise NotImplementedError
