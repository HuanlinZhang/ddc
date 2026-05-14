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
        raise NotImplementedError

    def from_dict(self, data: Dict[str, Any]) -> None:
        raise NotImplementedError


def sample_world(seed: int) -> World:
    raise NotImplementedError


def sample_initial_state(cell_seed: int, world: World) -> Tuple[Tensor, Tensor]:
    raise NotImplementedError


def compute_regulatory_response(P: Tensor, world: World) -> Tensor:
    raise NotImplementedError


def update_mRNA(X: Tensor, P: Tensor, world: World) -> Tensor:
    raise NotImplementedError


def update_protein(P: Tensor, X: Tensor, world: World) -> Tensor:
    raise NotImplementedError


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
