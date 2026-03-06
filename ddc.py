"""
DDC Phase 0 - Designed Digital-Cell Model
==========================================

A computational simulation framework for gene regulatory network dynamics.
Based on the Phase 0 specification documents.

Author: zhanghl
Version: v1.1
Status: Active

================================================================================
Version History & Update Notes
================================================================================

v1.1 (2026-03-06):
    - Modified normalize_protein() function: changed protein normalization
      from sum-based to mean-based calculation
    - Previous:  tilde_P = P / (torch.sum(P) + epsilon)
    - Current:   tilde_P = P / (torch.mean(P) + epsilon)
    - Reason:    To address data decay issue where X/P values rapidly approach
      zero. The mean-based normalization compensates for the signal dilution
      effect caused by network scale (G=50 genes), allowing TF input to
      cross Hill function activation threshold.

v1.0 (2026-03-03):
    - Initial release
    - Status: Frozen
"""
import json
import os
import torch
import copy
from typing import Dict, Tuple, List, Any

# ==========================================
# Global Constants
# ==========================================
G: int = 50
T: int = 200
R_TOTAL: float = 1.0
# ### 0305 方案2：增大R_total，其余不变
# R_TOTAL: float = 25.0
# ### 经测试，无效果
EPSILON: float = 1e-8
K_POP: float = 1.0
DTYPE: torch.dtype = torch.float64

Tensor = torch.Tensor
State = Dict[str, Any]

# Gene Categories Indices
TF_GENES: List[int] = list(range(0, 6))
EPI_GENES: List[int] = list(range(17, 20))

### 0304 ADD: Other Gene Categories Indices
RBP_GENES: List[int] = list(range(6, 11))
KINASE_GENES: List[int] = list(range(11, 14))
PHOSPHATASE_GENES: List[int] = list(range(14, 17))

CELL_CYCLE_GENES: List[int] = list(range(20, 23))
APOPTOSIS_GENES: List[int] = list(range(23, 26))

BACKGROUND_GENES: List[int] = list(range(26, 50))

GENE_CATEGORIES: Dict[str, Dict[str, list[int]]] = {
    "Core Regulatory": {
        "TF": TF_GENES,
        "RBP": RBP_GENES,
        "Kinase": KINASE_GENES,
        "Phosphatase": PHOSPHATASE_GENES,
        "Epigenetic": EPI_GENES
    },
    "Cell Fate": {
        "Cell Cycle": CELL_CYCLE_GENES,
        "Apoptosis": APOPTOSIS_GENES
    },
    "Background": {
        "Background": BACKGROUND_GENES
    }
}

# reverse mapping dictionary for plotting
GENE_TO_MACRO: Dict[int, str] = {}
GENE_TO_MICRO: Dict[int, str] = {}
for macro, micro_dict in GENE_CATEGORIES.items():
    for micro, gene_list in micro_dict.items():
        for g in gene_list:
            GENE_TO_MACRO[g] = macro
            GENE_TO_MICRO[g] = micro
### 0304 ADD END

def stable_sigmoid(x: Tensor) -> Tensor:
    return torch.where(x >= 0,
                       1.0 / (1.0 + torch.exp(-x)),
                       torch.exp(x) / (1.0 + torch.exp(x)))

# ==========================================
# World Object Definition
# ==========================================
class World():
    def __init__(self, seed: int):
        self.seed: int = seed

        # Regulation Graphs
        self.P_graph: Dict[int, List[int]] = {}
        self.E_graph: Dict[int, List[int]] = {}

        # Gene-level Parameters
        self.alpha: Tensor = torch.zeros(G, dtype=DTYPE)
        self.rho: Tensor = torch.zeros(G, dtype=DTYPE)
        self.K: Tensor = torch.zeros(G, dtype=DTYPE)
        self.n: Tensor = torch.zeros(G, dtype=DTYPE)
        self.delta_x: Tensor = torch.zeros(G, dtype=DTYPE)
        self.delta_p: Tensor = torch.zeros(G, dtype=DTYPE)
        self.gamma: Tensor = torch.zeros(G, dtype=DTYPE)

        # Edge-level Parameters
        self.a_ij: Dict[int, Dict[int, float]] = {}
        self.beta_ij: Dict[int, Dict[int, float]] = {}

        # Run-level Parameters
        self.r: float = 0.0
        self.K_pop: float = K_POP

        # Global constants
        self.R_total: float = R_TOTAL
        self.epsilon: float = EPSILON

        ### 0304 ADD
        self.gene_categories: Dict[str, Dict[str, list[int]]] = GENE_CATEGORIES
        self.gene_to_macro: Dict[int, str] = GENE_TO_MACRO
        self.gene_to_micro: Dict[int, str] = GENE_TO_MICRO
        ###
    
    def to_dict(self) -> Dict[str, Any]:
        # Convert tensors to lists for JSON serialization
        return {
            'seed': self.seed,
            'P_graph': self.P_graph,
            'E_graph': self.E_graph,
            'parameters': {
                'alpha': self.alpha.tolist(),
                'rho': self.rho.tolist(),
                'K': self.K.tolist(),
                'n': self.n.tolist(),
                'delta_x': self.delta_x.tolist(),
                'delta_p': self.delta_p.tolist(),
                'gamma': self.gamma.tolist(),
                'a_ij': self.a_ij,
                'beta_ij': self.beta_ij,
                'r': self.r,
                'K_pop': self.K_pop,
                'R_total': self.R_total,
                'epsilon': self.epsilon,
            },
            ### 0304 ADD
            'gene_annotation': {
                'category': self.gene_categories,
                'to_macro': self.gene_to_macro,
                'to_micro': self.gene_to_micro
            }
            ###
        }

    def from_dict(self, data: Dict[str, Any]) -> None:
        self.seed = data['seed']
        
        self.P_graph = {int(k): v for k, v in data['P_graph'].items()}
        self.E_graph = {int(k): v for k, v in data['E_graph'].items()}

        params = data['parameters']
        self.alpha = torch.tensor(params['alpha'], dtype=DTYPE)
        self.rho = torch.tensor(params['rho'], dtype=DTYPE)
        self.K = torch.tensor(params['K'], dtype=DTYPE)
        self.n = torch.tensor(params['n'], dtype=DTYPE)
        self.delta_x = torch.tensor(params['delta_x'], dtype=DTYPE)
        self.delta_p = torch.tensor(params['delta_p'], dtype=DTYPE)
        self.gamma = torch.tensor(params['gamma'], dtype=DTYPE)

        self.a_ij = {int(k): {int(sub_k): float(sub_v) for sub_k, sub_v in v.items()}
                     for k, v in params['a_ij'].items()}
        self.beta_ij = {int(k): {int(sub_k): float(sub_v) for sub_k, sub_v in v.items()}
                        for k, v in params['beta_ij'].items()}
        
        self.r = params['r']
        self.K_pop = params['K_pop']
        self.R_total = params['R_total']
        self.epsilon = params['epsilon']

        ### 0304 ADD
        self.gene_categories = data['gene_annotation']['category']
        self.gene_to_macro = {int(k): v for k, v in data['gene_annotation']['to_macro'].items()}
        self.gene_to_micro = {int(k): v for k, v in data['gene_annotation']['to_micro'].items()}
        ###

# ==========================================
# Monte Carlo World Sampler
# ==========================================
def sample_world(seed: int) -> World:
    rng: torch.Generator = torch.Generator()
    rng.manual_seed(seed)
    
    world: World = World(seed)

    world.alpha = torch.empty(G, dtype=DTYPE).normal_(0, 1, generator=rng)
    world.rho = torch.empty(G, dtype=DTYPE).uniform_(0.5, 2.0, generator=rng)
    world.K = torch.empty(G, dtype=DTYPE).uniform_(0.1, 1.0, generator=rng)
    # ### 0305 adjust K to prevent X, P from decreasing too fast 
    # world.K = torch.empty(G, dtype=DTYPE).uniform_(0.01, 0.1, generator=rng)
    # ###
    world.n = torch.full((G,), 2.0, dtype=DTYPE)
    world.delta_x = torch.empty(G, dtype=DTYPE).uniform_(0.1, 0.5, generator=rng)
    world.delta_p = torch.empty(G, dtype=DTYPE).uniform_(0.05, 0.3, generator=rng)
    world.gamma = torch.full((G,), 1.0, dtype=DTYPE)
    world.r = float(torch.empty(1, dtype=DTYPE).uniform_(0.05, 0.2, generator=rng).item())
    
    for i in range(G):
        # Transcription Graph P(i)
        d_i: int = int(torch.randint(1, 4, (1,), generator=rng).item())
        tf_pool: Tensor = torch.tensor([tf for tf in TF_GENES if tf != i], dtype=torch.int64)
        idx: Tensor = torch.randperm(len(tf_pool), generator=rng)[:d_i]
        p_nodes: List[int] = tf_pool[idx].tolist()
        world.P_graph[i] = p_nodes

        world.a_ij[i] = {}
        for j in p_nodes:
            world.a_ij[i][j] = torch.empty(1, dtype=DTYPE).uniform_(0.5, 2.0, generator=rng).item()

        # Chromatin Graph E(i)
        epi_pool: Tensor = torch.tensor(EPI_GENES, dtype=torch.int64)
        idx_e: Tensor = torch.randperm(len(epi_pool), generator=rng)[:2]
        e_nodes: List[int] = epi_pool[idx_e].tolist()
        world.E_graph[i] = e_nodes

        world.beta_ij[i] = {}
        for j in e_nodes:
            world.beta_ij[i][j] = torch.empty(1, dtype=DTYPE).normal_(0, 1.5, generator=rng).item()

    return world

# ==========================================
# Module Interface Contracts
# ==========================================
def normalize_protein(P: Tensor, world: World) -> Tensor:
    # return P / (torch.sum(P) + world.epsilon)
    ### 0305 方案3：思路与方案1相同，但改为在生成tilde_P就乘G
    ### 以便update_chromatin时避免尺度失衡
    return P / (torch.mean(P) + world.epsilon)
    ### 采用此方案

def compute_TFinput(tilde_P: Tensor, world: World) -> Tensor:
    TFinput: Tensor = torch.zeros(G, dtype=DTYPE)
    for i in range(G):
        d_i = len(world.P_graph[i])
        prod = 1.0
        for j in world.P_graph[i]:
            # ### 0305 方案1：乘以 G，将均值拉回 1.0 附近，抵消网络规模带来的稀释，其余不变
            # scaled_P = tilde_P[j] * G
            # prod *= (scaled_P ** world.a_ij[i][j])
            # ### 这样会使下游的Z计算受影响，导致Z几乎不变，故弃用
            prod *= (tilde_P[j] ** world.a_ij[i][j])
        TFinput[i] = prod ** (1.0 / d_i)
    return TFinput

def update_chromatin(tilde_P: Tensor, world: World) -> Tensor:
    Z_next: Tensor = torch.zeros(G, dtype=DTYPE)
    for i in range(G):
        epi_sum: float = sum(world.beta_ij[i][j] * tilde_P[j] for j in world.E_graph[i])
        Z_next[i] = world.alpha[i] + epi_sum
    return stable_sigmoid(Z_next)

def update_mRNA(X: Tensor, Z: Tensor, TFinput: Tensor, world: World) -> Tensor:
    hill_term = (TFinput ** world.n) / (world.K ** world.n + TFinput ** world.n)
    X_next = (1.0 - world.delta_x) * X + Z * world.rho * hill_term

    ### 0305 ADD
    print(f'Hill term avg: {hill_term.mean().item()}')
    ###
    return X_next

def update_protein_raw(P: Tensor, X: Tensor, world: World) -> Tensor:
    P_next = (1.0 - world.delta_p) * P + world.gamma * X
    return P_next

def apply_resource_projection(P_raw: Tensor, world: World) -> Tensor:
    total_P = torch.sum(P_raw)
    if total_P > world.R_total:
        return P_raw * (world.R_total / total_P)
    return P_raw

def update_fate(N: float, world: World) -> float:
    return N + world.r * N * (1.0 - N / world.K_pop)

# ==========================================
# Main Time Loop
# ==========================================
def simulate_single_cell(world: World, X0: Tensor, P0: Tensor, Z0: Tensor, N0: float, t_steps: int = T) -> Dict[str, Tensor]:
    """
    Returns:
        X_traj: Tensor[t_steps+1, G]
        P_traj: Tensor[t_steps+1, G]
        Z_traj: Tensor[t_steps+1, G]
        N_traj: Tensor[t_steps+1, ]
    Convention: index 0 stores initial state at t=0.
    """

    X, P, Z, N = X0, P0, Z0, N0

    X_traj: Tensor = torch.zeros((t_steps+1, G), dtype=DTYPE)
    P_traj: Tensor = torch.zeros((t_steps+1, G), dtype=DTYPE)
    Z_traj: Tensor = torch.zeros((t_steps+1, G), dtype=DTYPE)
    N_traj: Tensor = torch.zeros((t_steps+1,), dtype=DTYPE)

    X_traj[0], P_traj[0], Z_traj[0], N_traj[0] = X, P, Z, N

    for t in range(t_steps):
        tilde_P: Tensor = normalize_protein(P, world)
        TFinput: Tensor = compute_TFinput(tilde_P, world)
        Z_next: Tensor = update_chromatin(tilde_P, world)
        X_next: Tensor = update_mRNA(X, Z, TFinput, world)
        P_raw: Tensor = update_protein_raw(P, X, world)
        P_next: Tensor = apply_resource_projection(P_raw, world)
        N_next: float = update_fate(N, world)

        X, P, Z, N = X_next, P_next, Z_next, N_next
        X_traj[t+1], P_traj[t+1], Z_traj[t+1], N_traj[t+1] = X, P, Z, N
    return  {'X_traj': X_traj, 'P_traj': P_traj, 'Z_traj': Z_traj, 'N_traj': N_traj}

# ==========================================
# Monte Carlo Wrappers (scRNA mode)
# ==========================================
def sample_initial_state(cell_seed: int, world: World) -> Tuple[Tensor, Tensor, Tensor, float]:
    rng: torch.Generator = torch.Generator()
    rng.manual_seed(cell_seed)

    X0: Tensor = torch.empty(G, dtype=DTYPE).uniform_(0, 1, generator=rng)
    P0_raw: Tensor = world.gamma * X0
    P0: Tensor = apply_resource_projection(P0_raw, world)
    Z0: Tensor = stable_sigmoid(world.alpha)
    N0: float = 1.0

    return X0, P0, Z0, N0

def run_simulation(seed: int, save_path: str = None) -> Dict[str, Tensor]:
    world_seed: int = seed
    cell_seed: int = seed + 1

    world: World = sample_world(world_seed)
    X0, P0, Z0, N0 = sample_initial_state(cell_seed, world)

    traj: Dict[str, Tensor] = simulate_single_cell(world, X0, P0, Z0, N0, T)

    if save_path is not None:
        torch.save({
            'X_traj': traj['X_traj'],
            'P_traj': traj['P_traj'],
            'Z_traj': traj['Z_traj'],
            'N_traj': traj['N_traj'],
            'world': world.to_dict(),
        }, save_path)

    return traj

def generate_dataset(world_seed: int, M: int, save_path: str = None) -> Tuple[Tensor, World]:
    world: World = sample_world(world_seed)
    C: Tensor = torch.zeros((M, G), dtype=DTYPE)

    for c in range(M):
        cell_seed_c: int = world_seed + 1 + c
        X0, P0, Z0, N0 = sample_initial_state(cell_seed_c, world)
        traj: Dict[str, Tensor] = simulate_single_cell(world, X0, P0, Z0, N0, T)
        C[c, :] = traj['X_traj'][T, :]

    if save_path is not None:
        torch.save({
            'expression': C,
            'world': world.to_dict(),
        }, save_path)

    return C, world

# ==========================================
# Perturbation Interface
# ==========================================
def apply_perturbation(world: World, state: State, config: Dict[str, Any]) -> Tuple[World, State]:
    world_perturbed: World = copy.deepcopy(world)
    state_perturbed: State = copy.deepcopy(state)

    if 'knockout' in config:
        for i in config['knockout']:
            state_perturbed['X'][i] = 0.0
    
    if 'override_rho' in config:
        for i, val in config['override_rho'].items():
            world_perturbed.rho[i] = float(val)

    if 'override_a_ij' in config:
        for i, j, val in config['override_a_ij']:
            if j in world_perturbed.a_ij[i]:
                world_perturbed.a_ij[i][j] = float(val)
    
    if 'override_alpha' in config:
        for i, val in config['override_alpha'].items():
            world_perturbed.alpha[i] = float(val)

    if 'R_total' in config:
        world_perturbed.R_total = float(config['R_total'])

    return world_perturbed, state_perturbed


# ==========================================
# Execution Sanity Tests
# ==========================================
def run_smoke_test(seed: int, T: int = 10) -> Dict[str, Tensor]:
    """
    Run a quick smoke test with T time steps.
    
    Tests basic functionality without full validation.
    """
    print(f"Running smoke test with T={T}...")
    
    world: World = sample_world(seed)
    cell_seed: int = seed + 1
    X0, P0, Z0, N0 = sample_initial_state(cell_seed, world)
    
    traj: Dict[str, Tensor] = simulate_single_cell(world, X0, P0, Z0, N0, T)
    
    print(f"  X shape: {traj['X_traj'].shape}")
    print(f"  P shape: {traj['P_traj'].shape}")
    print(f"  Z shape: {traj['Z_traj'].shape}")
    print(f"  N shape: {traj['N_traj'].shape}")
    print("Smoke test passed!")
    
    return traj


def run_sanity_tests(seed: int) -> None:
    """
    Run sanity tests to verify simulation correctness.
    
    Tests:
        1. Reproducibility: same seed produces identical results
        2. Non-negativity: X, P >= 0, Z in [0, 1]
        3. Resource bound: sum(P) <= R_total
        4. Stability: all values finite for T=200
    """
    print("Running sanity tests...")
    
    run1: Dict[str, Tensor] = run_simulation(seed)
    run2: Dict[str, Tensor] = run_simulation(seed)
    
    assert torch.equal(run1['X_traj'], run2['X_traj']), "Reproducibility check failed."
    print('Reproducibility check passed.')
    
    assert torch.all(run1['X_traj'] >= 0), "X must be >= 0."
    assert torch.all(run1['P_traj'] >= 0), "P must be >= 0."
    assert torch.all((run1['Z_traj'] >= 0) & (run1['Z_traj'] <= 1)), "Z must be in [0, 1]."    
    assert torch.all(torch.sum(run1['P_traj'], dim=1) <= R_TOTAL + EPSILON), "Resource limit exceeded."
    print('Non-negativity and Resource bound checks passed.')
    
    assert torch.all(torch.isfinite(run1['X_traj'])), "X must be finite."
    assert torch.all(torch.isfinite(run1['P_traj'])), "P must be finite."
    assert torch.all(torch.isfinite(run1['Z_traj'])), "Z must be finite."
    assert torch.all(torch.isfinite(run1['N_traj'])), "N must be finite."
    print('Stability check passed.')
    
    print("All sanity tests passed!")


if __name__ == "__main__":
    os.chdir(os.path.dirname(os.path.abspath(__file__)) or '.')
    os.makedirs('./data', exist_ok=True)
    print("Running DDC Phase 0 Standard Pipeline...")
    SEED: int = 42
    
    print("\n--- T=10 Smoke Test ---")
    run_smoke_test(SEED, T=10)
    
    # print("\n--- T=200 Stability Test ---")
    # run_sanity_tests(SEED)
    
    # print(f"\n--- Multi-cell Dataset, random seed: {SEED} ---")
    # M_cells: int = 100
    # multi_cell_sim_path = './data/multi_cell_trajectory.pt'
    # dataset, world = generate_dataset(SEED, M_cells, save_path = multi_cell_sim_path)
    # print(f'Dataset generated successfully: {dataset.shape}, data saved to {multi_cell_sim_path}')

    # print(f'\n--- Single-cell Dataset, random seed: {SEED} ---')
    # sim_path = './data/trajectory.pt'
    # result_dict = run_simulation(SEED, save_path = sim_path)
    # print(f'Simulation completed. Trajectory saved to {sim_path}')