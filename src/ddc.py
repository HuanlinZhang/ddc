"""
DDC Phase 0 - Designed Digital-Cell Model
==========================================

A computational simulation framework for gene regulatory network dynamics.
Based on the Phase 0 specification documents.

Author: zhanghl
Version: v1.1.1
Status: Active

================================================================================
Version History & Update Notes
================================================================================

Note: v1.1.0 and v1.1.1 are parallel branches, each modifying only the stated
change on top of v1.0. All other code remains identical to v1.0.

v1.1.1 (2026-03-06):
    - File: ddc.py
    - Modified K parameter sampling: K values are now divided by G (gene count)
    - Previous:  world.K = uniform(0.1, 1.0)
    - Current:   world.K = uniform(0.1, 1.0) / G
    - Reason:    To address convergence issue. The Hill function half-saturation
      constant K needs to be scaled relative to network size for proper
      activation threshold crossing.

v1.1.0 (2026-03-06):
    - File: ddc_v_1_1_0.py
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

# ==========================================
# Network Sampling Functions
# ==========================================
def sample_network_graphs(seed: int) -> Tuple[Dict[int, List[int]], Dict[int, List[int]]]:
    """
    Sample P_graph and E_graph for the gene regulatory network.
    
    P_graph: TF (0-5) → All genes (excluding self)
    E_graph: EPI (17-19) → All genes
    """
    rng = torch.Generator()
    rng.manual_seed(seed)
    
    P_graph = {}
    E_graph = {}
    
    # P_graph: TF (0-5) regulate all genes except themselves
    for i in range(6):  # TF genes
        targets = list(range(G))
        targets.remove(i)  # Exclude self
        
        # Sample 1-3 targets per TF
        k = torch.randint(1, 4, (1,), generator=rng).item()
        selected_targets = torch.randperm(len(targets), generator=rng)[:k].tolist()
        P_graph[i] = [targets[j] for j in selected_targets]
    
    # E_graph: EPI (17-19) regulate all genes
    for i in range(17, 20):  # EPI genes
        # Sample 2 targets per EPI
        targets = list(range(G))
        k = 2
        selected_targets = torch.randperm(len(targets), generator=rng)[:k].tolist()
        E_graph[i] = [targets[j] for j in selected_targets]
    
    return P_graph, E_graph

# ==========================================
# Parameter Sampling Functions
# ==========================================
def sample_gene_parameters(seed: int) -> Tuple[Tensor, Tensor, Tensor, Tensor, Tensor, Tensor, Tensor]:
    """
    Sample gene-level parameters.
    
    Returns: alpha, rho, K, n, delta_x, delta_p, gamma
    """
    rng = torch.Generator()
    rng.manual_seed(seed)
    
    # alpha: basal chromatin activation ~ N(0, 1)
    alpha = torch.normal(0.0, 1.0, (G,), generator=rng, dtype=DTYPE)
    
    # rho: mRNA production rate ~ U(0.5, 2.0)
    rho = torch.empty(G, dtype=DTYPE).uniform_(0.5, 2.0, generator=rng)
    
    # K: Hill coefficient denominator ~ U(0.1, 1.0) / G
    K = torch.empty(G, dtype=DTYPE).uniform_(0.1, 1.0, generator=rng) / G
    
    # n: Hill coefficient (fixed at 2.0)
    n = torch.full((G,), 2.0, dtype=DTYPE)
    
    # delta_x: mRNA degradation rate ~ U(0.1, 0.5)
    delta_x = torch.empty(G, dtype=DTYPE).uniform_(0.1, 0.5, generator=rng)
    
    # delta_p: protein degradation rate ~ U(0.05, 0.3)
    delta_p = torch.empty(G, dtype=DTYPE).uniform_(0.05, 0.3, generator=rng)
    
    # gamma: translation rate (fixed at 1.0)
    gamma = torch.full((G,), 1.0, dtype=DTYPE)
    
    return alpha, rho, K, n, delta_x, delta_p, gamma

def sample_edge_parameters(seed: int, P_graph: Dict[int, List[int]], E_graph: Dict[int, List[int]]) -> Tuple[Dict[int, Dict[int, float]], Dict[int, Dict[int, float]]]:
    """
    Sample edge-level parameters.
    
    Returns: a_ij (TF regulation strength), beta_ij (epigenetic regulation strength)
    """
    rng = torch.Generator()
    rng.manual_seed(seed)
    
    a_ij = {}
    beta_ij = {}
    
    # a_ij: TF regulation strength ~ U(0.5, 2.0)
    for i, targets in P_graph.items():
        a_ij[i] = {}
        for j in targets:
            a_ij[i][j] = torch.empty(1, dtype=DTYPE).uniform_(0.5, 2.0, generator=rng).item()
    
    # beta_ij: epigenetic regulation strength ~ N(0, 1.5)
    for i, targets in E_graph.items():
        beta_ij[i] = {}
        for j in targets:
            beta_ij[i][j] = torch.normal(0.0, 1.5, (1,), generator=rng, dtype=DTYPE).item()
    
    return a_ij, beta_ij

def sample_run_parameters(seed: int) -> float:
    """
    Sample run-level parameters.
    
    Returns: r (cell growth rate) ~ U(0.05, 0.2)
    """
    rng = torch.Generator()
    rng.manual_seed(seed)
    
    r = torch.empty(1, dtype=DTYPE).uniform_(0.05, 0.2, generator=rng).item()
    return r

# ==========================================
# World Construction Function
# ==========================================
def sample_world(seed: int) -> World:
    """
    Generate a complete World object with all parameters sampled.
    """
    world = World(seed)
    
    # Sample network topology
    P_graph, E_graph = sample_network_graphs(seed)
    world.P_graph = P_graph
    world.E_graph = E_graph
    
    # Sample gene-level parameters
    alpha, rho, K, n, delta_x, delta_p, gamma = sample_gene_parameters(seed)
    world.alpha = alpha
    world.rho = rho
    world.K = K
    world.n = n
    world.delta_x = delta_x
    world.delta_p = delta_p
    world.gamma = gamma
    
    # Sample edge-level parameters
    a_ij, beta_ij = sample_edge_parameters(seed, P_graph, E_graph)
    world.a_ij = a_ij
    world.beta_ij = beta_ij
    
    # Sample run-level parameters
    world.r = sample_run_parameters(seed)
    
    return world

# ==========================================
# Core Simulation Functions
# ==========================================
def normalize_protein(P: Tensor, world: World) -> Tensor:
    """
    Normalize protein concentrations to sum to 1.
    
    Returns: tilde_P = P / (sum(P) + epsilon)
    """
    return P / (torch.sum(P) + world.epsilon)

def compute_TF_input(tilde_P: Tensor, world: World, gene: int) -> float:
    """
    Compute TF input for a given gene using geometric mean.
    
    TF_input = (∏ tilde_P_j^a_ij)^(1/d_i)
    where d_i is the degree of gene i in P_graph
    """
    if gene in world.P_graph:
        regulators = world.P_graph[gene]
        if not regulators:
            return 1.0
        
        product = 1.0
        for j in regulators:
            if j in world.a_ij[gene]:
                a_ij = world.a_ij[gene][j]
                product *= tilde_P[j] ** a_ij
        
        degree = len(regulators)
        return product ** (1.0 / degree) if degree > 0 else 1.0
    else:
        return 1.0

def hill_function(x: float, K: float, n: float) -> float:
    """
    Hill function: x^n / (K^n + x^n)
    """
    return (x ** n) / ((K ** n) + (x ** n))

def update_chromatin(tilde_P: Tensor, world: World) -> Tensor:
    """
    Update chromatin states using sigmoid activation.
    
    Z = σ(alpha + Σ beta_ij * tilde_P_j)
    """
    Z_input = world.alpha.clone()
    
    for i in range(G):
        if i in world.E_graph:
            for j in world.E_graph[i]:
                if j in world.beta_ij[i]:
                    Z_input[i] += world.beta_ij[i][j] * tilde_P[j]
    
    return stable_sigmoid(Z_input)

def update_mRNA(X: Tensor, Z: Tensor, tilde_P: Tensor, world: World) -> Tensor:
    """
    Update mRNA concentrations.
    
    X' = (1 - delta_x) * X + Z * rho * hill(TF_input)
    """
    X_new = torch.zeros_like(X)
    
    for i in range(G):
        TF_input = compute_TF_input(tilde_P, world, i)
        hill_val = hill_function(TF_input, world.K[i].item(), world.n[i].item())
        
        X_new[i] = (1.0 - world.delta_x[i]) * X[i] + Z[i] * world.rho[i] * hill_val
    
    return X_new

def update_protein_raw(P: Tensor, X: Tensor, world: World) -> Tensor:
    """
    Update protein concentrations without resource constraint.
    
    P_raw = (1 - delta_p) * P + gamma * X
    """
    return (1.0 - world.delta_p) * P + world.gamma * X

def apply_resource_projection(P_raw: Tensor, world: World) -> Tensor:
    """
    Apply resource constraint to protein concentrations.
    
    If sum(P_raw) > R_total, scale proportionally.
    """
    total_protein = torch.sum(P_raw)
    
    if total_protein > world.R_total:
        return P_raw * (world.R_total / total_protein)
    else:
        return P_raw

def update_cell_fate(N: float, world: World) -> float:
    """
    Update cell count using logistic growth.
    
    N' = N + r * N * (1 - N / K_pop)
    """
    return N + world.r * N * (1.0 - N / world.K_pop)

def simulate_single_cell(world: World, X0: Tensor, P0: Tensor, Z0: Tensor, N0: float, t_steps: int = T) -> Dict[str, Tensor]:
    """
    Simulate single cell time evolution.
    
    Returns: trajectories for X, P, Z, N
    """
    # Initialize trajectories
    X_traj = torch.zeros(t_steps + 1, G, dtype=DTYPE)
    P_traj = torch.zeros(t_steps + 1, G, dtype=DTYPE)
    Z_traj = torch.zeros(t_steps + 1, G, dtype=DTYPE)
    N_traj = torch.zeros(t_steps + 1, dtype=DTYPE)
    
    # Set initial state
    X_traj[0] = X0
    P_traj[0] = P0
    Z_traj[0] = Z0
    N_traj[0] = N0
    
    # Time evolution
    for t in range(t_steps):
        X = X_traj[t]
        P = P_traj[t]
        Z = Z_traj[t]
        N = N_traj[t]
        
        # Normalize protein
        tilde_P = normalize_protein(P, world)
        
        # Update chromatin
        Z_new = update_chromatin(tilde_P, world)
        
        # Update mRNA
        X_new = update_mRNA(X, Z_new, tilde_P, world)
        
        # Update protein (raw)
        P_raw = update_protein_raw(P, X_new, world)
        
        # Apply resource constraint
        P_new = apply_resource_projection(P_raw, world)
        
        # Update cell fate
        N_new = update_cell_fate(N, world)
        
        # Store new state
        X_traj[t + 1] = X_new
        P_traj[t + 1] = P_new
        Z_traj[t + 1] = Z_new
        N_traj[t + 1] = N_new
    
    return {
        'X_traj': X_traj,
        'P_traj': P_traj,
        'Z_traj': Z_traj,
        'N_traj': N_traj
    }

# ==========================================
# Initial State Sampling
# ==========================================
def sample_initial_state(cell_seed: int, world: World) -> Tuple[Tensor, Tensor, Tensor, float]:
    """
    Sample initial state for a single cell.
    """
    rng = torch.Generator()
    rng.manual_seed(cell_seed)
    
    # X0: random initial mRNA ~ U(0, 1)
    X0 = torch.empty(G, dtype=DTYPE).uniform_(0.0, 1.0, generator=rng)
    
    # P0: initial protein = gamma * X0 with resource projection
    P0_raw = world.gamma * X0
    P0 = apply_resource_projection(P0_raw, world)
    
    # Z0: chromatin initial state = sigmoid(alpha)
    Z0 = stable_sigmoid(world.alpha)
    
    # N0: initial cell count = 1.0
    N0 = 1.0
    
    return X0, P0, Z0, N0

# ==========================================
# High-level API Functions
# ==========================================
def run_simulation(seed: int, save_path: str = None) -> Dict[str, Tensor]:
    """
    Run single-cell simulation with given seed.
    
    Returns: trajectories for X, P, Z, N
    """
    # Sample world and initial state
    world = sample_world(seed)
    X0, P0, Z0, N0 = sample_initial_state(seed + 1, world)
    
    # Run simulation
    result = simulate_single_cell(world, X0, P0, Z0, N0)
    
    # Save if path provided
    if save_path:
        data_to_save = {
            'X_traj': result['X_traj'],
            'P_traj': result['P_traj'],
            'Z_traj': result['Z_traj'],
            'N_traj': result['N_traj'],
            'world': world.to_dict()
        }
        torch.save(data_to_save, save_path)
    
    return result

def generate_dataset(world_seed: int, M: int, save_path: str = None) -> Tuple[Tensor, World]:
    """
    Generate multi-cell dataset (scRNA-seq style).
    
    Returns: expression matrix (M, G), World object
    """
    # Sample world
    world = sample_world(world_seed)
    
    # Generate M cells
    expression_matrix = torch.zeros(M, G, dtype=DTYPE)
    
    for m in range(M):
        cell_seed = world_seed + m + 1
        X0, P0, Z0, N0 = sample_initial_state(cell_seed, world)
        
        # Run short simulation to get final state
        result = simulate_single_cell(world, X0, P0, Z0, N0, t_steps=10)
        
        # Use final mRNA as expression
        expression_matrix[m] = result['X_traj'][-1]
    
    # Save if path provided
    if save_path:
        data_to_save = {
            'expression': expression_matrix,
            'world': world.to_dict()
        }
        torch.save(data_to_save, save_path)
    
    return expression_matrix, world

def apply_perturbation(world: World, state: State, config: Dict[str, Any]) -> Tuple[World, State]:
    """
    Apply gene perturbation to world and state.
    
    Config options:
    - knockout: List of gene indices
    - override_rho: Dict {gene: value}
    - override_a_ij: List of (from, to, value)
    - override_alpha: Dict {gene: value}
    - R_total: float
    """
    world_perturbed = copy.deepcopy(world)
    state_perturbed = copy.deepcopy(state)
    
    # Knockout: set rho=0 for specified genes
    if 'knockout' in config:
        for gene in config['knockout']:
            world_perturbed.rho[gene] = 0.0
    
    # Override rho values
    if 'override_rho' in config:
        for gene, value in config['override_rho'].items():
            world_perturbed.rho[gene] = value
    
    # Override a_ij values
    if 'override_a_ij' in config:
        for from_gene, to_gene, value in config['override_a_ij']:
            if from_gene in world_perturbed.a_ij and to_gene in world_perturbed.a_ij[from_gene]:
                world_perturbed.a_ij[from_gene][to_gene] = value
    
    # Override alpha values
    if 'override_alpha' in config:
        for gene, value in config['override_alpha'].items():
            world_perturbed.alpha[gene] = value
    
    # Override R_total
    if 'R_total' in config:
        world_perturbed.R_total = config['R_total']
    
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
    Run comprehensive sanity tests.
    
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

# ==========================================
# Perturbation Helper Functions
# ==========================================
def perturb_world_and_state(world: World, state: State, config: Dict[str, Any]) -> Tuple[World, State]:
    """
    Apply perturbation to world and state based on config.
    
    Returns perturbed copies of world and state.
    """
    world_perturbed = copy.deepcopy(world)
    state_perturbed = copy.deepcopy(state)
    
    # Apply knockout
    if 'knockout' in config:
        for i, v in config['knockout'].items():
            world_perturbed.rho[i] = float(v)
    
    # Apply override_rho
    if 'override_rho' in config:
        for i, v in config['override_rho'].items():
            world_perturbed.rho[i] = float(v)
    
    # Apply override_a_ij
    if 'override_a_ij' in config:
        for i, j, v in config['override_a_ij']:
            if j in world_perturbed.a_ij[i]:
                world_perturbed.a_ij[i][j] = float(v)
    
    # Apply override_alpha
    if 'override_alpha' in config:
        for i, v in config['override_alpha'].items():
            world_perturbed.alpha[i] = float(v)
    
    # Apply R_total override
    if 'R_total' in config:
        world_perturbed.R_total = float(config['R_total'])
    
    return world_perturbed, state_perturbed

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
    
    # print(f"\n--- Single-cell Dataset, random seed: {SEED} ---")
    # sim_path = './data/trajectory.pt'
    # result_dict = run_simulation(SEED, save_path = sim_path)
    # print(f'Simulation completed. Trajectory saved to {sim_path}')