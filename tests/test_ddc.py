"""
Unit tests for DDC simulation engine.
Run with: pytest tests/
"""

import pytest
import torch
from ddc import (
    sample_world,
    sample_initial_state,
    simulate_single_cell,
    run_simulation,
    generate_dataset,
    apply_perturbation,
    apply_intervention,
    run_smoke_test,
    World,
    G,
    T,
    R_TOTAL,
    EPSILON,
    DTYPE,
)


class TestReproducibility:
    """Test that same seed produces identical results."""

    def test_run_simulation_reproducibility(self):
        traj1 = run_simulation(seed=42)
        traj2 = run_simulation(seed=42)
        for key in ["X_traj", "P_traj", "Z_traj", "N_traj"]:
            assert torch.equal(traj1[key], traj2[key]), f"{key} differs between runs"

    def test_sample_world_reproducibility(self):
        w1 = sample_world(seed=123)
        w2 = sample_world(seed=123)
        assert w1.to_dict() == w2.to_dict()


class TestNonNegativity:
    """Test that state variables stay within valid ranges."""

    def test_single_cell_trajectory(self):
        world = sample_world(seed=42)
        X0, P0, Z0, N0 = sample_initial_state(seed=42, world=world)
        traj = simulate_single_cell(world, X0, P0, Z0, N0, t_steps=200)

        assert torch.all(traj["X_traj"] >= 0), "X must be >= 0"
        assert torch.all(traj["P_traj"] >= 0), "P must be >= 0"
        assert torch.all((traj["Z_traj"] >= 0) & (traj["Z_traj"] <= 1)), "Z must be in [0, 1]"

    def test_dataset_cells(self):
        dataset, world = generate_dataset(world_seed=99, M=10)
        assert torch.all(dataset >= 0), "dataset values must be >= 0"


class TestResourceBound:
    """Test that protein resource constraint is respected."""

    def test_resource_projection(self):
        world = sample_world(seed=42)
        X0, P0, Z0, N0 = sample_initial_state(seed=42, world=world)
        traj = simulate_single_cell(world, X0, P0, Z0, N0, t_steps=200)

        for t in range(1, T + 1):
            total_P = torch.sum(traj["P_traj"][t])
            assert total_P <= R_TOTAL + EPSILON, f"Resource bound violated at t={t}: sum(P)={total_P}"


class TestStability:
    """Test that values remain finite over full simulation."""

    def test_full_trajectory_finite(self):
        traj = run_simulation(seed=99)
        for key in ["X_traj", "P_traj", "Z_traj", "N_traj"]:
            assert torch.all(torch.isfinite(traj[key])), f"{key} contains non-finite values"


class TestWorld:
    """Test World object construction and serialization."""

    def test_world_gene_count(self):
        world = sample_world(seed=42)
        assert world.alpha.shape[0] == G
        assert world.rho.shape[0] == G
        assert world.K.shape[0] == G
        assert world.n.shape[0] == G
        assert world.delta_x.shape[0] == G
        assert world.delta_p.shape[0] == G
        assert world.gamma.shape[0] == G

    def test_world_serialization(self):
        world = sample_world(seed=42)
        data = world.to_dict()
        assert "parameters" in data
        assert "gene_annotation" in data

        new_world = World(seed=0)
        new_world.from_dict(data)
        assert new_world.to_dict() == data


class TestInitialConditions:
    """Test initial state sampling."""

    def test_initial_Z_bounds(self):
        world = sample_world(seed=42)
        X0, P0, Z0, N0 = sample_initial_state(seed=42, world=world)
        assert torch.all((Z0 >= 0) & (Z0 <= 1)), "Z0 must be in [0, 1]"

    def test_initial_N(self):
        world = sample_world(seed=42)
        X0, P0, Z0, N0 = sample_initial_state(seed=42, world=world)
        assert N0 == 1.0, "N0 must be 1.0"


class TestPerturbation:
    """Test perturbation (parameter-level) interface."""

    def test_perturbation_creates_copy(self):
        world = sample_world(seed=42)
        state = (torch.zeros(G), torch.zeros(G), torch.zeros(G), 1.0)
        world2, state2 = apply_perturbation(world, state, {"knockout": [0]})
        assert world is not world2, "apply_perturbation must return a new world copy"
        assert world.rho[0] != 0, "original world must not be mutated"

    def test_knockout_effect(self):
        world = sample_world(seed=42)
        state = (torch.zeros(G), torch.zeros(G), torch.zeros(G), 1.0)
        world2, _ = apply_perturbation(world, state, {"knockout": [0]})
        assert world2.rho[0] == 0.0
        assert world.rho[0] != 0.0


class TestIntervention:
    """Test intervention (state-level) interface."""

    def test_intervention_knockdown(self):
        state = (torch.ones(G), torch.ones(G), torch.ones(G) * 0.5, 1.0)
        X, P, Z, N = state
        new_state = apply_intervention(state, {"knockdown_X": [0]})
        assert new_state[0][0] == 0.0, "X[0] must be set to 0"
        assert new_state[1][0] == 1.0, "P[0] unchanged"
        assert state[0][0] == 1.0, "original state must not be mutated"

    def test_intervention_set_P(self):
        state = (torch.ones(G), torch.ones(G), torch.ones(G) * 0.5, 1.0)
        new_state = apply_intervention(state, {"set_P": [(1, 2.5)]})
        assert new_state[1][1] == 2.5


class TestSmokeTest:
    """Test that smoke test function runs without error."""

    def test_smoke_test(self):
        traj = run_smoke_test(seed=42, T=10)
        assert traj["X_traj"].shape == (11, 50)
        assert traj["P_traj"].shape == (11, 50)
        assert traj["Z_traj"].shape == (11, 50)
        assert traj["N_traj"].shape == (11,)
