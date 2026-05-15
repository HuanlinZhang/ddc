"""
DDC CLI
=======
Command-line interface for the DDC simulator.

Usage:
    ddc run --seed 42 --output ./runs/exp1/
    ddc dataset --world-seed 0 --M 100 --output dataset.pt
    ddc world --seed 42 --output world.json
    ddc smoke-test
"""

import argparse
import json
import os
import sys

from ddc.core import (
    run_simulation,
    generate_dataset,
    sample_world,
    run_smoke_test,
    run_sanity_tests,
)


def cmd_run(args):
    """Run single-cell simulation."""
    intervention_config = None
    if args.knockdown_gene is not None:
        intervention_config = {"knockdown_X": args.knockdown_gene}

    perturbation_config = None
    if args.knockout_gene is not None:
        perturbation_config = {"knockout": args.knockout_gene}

    traj = run_simulation(
        world_seed=args.seed,
        save_path=args.output,
        T=args.T,
        cell_seed=args.cell_seed,
        intervention_time=args.intervention_time,
        intervention_config=intervention_config,
        perturbation_config=perturbation_config,
        enable_resource_projection=not args.disable_resource_projection,
    )
    print(f"Trajectory saved to {args.output}")
    print(f"  X_traj: {traj['X_traj'].shape}")
    print(f"  P_traj: {traj['P_traj'].shape}")
    print(f"  Z_traj: {traj['Z_traj'].shape}")
    print(f"  N_traj: {traj['N_traj'].shape}")


def cmd_dataset(args):
    """Generate multi-cell dataset."""
    dataset, world = generate_dataset(
        world_seed=args.world_seed,
        M=args.M,
        save_path=args.output,
    )
    print(f"Dataset saved to {args.output}")
    print(f"  Shape: {dataset.shape}")
    print(f"  World seed: {world.seed}")


def cmd_world(args):
    """Sample and save a world (gene regulatory network)."""
    world = sample_world(args.seed)
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        json.dump(world.to_dict(), f, indent=2)
    print(f"World saved to {args.output}")
    print(f"  Gene count: 50")
    print(f"  World seed: {world.seed}")


def cmd_smoke_test(args):
    """Run smoke test."""
    run_smoke_test(seed=args.seed, T=args.T)
    print("Smoke test passed.")


def cmd_sanity(args):
    """Run sanity tests."""
    run_sanity_tests(seed=args.seed)
    print("All sanity tests passed.")


def main():
    parser = argparse.ArgumentParser(
        prog="ddc",
        description="DDC - Designed Digital Cell simulator",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ddc run
    p_run = sub.add_parser("run", help="Run single-cell simulation")
    p_run.add_argument("--seed", type=int, required=True, help="Random seed for gene network world")
    p_run.add_argument("--output", type=str, required=True, help="Output .pt file path")
    p_run.add_argument(
        "--T",
        type=int,
        default=200,
        help="Number of simulation timesteps (default: 200)",
    )
    p_run.add_argument(
        "--cell-seed",
        type=int,
        default=None,
        help="Random seed for cell initial state (default: --seed + 1)",
    )
    p_run.add_argument(
        "--intervention-time",
        type=int,
        default=None,
        help="Apply intervention at this time step",
    )
    p_run.add_argument(
        "--knockdown-gene",
        type=int,
        nargs="+",
        default=None,
        help="Gene indices to knockdown (state-level, single-step)",
    )
    p_run.add_argument(
        "--knockout-gene",
        type=int,
        nargs="+",
        default=None,
        help="Gene indices to knockout (parameter-level, persistent, rho_i=0)",
    )
    p_run.add_argument(
        "--disable-resource-projection",
        action="store_true",
        default=False,
        help="Disable resource constraint (ΣP ≤ R_total)",
    )
    p_run.set_defaults(func=cmd_run)

    # ddc dataset
    p_ds = sub.add_parser("dataset", help="Generate multi-cell dataset")
    p_ds.add_argument("--world-seed", type=int, required=True, help="World random seed")
    p_ds.add_argument("--M", type=int, required=True, help="Number of cells")
    p_ds.add_argument("--output", type=str, required=True, help="Output .pt file path")
    p_ds.set_defaults(func=cmd_dataset)

    # ddc world
    p_w = sub.add_parser("world", help="Sample and save a world")
    p_w.add_argument("--seed", type=int, required=True, help="Random seed")
    p_w.add_argument("--output", type=str, required=True, help="Output .json file path")
    p_w.set_defaults(func=cmd_world)

    # ddc smoke-test
    p_st = sub.add_parser("smoke-test", help="Run smoke test (T=10)")
    p_st.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p_st.add_argument("--T", type=int, default=10, help="Time steps (default: 10)")
    p_st.set_defaults(func=cmd_smoke_test)

    # ddc sanity
    p_san = sub.add_parser("sanity", help="Run sanity tests (T=200)")
    p_san.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    p_san.set_defaults(func=cmd_sanity)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
