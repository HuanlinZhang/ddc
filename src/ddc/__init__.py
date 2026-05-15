"""
DDC - Designed Digital Cell
===========================
A computational framework for simulating gene regulatory network dynamics.

Usage:
    from ddc import run_simulation, generate_dataset, sample_world
"""

from ddc.core import (
    sample_world,
    sample_initial_state,
    simulate_single_cell,
    run_simulation,
    generate_dataset,
    apply_perturbation,
    apply_intervention,
    run_smoke_test,
    run_sanity_tests,
    run_intervention_sanity_test,
    World,
    G,
    T,
    R_TOTAL,
    DTYPE,
    ENABLE_RESOURCE_PROJECTION,
)

__version__ = "0.1.0"
__author__ = "zhanghl"

__all__ = [
    "sample_world",
    "sample_initial_state",
    "simulate_single_cell",
    "run_simulation",
    "generate_dataset",
    "apply_perturbation",
    "apply_intervention",
    "run_smoke_test",
    "run_sanity_tests",
    "run_intervention_sanity_test",
    "World",
    "G",
    "T",
    "R_TOTAL",
    "DTYPE",
    "ENABLE_RESOURCE_PROJECTION",
]
