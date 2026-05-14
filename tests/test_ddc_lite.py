"""
DDC Lite Level 0 — Sanity Tests
================================

Tests for src/ddc_lite.py per:
  docs/01_DDC_Lite_Curriculum/Level_0_Minimal_TF/02_Implementation_Guidance.md §12
  docs/01_DDC_Lite_Curriculum/Level_0_Minimal_TF/03_Acceptance_Criteria.md §11

Author: zhanghl
"""
import os
import sys
import json
import torch

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import ddc_lite as dl


SEED: int = 42


def test_world_reproducibility():
    """same seed -> identical World"""
    w1 = dl.sample_world(SEED)
    w2 = dl.sample_world(SEED)
    d1 = json.dumps(w1.to_dict(), sort_keys=True)
    d2 = json.dumps(w2.to_dict(), sort_keys=True)
    assert d1 == d2

    w3 = dl.sample_world(SEED + 1)
    d3 = json.dumps(w3.to_dict(), sort_keys=True)
    assert d1 != d3


def test_trajectory_reproducibility():
    """same seed -> identical trajectory"""
    r1 = dl.run_simulation(SEED)
    r2 = dl.run_simulation(SEED)
    assert torch.equal(r1['X_traj'], r2['X_traj'])
    assert torch.equal(r1['P_traj'], r2['P_traj'])


def test_output_shapes():
    """X_traj: (T+1, G),  P_traj: (T+1, G)"""
    traj = dl.run_simulation(SEED)
    assert traj['X_traj'].shape == (dl.T + 1, dl.G)
    assert traj['P_traj'].shape == (dl.T + 1, dl.G)


def test_nonnegative_X_P():
    traj = dl.run_simulation(SEED)
    assert torch.all(traj['X_traj'] >= 0)
    assert torch.all(traj['P_traj'] >= 0)


def test_activation_monotonicity():
    """higher P_j -> higher R_i for activation edges"""
    w = dl.sample_world(SEED)
    P_lo = torch.ones(dl.G, dtype=dl.DTYPE) * 0.5
    P_hi = torch.ones(dl.G, dtype=dl.DTYPE) * 5.0
    R_lo = dl.compute_regulatory_response(P_lo, w)
    R_hi = dl.compute_regulatory_response(P_hi, w)
    for i in range(dl.G):
        s = w.edge_signs[i][w.P_graph[i][0]]
        if s == dl.ACTIVATION:
            assert R_lo[i] <= R_hi[i], f"activation monotonicity failed at gene {i}"
        else:
            assert R_lo[i] >= R_hi[i], f"repression monotonicity failed at gene {i}"


def test_no_resource_projection():
    """sum(P) is not capped; no R_TOTAL constant"""
    traj = dl.run_simulation(SEED)
    P_sum = traj['P_traj'].sum(dim=1)
    assert float(P_sum[-1]) > 0.0
    assert not hasattr(dl, 'R_TOTAL')


def test_no_chromatin_dynamics():
    """no Z in trajectory or World"""
    traj = dl.run_simulation(SEED)
    assert 'Z' not in traj
    assert 'Z_traj' not in traj
    assert 'N_traj' not in traj

    w = dl.sample_world(SEED)
    X0, P0 = dl.sample_initial_state(SEED + 1, w)
    r = dl.simulate_single_cell(w, X0, P0, 5)
    assert 'Z_traj' not in r
    assert 'N_traj' not in r


def test_ground_truth_export():
    """world.to_dict() contains all required fields and is JSON-serializable"""
    w = dl.sample_world(SEED)
    d = w.to_dict()

    assert 'seed' in d
    assert 'P_graph' in d
    assert 'edge_signs' in d
    assert 'parameters' in d

    params = d['parameters']
    for key in ('rho', 'K', 'n', 'delta_x', 'delta_p', 'gamma'):
        assert key in params

    json_str = json.dumps(d)
    d2 = json.loads(json_str)
    assert d2['seed'] == d['seed']


def test_phase0_file_untouched():
    """src/ddc.py contains Phase 0 markers"""
    phase0_path = os.path.join(os.path.dirname(__file__), '..', 'src', 'ddc.py')
    with open(phase0_path, 'rb') as f:
        content = f.read()
    assert b'DDC Phase 0' in content
    assert b'v1.2' in content
    assert b'G: int = 50' in content


def test_no_phase0_mechanisms():
    """ddc_lite must not import or expose excluded Phase 0 mechanisms"""
    for name in ('normalize_protein', 'update_chromatin',
                 'apply_resource_projection', 'update_fate',
                 'compute_TFinput', 'stable_sigmoid'):
        assert not hasattr(dl, name), f"excluded: {name}"


if __name__ == '__main__':
    tests = [
        ('world_reproducibility',    test_world_reproducibility),
        ('trajectory_reproducibility', test_trajectory_reproducibility),
        ('output_shapes',            test_output_shapes),
        ('nonnegative_X_P',          test_nonnegative_X_P),
        ('activation_monotonicity',  test_activation_monotonicity),
        ('no_resource_projection',   test_no_resource_projection),
        ('no_chromatin_dynamics',    test_no_chromatin_dynamics),
        ('ground_truth_export',      test_ground_truth_export),
        ('phase0_file_untouched',    test_phase0_file_untouched),
        ('no_phase0_mechanisms',     test_no_phase0_mechanisms),
    ]

    passed = 0
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f'  PASS: {name}')
            passed += 1
        except AssertionError as e:
            print(f'  FAIL: {name} — {e}')
            failed += 1
        except Exception as e:
            print(f'  ERROR: {name} — {e}')
            failed += 1

    print(f'\n{passed}/{len(tests)} passed, {failed} failed')
    if failed > 0:
        sys.exit(1)
