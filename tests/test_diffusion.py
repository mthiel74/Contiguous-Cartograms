"""Basic correctness tests for :class:`DiffusionSolver`."""

import numpy as np
import pytest

from cartogram import DiffusionSolver


def _random_rho(rng, ny=32, nx=48):
    rho = 1.0 + 0.5 * rng.standard_normal((ny, nx))
    return np.clip(rho, 0.1, None)


def test_t0_returns_input_density():
    rng = np.random.default_rng(0)
    rho0 = _random_rho(rng)
    solver = DiffusionSolver(rho0, bbox=(0, 0, 2.0, 1.0))
    assert np.allclose(solver.density_at(0.0), rho0, atol=1e-10)


def test_mean_is_conserved_for_all_t():
    """The heat equation on a Neumann domain conserves the spatial mean."""
    rng = np.random.default_rng(1)
    rho0 = _random_rho(rng)
    solver = DiffusionSolver(rho0, bbox=(0, 0, 1.0, 1.0))
    mean0 = rho0.mean()
    for t in [0.0, 0.01, 0.1, 1.0, 10.0]:
        assert np.isclose(solver.density_at(t).mean(), mean0, atol=1e-10)


def test_density_approaches_uniform_for_large_t():
    rng = np.random.default_rng(2)
    rho0 = _random_rho(rng)
    solver = DiffusionSolver(rho0, bbox=(0, 0, 1.0, 1.0))
    rho_inf = solver.density_at(solver.convergence_time(tol=1e-6))
    assert rho_inf.std() < 1e-4


def test_nonnegative_input_required():
    rho = np.ones((8, 8))
    rho[0, 0] = -0.01
    with pytest.raises(ValueError):
        DiffusionSolver(rho, bbox=(0, 0, 1, 1))


def test_gradient_has_correct_shape_and_vanishes_at_infinity():
    rng = np.random.default_rng(3)
    rho0 = _random_rho(rng)
    solver = DiffusionSolver(rho0, bbox=(0, 0, 1.0, 1.0))
    rho, gx, gy = solver.density_and_gradient_at(solver.convergence_time(1e-6))
    assert rho.shape == gx.shape == gy.shape == rho0.shape
    assert np.max(np.abs(gx)) < 1e-3
    assert np.max(np.abs(gy)) < 1e-3
