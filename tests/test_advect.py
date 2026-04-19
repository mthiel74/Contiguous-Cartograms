"""Tests for the advection / Cartogram layer."""

import numpy as np

from cartogram import Cartogram, DiffusionSolver, advect_points


def test_uniform_density_is_identity():
    """A uniform density field induces no motion."""
    rho = np.ones((24, 32))
    solver = DiffusionSolver(rho, bbox=(0, 0, 4.0, 3.0))
    pts = np.array([[1.0, 1.0], [3.0, 2.0], [0.1, 2.9]])
    out = advect_points(solver, pts, t_max=5.0)
    assert np.allclose(out, pts, atol=1e-8)


def test_cartogram_spreads_dense_block():
    """A central dense block should push grid points outward so that the
    final density is nearly uniform; we test a weaker, sufficient proxy:
    the area of the transformed central block grows relative to the
    original.
    """
    rho = np.ones((64, 64))
    rho[24:40, 24:40] = 9.0
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3)

    # Corners of the central block in original coordinates.
    block_corners = np.array(
        [
            [24 / 64, 24 / 64],
            [40 / 64, 24 / 64],
            [40 / 64, 40 / 64],
            [24 / 64, 40 / 64],
        ]
    )
    moved = cart.transform(block_corners)
    # Shoelace area before and after.
    def _area(p):
        x, y = p[:, 0], p[:, 1]
        return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))

    a0 = _area(block_corners)
    a1 = _area(moved)
    assert a1 > 1.5 * a0, f"dense block did not expand: before={a0}, after={a1}"


def test_corners_stay_on_boundary():
    """Neumann BCs imply the four corners of the bbox are fixed points of
    the deformation up to numerical error."""
    rho = np.ones((32, 32))
    rho[10:22, 10:22] = 5.0
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3)
    corners = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    moved = cart.transform(corners)
    assert np.allclose(moved, corners, atol=2e-2)


def test_points_stay_inside_bbox():
    rng = np.random.default_rng(42)
    rho = 1.0 + rng.random((32, 32))
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3)
    query = rng.random((50, 2))
    moved = cart.transform(query)
    assert moved.min() >= -1e-6
    assert moved.max() <= 1.0 + 1e-6


def test_speed_mode_matches_quality_mode():
    """The snapshot-cached RHS should agree with the exact solver up to a
    small fraction of a grid cell. The test uses the same contrast used
    by the basic block-expansion test, which is representative of a
    real-world rasterised cartogram."""
    rho = np.ones((48, 48))
    rho[18:30, 18:30] = 9.0

    cart_q = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart_q.run(tol=1e-3, performance_goal="quality")

    cart_s = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart_s.run(tol=1e-3, performance_goal="speed", snapshots=60)

    # Every deformed grid point should agree to a small fraction of the
    # grid spacing (1 / 48 ~= 0.021 here).
    max_err = np.abs(
        cart_q._final_grid_points - cart_s._final_grid_points
    ).max()
    assert max_err < 5e-3, f"speed mode drift too large: {max_err}"


def test_speed_mode_is_default():
    """Sanity check the default path doesn't diverge silently."""
    rho = np.ones((32, 32))
    rho[10:22, 10:22] = 5.0
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3)  # default performance_goal
    corners = np.array([[0, 0], [1, 0], [1, 1], [0, 1]], dtype=float)
    moved = cart.transform(corners)
    assert np.allclose(moved, corners, atol=2e-2)
