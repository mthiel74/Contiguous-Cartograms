"""Integrate points through the Gastner–Newman velocity field.

Given a :class:`DiffusionSolver`, the velocity field at time ``t`` is

    v(x, t) = - ∇ρ(x, t) / ρ(x, t),

sampled on the grid and bilinearly interpolated elsewhere. The cartogram
positions are the limit as t → ∞ of the advected trajectories

    dr/dt = v(r(t), t),    r(0) = r₀.

Because ρ → ⟨ρ⟩ everywhere as t → ∞, the velocity decays to zero and
trajectories converge. In practice we integrate up to
``solver.convergence_time(tol)``, by which time residual motion is below
the tolerance.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .diffusion import DiffusionSolver


def _bilinear_sample(
    field: np.ndarray,
    xs_grid: np.ndarray,
    ys_grid: np.ndarray,
    xq: np.ndarray,
    yq: np.ndarray,
) -> np.ndarray:
    """Bilinear interpolation of ``field[i, j] = f(xs_grid[j], ys_grid[i])``.

    Queries outside the grid are clamped to the boundary. The grid is
    assumed to be uniformly spaced, which is guaranteed by the solver.
    """
    ny, nx = field.shape
    dx = xs_grid[1] - xs_grid[0]
    dy = ys_grid[1] - ys_grid[0]

    fx = (xq - xs_grid[0]) / dx
    fy = (yq - ys_grid[0]) / dy

    fx = np.clip(fx, 0.0, nx - 1.0 - 1e-12)
    fy = np.clip(fy, 0.0, ny - 1.0 - 1e-12)

    j0 = np.floor(fx).astype(np.intp)
    i0 = np.floor(fy).astype(np.intp)
    j1 = j0 + 1
    i1 = i0 + 1
    # Safe clamp — j1/i1 can exceed the last valid index after the epsilon
    # trick above for points that lie exactly at the edge.
    j1 = np.minimum(j1, nx - 1)
    i1 = np.minimum(i1, ny - 1)

    tx = fx - j0
    ty = fy - i0

    f00 = field[i0, j0]
    f01 = field[i0, j1]
    f10 = field[i1, j0]
    f11 = field[i1, j1]

    return (
        (1 - ty) * ((1 - tx) * f00 + tx * f01)
        + ty * ((1 - tx) * f10 + tx * f11)
    )


def advect_points(
    solver: DiffusionSolver,
    points: np.ndarray,
    t_max: float | None = None,
    rtol: float = 1e-6,
    atol: float = 1e-9,
    max_step: float | None = None,
) -> np.ndarray:
    """Flow ``points`` through the diffusion-induced velocity field.

    Parameters
    ----------
    solver : DiffusionSolver
        Precomputed diffusion on the density grid.
    points : np.ndarray, shape (N, 2)
        Initial ``(x, y)`` positions in physical coordinates.
    t_max : float, optional
        Upper integration time. Defaults to ``solver.convergence_time()``.
    rtol, atol : float
        Tolerances passed to ``scipy.integrate.solve_ivp``.
    max_step : float, optional
        Maximum RK step. Defaults to ``t_max / 64``, which keeps the time
        resolution comparable to the grid’s finest resolvable mode.

    Returns
    -------
    np.ndarray, shape (N, 2)
        Final positions after advection.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")

    if t_max is None:
        t_max = solver.convergence_time()

    if max_step is None:
        max_step = t_max / 64.0

    xs_grid, ys_grid = solver.grid_coords()
    n = points.shape[0]

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        rho, drho_dx, drho_dy = solver.density_and_gradient_at(t)
        xq = y[:n]
        yq = y[n:]
        rho_q = _bilinear_sample(rho, xs_grid, ys_grid, xq, yq)
        gx_q = _bilinear_sample(drho_dx, xs_grid, ys_grid, xq, yq)
        gy_q = _bilinear_sample(drho_dy, xs_grid, ys_grid, xq, yq)
        # Guard: ρ stays positive for all finite t because the DCT-II
        # expansion is analytic, but round-off near the final uniform
        # state can produce tiny values. Clamp defensively.
        rho_q = np.maximum(rho_q, 1e-300)
        vx = -gx_q / rho_q
        vy = -gy_q / rho_q
        return np.concatenate([vx, vy])

    y0 = np.concatenate([points[:, 0], points[:, 1]])
    sol = solve_ivp(
        rhs,
        (0.0, t_max),
        y0,
        method="RK45",
        rtol=rtol,
        atol=atol,
        max_step=max_step,
        dense_output=False,
        vectorized=False,
    )
    if not sol.success:
        raise RuntimeError(f"advection integration failed: {sol.message}")

    final = sol.y[:, -1]
    out = np.empty_like(points)
    out[:, 0] = final[:n]
    out[:, 1] = final[n:]
    return out


def clamp_to_bbox(points: np.ndarray, bbox: Tuple[float, float, float, float]) -> np.ndarray:
    """Clamp points to lie inside ``bbox``; useful after advection round-off."""
    xmin, ymin, xmax, ymax = bbox
    out = points.copy()
    out[:, 0] = np.clip(out[:, 0], xmin, xmax)
    out[:, 1] = np.clip(out[:, 1], ymin, ymax)
    return out
