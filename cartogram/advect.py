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

Performance goal
----------------
A naive implementation recomputes the full inverse DCT of ρ(·, t) on
every RK45 sub-step — hundreds of inverse DCTs per Cartogram.run() on
a 256×256 grid. The ``performance_goal="speed"`` path (default) caches
``(ρ, ∂xρ, ∂yρ)`` at a small number of geometrically-spaced time
samples and linearly interpolates in t inside the RHS, matching the
Gastner–Newman exact-decay solution up to a few parts in 10⁻³ (the
exponential modes are smooth in t, so linear interpolation between
nearby samples is very accurate). ``performance_goal="quality"``
preserves the original behaviour for bit-exact reproduction.
"""

from __future__ import annotations

from typing import Literal, Tuple

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


class _SnapshotCache:
    """Geometrically-spaced time samples of ``(ρ, ∂xρ, ∂yρ)`` for fast RHS eval.

    The modes of the Gastner–Newman diffusion decay as ``exp(-λ_mn t)``,
    so the densities change fastest near ``t = 0`` and barely at all near
    ``t_max``. We place snapshots geometrically from just before the
    fastest mode collapses out to ``t_max`` (with ``t = 0`` prepended),
    matching the Wolfram port's ``PerformanceGoal -> "Speed"`` scheduling.
    """

    def __init__(self, solver: DiffusionSolver, t_max: float, n_snapshots: int):
        if n_snapshots < 2:
            raise ValueError("n_snapshots must be >= 2")

        lam_max = float(np.max(solver._lam))
        if lam_max <= 0:
            t0 = t_max / 1000.0
        else:
            t0 = min(0.1 / lam_max, t_max / 1000.0)
        ratio = (t_max / t0) ** (1.0 / (n_snapshots - 1))
        times = np.empty(n_snapshots + 1, dtype=float)
        times[0] = 0.0
        times[1:] = t0 * ratio ** np.arange(n_snapshots)

        rhos = np.empty((len(times), *solver.shape), dtype=float)
        gxs = np.empty_like(rhos)
        gys = np.empty_like(rhos)
        for k, t in enumerate(times):
            rho, gx, gy = solver.density_and_gradient_at(float(t))
            rhos[k] = rho
            gxs[k] = gx
            gys[k] = gy

        self.times = times
        self.rhos = rhos
        self.gxs = gxs
        self.gys = gys

    def sample(self, t: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(ρ, ∂xρ, ∂yρ)`` at ``t`` by linear interpolation."""
        t = float(t)
        times = self.times
        # searchsorted finds the right bracket.
        k = int(np.searchsorted(times, t, side="right")) - 1
        k = max(0, min(k, len(times) - 2))
        t_lo, t_hi = times[k], times[k + 1]
        alpha = (t - t_lo) / (t_hi - t_lo)
        alpha = float(np.clip(alpha, 0.0, 1.0))
        rho = (1.0 - alpha) * self.rhos[k] + alpha * self.rhos[k + 1]
        gx = (1.0 - alpha) * self.gxs[k] + alpha * self.gxs[k + 1]
        gy = (1.0 - alpha) * self.gys[k] + alpha * self.gys[k + 1]
        return rho, gx, gy


def advect_points(
    solver: DiffusionSolver,
    points: np.ndarray,
    t_max: float | None = None,
    rtol: float = 1e-5,
    atol: float = 1e-7,
    max_step: float | None = None,
    method: str = "RK45",
    performance_goal: Literal["speed", "quality"] = "speed",
    snapshots: int = 60,
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
        resolution comparable to the grid's finest resolvable mode.
    method : str
        Integrator passed to ``scipy.integrate.solve_ivp``.
    performance_goal : {"speed", "quality"}
        ``"speed"`` (default) pre-computes ``snapshots`` geometrically-spaced
        density/gradient snapshots and linearly interpolates between them
        inside the RHS. Matches ``"quality"`` to sub-grid-cell accuracy on
        smooth densities and is roughly 5–7× faster on world-scale
        cartograms. ``"quality"`` re-evaluates the inverse DCT at every RK
        sub-step for bit-exact reproduction with the reference solver.
    snapshots : int
        Snapshot count used when ``performance_goal="speed"``. More is
        more accurate but diminishing returns past ~60.

    Returns
    -------
    np.ndarray, shape (N, 2)
        Final positions after advection.
    """
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")

    if performance_goal not in ("speed", "quality"):
        raise ValueError(
            f"performance_goal must be 'speed' or 'quality'; got {performance_goal!r}"
        )

    if t_max is None:
        t_max = solver.convergence_time()

    if max_step is None:
        max_step = t_max / 64.0

    xs_grid, ys_grid = solver.grid_coords()
    n = points.shape[0]

    cache: _SnapshotCache | None = None
    if performance_goal == "speed":
        cache = _SnapshotCache(solver, t_max, snapshots)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        if cache is not None:
            rho, drho_dx, drho_dy = cache.sample(t)
        else:
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
        method=method,
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
