"""Analytic diffusion on a rectangle with Neumann (no-flux) boundary conditions.

The 2-D heat equation

    dρ/dt = ∇²ρ,       x ∈ [0, Lx],  y ∈ [0, Ly]
    (∂ρ/∂n)|boundary = 0

is solved by cosine expansion. On a cell-centred grid of shape (ny, nx) the
DCT-II basis functions are

    φ_{mn}(x, y) = cos(m π x / Lx) · cos(n π y / Ly),    m, n ≥ 0

with eigenvalues

    λ_{mn} = π² (m² / Lx² + n² / Ly²).

The density at time t is therefore

    ρ(x, y, t) = Σ_{m, n} A_{mn} · exp(-λ_{mn} · t) · φ_{mn}(x, y).

This file exposes :class:`DiffusionSolver` which precomputes the coefficients
``A_{mn}`` and the eigenvalue grid ``λ_{mn}`` once, then provides O(N log N)
evaluation of ρ(t) and of its gradient for arbitrary t.

Grid convention
---------------
All density arrays have shape ``(ny, nx)``. Index ``(i, j)`` corresponds to
physical coordinate

    x_j = xmin + (j + 0.5) · dx,     dx = (xmax - xmin) / nx,
    y_i = ymin + (i + 0.5) · dy,     dy = (ymax - ymin) / ny.

(Cell-centred placement is what DCT-II expects.)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.fft import dctn, idctn


BBox = Tuple[float, float, float, float]


@dataclass
class DiffusionSolver:
    """Analytic heat-equation solver on a cell-centred rectangular grid.

    Parameters
    ----------
    rho0 : np.ndarray, shape (ny, nx)
        Initial density. Must be strictly positive after the optional offset
        (see :meth:`Cartogram` for the standard mean-floor trick).
    bbox : (xmin, ymin, xmax, ymax)
        Physical extent of the grid.
    """

    rho0: np.ndarray
    bbox: BBox

    def __post_init__(self) -> None:
        if self.rho0.ndim != 2:
            raise ValueError("rho0 must be a 2-D array of shape (ny, nx)")
        if np.any(self.rho0 <= 0):
            raise ValueError(
                "rho0 must be strictly positive; add a background floor "
                "before constructing the solver."
            )

        xmin, ymin, xmax, ymax = self.bbox
        if not (xmax > xmin and ymax > ymin):
            raise ValueError(f"invalid bbox {self.bbox}")

        ny, nx = self.rho0.shape
        self._ny, self._nx = ny, nx
        self._Lx = xmax - xmin
        self._Ly = ymax - ymin
        self._dx = self._Lx / nx
        self._dy = self._Ly / ny

        # Cosine coefficients (orthonormal DCT-II, so forward/inverse agree).
        self._coeffs = dctn(self.rho0, type=2, norm="ortho")

        # Eigenvalues of -Δ on the Neumann box. Row index -> y-mode n,
        # column index -> x-mode m.
        m = np.arange(nx)
        n = np.arange(ny)
        lam_x = (np.pi * m / self._Lx) ** 2
        lam_y = (np.pi * n / self._Ly) ** 2
        self._lam = lam_y[:, None] + lam_x[None, :]

        self._rho_mean = float(self.rho0.mean())

    # -- accessors --------------------------------------------------------
    @property
    def shape(self) -> Tuple[int, int]:
        return self._ny, self._nx

    @property
    def dx(self) -> float:
        return self._dx

    @property
    def dy(self) -> float:
        return self._dy

    @property
    def mean_density(self) -> float:
        return self._rho_mean

    # -- main API ---------------------------------------------------------
    def density_at(self, t: float) -> np.ndarray:
        """Return ρ(·, t) sampled on the original grid."""
        decayed = self._coeffs * np.exp(-self._lam * float(t))
        return idctn(decayed, type=2, norm="ortho")

    def density_and_gradient_at(
        self, t: float
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(ρ, ∂ρ/∂x, ∂ρ/∂y)`` on the grid at time ``t``.

        Gradients are computed by second-order central differences with
        one-sided stencils at the boundary (consistent with the Neumann BC,
        which sends the normal derivative to zero as t → ∞ anyway).
        """
        rho = self.density_at(t)
        drho_dy, drho_dx = np.gradient(rho, self._dy, self._dx)
        return rho, drho_dx, drho_dy

    # -- convenience helpers ----------------------------------------------
    def grid_coords(self) -> Tuple[np.ndarray, np.ndarray]:
        """Return the 1-D arrays of x and y grid-point coordinates."""
        xmin, ymin, _, _ = self.bbox
        xs = xmin + (np.arange(self._nx) + 0.5) * self._dx
        ys = ymin + (np.arange(self._ny) + 0.5) * self._dy
        return xs, ys

    def convergence_time(self, tol: float = 1e-3) -> float:
        """Time at which the slowest non-zero mode has decayed to ``tol``.

        This provides a sensible upper bound for the advection integration:
        running longer is wasted work because ρ is already uniform to
        within ``tol`` and the induced velocity field is essentially zero.
        """
        lam_min = min(
            (np.pi / self._Lx) ** 2,
            (np.pi / self._Ly) ** 2,
        )
        return float(-np.log(tol) / lam_min)
