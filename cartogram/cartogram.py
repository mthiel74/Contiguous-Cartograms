"""High-level :class:`Cartogram` wrapper tying diffusion to point advection."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, List, Sequence, Tuple

import numpy as np

from .advect import advect_points, clamp_to_bbox
from .diffusion import DiffusionSolver

BBox = Tuple[float, float, float, float]


@dataclass
class Cartogram:
    """Driver that turns a density grid into a contiguous cartogram.

    The Gastner–Newman prescription is:

    1. Add a uniform background of ``mean_floor · ⟨ρ₀⟩`` so that density
       is strictly positive everywhere, which keeps the velocity field
       ``v = -∇ρ / ρ`` finite. The paper calls this the "sea" trick.
    2. Solve the heat equation for ``ρ(·, t)``.
    3. Flow a collection of points through ``v``; the limiting positions
       are the cartogram coordinates.

    Parameters
    ----------
    rho : np.ndarray, shape (ny, nx)
        Raw density grid (for example a population raster). Non-negative.
    bbox : (xmin, ymin, xmax, ymax)
        Physical extent of the grid.
    mean_floor : float
        Background fraction added to every cell, expressed as a multiple
        of the mean of the input density. ``0.0`` would produce infinite
        velocities in empty cells. Default ``0.005`` follows the paper.
    blur_sigma : float
        Standard deviation (in grid cells) of a Gaussian pre-smoothing
        applied to the density. Mathematically equivalent to starting
        the heat equation at t = sigma² dx dy / 2; practically, it
        regularises the near-step discontinuities produced by polygon
        rasterisation so the advection ODE is non-stiff at t = 0.
        Default ``1.0``.
    sea_density : float | "auto" | None
        If set, cells with zero input density (i.e. the ocean) are
        filled with this value before diffusion. ``"auto"`` fills them
        with the mean density of the non-zero (land) cells — which
        approximately preserves ocean area, because then the spatial
        mean of `ρ` is the same on land and sea and the deformation
        only redistributes material *within* the land. A custom float
        interpolates between the two regimes: higher values compress
        land more, lower values expand land more.
        Default ``None`` (no ocean fill; legacy behaviour).
    """

    rho: np.ndarray
    bbox: BBox
    mean_floor: float = 0.005
    blur_sigma: float = 1.0
    sea_density: float | str | None = None

    solver: DiffusionSolver = field(init=False, repr=False)
    _final_grid_points: np.ndarray | None = field(default=None, init=False, repr=False)

    def __post_init__(self) -> None:
        rho = np.asarray(self.rho, dtype=float)
        if rho.ndim != 2:
            raise ValueError("rho must be a 2-D array (ny, nx)")
        if np.any(rho < 0):
            raise ValueError("rho must be non-negative")
        if self.mean_floor <= 0:
            raise ValueError("mean_floor must be strictly positive")

        if self.sea_density is not None:
            land_mask = rho > 0
            if not land_mask.any():
                raise ValueError("input density is identically zero")
            if isinstance(self.sea_density, str):
                if self.sea_density != "auto":
                    raise ValueError(
                        "sea_density must be a float, 'auto', or None"
                    )
                ocean_val = float(rho[land_mask].mean())
            else:
                ocean_val = float(self.sea_density)
            rho_prepared = np.where(land_mask, rho, ocean_val)
            # A tiny mean_floor is still useful to guarantee positivity
            # for solver robustness, but it can be much smaller here
            # because the ocean is no longer near-zero.
            rho_prepared = rho_prepared + self.mean_floor * float(rho_prepared.mean())
        else:
            floor = self.mean_floor * float(rho.mean())
            if floor == 0.0:
                raise ValueError("input density is identically zero")
            rho_prepared = rho + floor

        if self.blur_sigma > 0:
            from scipy.ndimage import gaussian_filter
            rho_prepared = gaussian_filter(
                rho_prepared, sigma=self.blur_sigma, mode="reflect"
            )

        self.rho = rho_prepared
        self.solver = DiffusionSolver(rho_prepared, self.bbox)

    # ------------------------------------------------------------------
    def run(
        self,
        tol: float = 1e-3,
        rtol: float = 1e-5,
        atol: float = 1e-7,
        method: str = "RK45",
    ) -> None:
        """Precompute the deformation of the grid itself.

        After :meth:`run`, :meth:`transform` becomes essentially free: we
        bilinearly interpolate the displacement field defined by the
        already-advected grid points, rather than re-integrating the ODE
        for every query.
        """
        xs, ys = self.solver.grid_coords()
        gx, gy = np.meshgrid(xs, ys)
        pts = np.column_stack([gx.ravel(), gy.ravel()])
        t_max = self.solver.convergence_time(tol=tol)
        moved = advect_points(
            self.solver,
            pts,
            t_max=t_max,
            rtol=rtol,
            atol=atol,
            method=method,
        )
        moved = clamp_to_bbox(moved, self.bbox)
        self._final_grid_points = moved

    # ------------------------------------------------------------------
    def inverse_transform(self, points: np.ndarray) -> np.ndarray:
        """Inverse of :meth:`transform`. Defers to :mod:`cartogram.warp`."""
        from .warp import _inverse_transform
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("points must have shape (N, 2)")
        return _inverse_transform(self, points)

    def transform(self, points: np.ndarray) -> np.ndarray:
        """Map ``points`` from original to cartogram coordinates.

        Requires a prior call to :meth:`run`.
        """
        if self._final_grid_points is None:
            raise RuntimeError("Cartogram.run() must be called before transform()")
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2:
            raise ValueError("points must have shape (N, 2)")

        ny, nx = self.solver.shape
        xs, ys = self.solver.grid_coords()
        disp_x = self._final_grid_points[:, 0].reshape(ny, nx)
        disp_y = self._final_grid_points[:, 1].reshape(ny, nx)

        from .advect import _bilinear_sample

        out = np.empty_like(points)
        out[:, 0] = _bilinear_sample(disp_x, xs, ys, points[:, 0], points[:, 1])
        out[:, 1] = _bilinear_sample(disp_y, xs, ys, points[:, 0], points[:, 1])
        return out

    # ------------------------------------------------------------------
    def transform_polygons(
        self,
        polygons: Sequence,
        max_edge: float | None = None,
    ) -> List:
        """Apply :meth:`transform` to a sequence of shapely polygons.

        Long boundary edges are first densified so that, after
        deformation, the polygon boundary follows the cartogram flow
        smoothly rather than jumping across high-gradient regions in a
        single straight step. ``max_edge`` sets the densification
        threshold in physical units; the default is one grid cell.

        Because the deformation is a homeomorphism, the topology of each
        polygon is preserved (no self-intersection on well-resolved
        inputs). Holes are preserved.
        """
        from shapely.geometry import Polygon, MultiPolygon
        from shapely.geometry.polygon import orient
        from shapely import segmentize

        if max_edge is None:
            max_edge = min(self.solver.dx, self.solver.dy)

        def _transform_ring(coords: Iterable[Tuple[float, float]]) -> list:
            arr = np.asarray(list(coords), dtype=float)
            if arr.ndim != 2 or arr.shape[1] < 2:
                return []
            arr = arr[:, :2]
            moved = self.transform(arr)
            return [tuple(p) for p in moved]

        def _transform_one(poly: Polygon) -> Polygon:
            shell = _transform_ring(poly.exterior.coords)
            holes = [_transform_ring(h.coords) for h in poly.interiors]
            return orient(Polygon(shell, holes))

        out = []
        for geom in polygons:
            dense = segmentize(geom, max_edge)
            if isinstance(dense, MultiPolygon):
                out.append(MultiPolygon([_transform_one(p) for p in dense.geoms]))
            elif isinstance(dense, Polygon):
                out.append(_transform_one(dense))
            else:
                raise TypeError(
                    f"expected Polygon or MultiPolygon, got {type(geom).__name__}"
                )
        return out
