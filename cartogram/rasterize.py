"""Convert a list of polygons + per-polygon values into a density grid.

Density = value / polygon_area, distributed uniformly over the polygon. We
use a subpixel-sampling rasteriser (default 3×3 samples per cell) to
mitigate the aliasing that would otherwise arise at region boundaries.

Dependencies are limited to ``numpy`` and ``shapely``; no GIS stack
required.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np

BBox = Tuple[float, float, float, float]


def rasterize_polygons(
    polygons: Sequence,
    values: Sequence[float],
    bbox: BBox,
    shape: Tuple[int, int],
    subpixel: int = 3,
) -> np.ndarray:
    """Rasterise ``(polygon, value)`` pairs onto a grid of shape ``(ny, nx)``.

    Each polygon contributes a density of ``value / area`` to the cells it
    covers, assessed by sampling ``subpixel × subpixel`` points per cell
    and counting how many fall inside the polygon (assigning a fractional
    cell coverage).

    Parameters
    ----------
    polygons : sequence of shapely (Multi)Polygon
    values   : sequence of numbers, same length as ``polygons``
    bbox     : (xmin, ymin, xmax, ymax)
    shape    : (ny, nx) of the target grid
    subpixel : integer oversampling factor per cell (≥1)

    Returns
    -------
    np.ndarray, shape ``(ny, nx)``
        Density grid in units of ``value / area``.
    """
    from shapely.geometry import Polygon, MultiPolygon
    from shapely.prepared import prep

    if len(polygons) != len(values):
        raise ValueError("polygons and values must have equal length")
    if subpixel < 1:
        raise ValueError("subpixel must be >= 1")

    ny, nx = shape
    xmin, ymin, xmax, ymax = bbox
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny

    # Subpixel offsets within a single cell, evenly spaced in (0, 1).
    offs = (np.arange(subpixel) + 0.5) / subpixel
    sx = offs * dx
    sy = offs * dy

    density = np.zeros((ny, nx), dtype=float)

    for poly, val in zip(polygons, values):
        if val == 0:
            continue
        if poly.is_empty:
            continue
        area = poly.area
        if area <= 0:
            continue

        # Only iterate cells overlapping the polygon's bounding box.
        pminx, pminy, pmaxx, pmaxy = poly.bounds
        j0 = max(0, int(np.floor((pminx - xmin) / dx)))
        j1 = min(nx, int(np.ceil((pmaxx - xmin) / dx)))
        i0 = max(0, int(np.floor((pminy - ymin) / dy)))
        i1 = min(ny, int(np.ceil((pmaxy - ymin) / dy)))
        if j0 >= j1 or i0 >= i1:
            continue

        prepared = prep(poly)
        cell_density = val / area
        inv_samples = 1.0 / (subpixel * subpixel)

        from shapely.geometry import Point

        for i in range(i0, i1):
            cy = ymin + i * dy + sy  # subpixel y-coords
            for j in range(j0, j1):
                cx = xmin + j * dx + sx  # subpixel x-coords
                hits = 0
                for yy in cy:
                    for xx in cx:
                        if prepared.contains(Point(xx, yy)):
                            hits += 1
                if hits:
                    density[i, j] += cell_density * hits * inv_samples

    return density


def build_grid_points(bbox: BBox, shape: Tuple[int, int]) -> np.ndarray:
    """Return a ``(ny*nx, 2)`` array of cell-centred grid-point coordinates."""
    ny, nx = shape
    xmin, ymin, xmax, ymax = bbox
    xs = xmin + (np.arange(nx) + 0.5) * (xmax - xmin) / nx
    ys = ymin + (np.arange(ny) + 0.5) * (ymax - ymin) / ny
    gx, gy = np.meshgrid(xs, ys)
    return np.column_stack([gx.ravel(), gy.ravel()])
