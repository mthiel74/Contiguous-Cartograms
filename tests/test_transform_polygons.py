"""Tests for the vectorised polygon-transform path."""

from __future__ import annotations

import numpy as np
import pytest


def test_transform_polygons_preserves_count_and_topology():
    """MultiPolygon with a hole round-trips through the vectorised
    transform with the same structure: same number of sub-polygons,
    same count of interior rings, same (lat, lon) orientation."""
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import MultiPolygon, Polygon

    from cartogram import Cartogram

    rho = np.ones((32, 32))
    rho[10:22, 10:22] = 5.0
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3, performance_goal="quality")

    outer_with_hole = Polygon(
        [(0.05, 0.05), (0.45, 0.05), (0.45, 0.45), (0.05, 0.45)],
        holes=[[(0.10, 0.10), (0.20, 0.10), (0.20, 0.20), (0.10, 0.20)]],
    )
    second_part = Polygon(
        [(0.55, 0.55), (0.95, 0.55), (0.95, 0.95), (0.55, 0.95)]
    )
    mp = MultiPolygon([outer_with_hole, second_part])

    warped = cart.transform_polygons([mp, second_part])
    assert len(warped) == 2

    mp_warped = warped[0]
    assert isinstance(mp_warped, MultiPolygon)
    # Same number of sub-polygons.
    assert len(mp_warped.geoms) == 2
    # Hole preserved on the first sub-polygon.
    with_hole = mp_warped.geoms[0]
    assert len(list(with_hole.interiors)) == 1

    # CCW-outer orientation preserved by the explicit orient() call.
    assert with_hole.exterior.is_ccw
    # Interior ring should be CW for a well-oriented polygon (shapely
    # orient() enforces this).
    assert not list(with_hole.interiors)[0].is_ccw


def test_transform_polygons_matches_reference_loop():
    """The vectorised implementation must agree with a per-ring
    reference to bilinear-sample precision."""
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import Polygon

    from cartogram import Cartogram

    rho = np.ones((24, 24))
    rho[8:16, 8:16] = 4.0
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3, performance_goal="quality")

    triangle = Polygon([(0.1, 0.1), (0.8, 0.2), (0.3, 0.7)])
    warped_vec = cart.transform_polygons([triangle])[0]

    # Reference: run self.transform on the exterior coords directly.
    coords = np.asarray(list(triangle.exterior.coords), dtype=float)
    # Match the vectorised path: segmentize to one-grid-cell edges.
    from shapely import segmentize
    dense = segmentize(triangle, min(cart.solver.dx, cart.solver.dy))
    dense_coords = np.asarray(list(dense.exterior.coords), dtype=float)
    moved_ref = cart.transform(dense_coords)

    warped_coords = np.asarray(list(warped_vec.exterior.coords), dtype=float)
    assert warped_coords.shape == moved_ref.shape
    # Shapely may reverse orientation after orient(); compare as sets.
    a = {tuple(p) for p in warped_coords}
    b = {tuple(p) for p in moved_ref}
    # Floating-point-identical for the same sample points.
    assert a == b


def test_transform_polygons_requires_run():
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import Polygon

    from cartogram import Cartogram

    rho = np.ones((8, 8))
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    tri = Polygon([(0.1, 0.1), (0.4, 0.1), (0.2, 0.4)])
    with pytest.raises(RuntimeError, match="run"):
        cart.transform_polygons([tri])


def test_transform_polygons_rejects_non_polygon_types():
    shapely = pytest.importorskip("shapely")
    from shapely.geometry import Point

    from cartogram import Cartogram

    rho = np.ones((8, 8))
    cart = Cartogram(rho, bbox=(0, 0, 1.0, 1.0))
    cart.run(tol=1e-3)
    with pytest.raises(TypeError, match="expected Polygon"):
        cart.transform_polygons([Point(0.1, 0.1)])
