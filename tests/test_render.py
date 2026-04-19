"""Tests for the side_by_side renderer."""

from __future__ import annotations

import numpy as np
import pytest


def _fake_gdf_and_deformed(gpd, Polygon):
    gdf = gpd.GeoDataFrame(
        {
            "value": [1.0, 9.0],
            "geometry": [
                Polygon([(0.05, 0.05), (0.45, 0.05),
                         (0.45, 0.95), (0.05, 0.95)]),
                Polygon([(0.55, 0.05), (0.95, 0.05),
                         (0.95, 0.95), (0.55, 0.95)]),
            ],
        },
        crs="EPSG:4326",
    )
    # Trivial "deformed" geometries: identical to geographic.
    deformed = list(gdf.geometry)
    return gdf, deformed


def test_side_by_side_returns_figure_with_two_axes():
    gpd = pytest.importorskip("geopandas")
    _ = pytest.importorskip("matplotlib")
    from shapely.geometry import Polygon
    from cartogram import side_by_side

    gdf, deformed = _fake_gdf_and_deformed(gpd, Polygon)
    fig = side_by_side(gdf, deformed, label="demo")
    assert len(fig.axes) == 2
    titles = [ax.get_title() for ax in fig.axes]
    assert "Geographic" in titles[0]
    assert "Cartogram" in titles[1]


def test_side_by_side_mask_shape_check():
    gpd = pytest.importorskip("geopandas")
    _ = pytest.importorskip("matplotlib")
    from shapely.geometry import Polygon
    from cartogram import side_by_side

    gdf, deformed = _fake_gdf_and_deformed(gpd, Polygon)
    with pytest.raises(ValueError, match="mask"):
        side_by_side(
            gdf, deformed,
            has_value_mask=np.array([True, False, True]),  # wrong length
        )


def test_side_by_side_grey_path_draws_both_sets():
    gpd = pytest.importorskip("geopandas")
    _ = pytest.importorskip("matplotlib")
    from shapely.geometry import Polygon
    from cartogram import side_by_side

    gdf, deformed = _fake_gdf_and_deformed(gpd, Polygon)
    mask = np.array([True, False])           # second country "missing"
    fig = side_by_side(gdf, deformed, has_value_mask=mask)
    # Each axis should carry at least two PolyCollection-like child
    # artists: one for the grey set, one for the valued set.
    for ax in fig.axes:
        collections = [c for c in ax.collections]
        assert len(collections) >= 2
