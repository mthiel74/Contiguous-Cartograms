"""Tests for the world_cartogram one-liner.

Network access is forbidden, so these tests monkeypatch
``WorldMap.load`` with a synthetic two-country world constructed in
memory. Skipped when geopandas is unavailable.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest


def _fake_world(gpd, Polygon):
    """Two disjoint square "countries" in a 2:1 box."""
    gdf = gpd.GeoDataFrame(
        {
            "ISO_A3": ["AAA", "BBB"],
            "SOV_A3": ["AAA", "BBB"],
            "NAME": ["Aland", "Bland"],
            "population": [1.0, 9.0],
            "gdp": [1.0, 9.0],
            "gdp_per_capita": [1.0, 1.0],
            "geometry": [
                Polygon([(0.05, 0.05), (0.45, 0.05),
                         (0.45, 0.95), (0.05, 0.95)]),
                Polygon([(0.55, 0.05), (0.95, 0.05),
                         (0.95, 0.95), (0.55, 0.95)]),
            ],
        },
        crs="EPSG:4326",
    )
    return gdf


@pytest.fixture
def fake_world(monkeypatch):
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import Polygon
    from cartogram.world_data import WorldMap

    gdf = _fake_world(gpd, Polygon)

    def _load(cls, *args, **kwargs):
        return WorldMap(
            gdf=gdf.copy(),
            name_column="NAME",
            iso_column="ISO_A3",
        )

    monkeypatch.setattr(WorldMap, "load", classmethod(_load))
    return gdf


def test_world_cartogram_column_metric(fake_world):
    from cartogram import world_cartogram

    result = world_cartogram(
        "population",
        grid_size=(64, 128),
        render=False,
        performance_goal="quality",
    )
    assert result.figure is None
    assert len(result.deformed_geometries) == 2

    # The dense country (B, value=9) should end up larger than the
    # light country (A, value=1) on the cartogram.
    areas = [g.area for g in result.deformed_geometries]
    assert areas[1] > areas[0], (
        "high-value country did not expand on cartogram: "
        f"A={areas[0]:.3g}, B={areas[1]:.3g}"
    )


def test_world_cartogram_dict_metric(fake_world):
    from cartogram import world_cartogram

    result = world_cartogram(
        {"AAA": 1.0, "BBB": 9.0},
        label="demo",
        grid_size=(64, 128),
        render=False,
        performance_goal="quality",
    )
    assert result.world_map.value_column == "value"
    vals = result.world_map.gdf["value"].tolist()
    assert vals == [1.0, 9.0]


def test_world_cartogram_dict_metric_mixed_aliases(monkeypatch):
    """Dict keys in arbitrary forms (Wolfram-style, ISO2, common name)
    should route correctly via country_lookup."""
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import Polygon
    from cartogram.world_data import WorldMap

    gdf = gpd.GeoDataFrame(
        {
            "ISO_A3": ["USA", "GBR"],
            "SOV_A3": ["USA", "GBR"],
            "NAME": ["United States", "United Kingdom"],
            "population": [332.0, 67.0],
            "gdp": [25e12, 3e12],
            "gdp_per_capita": [75000.0, 45000.0],
            "geometry": [
                Polygon([(0.05, 0.05), (0.45, 0.05),
                         (0.45, 0.95), (0.05, 0.95)]),
                Polygon([(0.55, 0.05), (0.95, 0.05),
                         (0.95, 0.95), (0.55, 0.95)]),
            ],
        },
        crs="EPSG:4326",
    )

    def _load(cls, *args, **kwargs):
        return WorldMap(
            gdf=gdf.copy(),
            name_column="NAME",
            iso_column="ISO_A3",
        )

    monkeypatch.setattr(WorldMap, "load", classmethod(_load))

    from cartogram import world_cartogram

    result = world_cartogram(
        {
            "UnitedStates": 10.0,     # Wolfram-style
            "UK": 1.0,                # alpha-2
        },
        label="aliases",
        grid_size=(64, 128),
        render=False,
        performance_goal="quality",
    )
    vals = result.world_map.gdf.set_index("ISO_A3")["value"].to_dict()
    assert vals == {"USA": 10.0, "GBR": 1.0}


def test_world_cartogram_callable_metric(fake_world):
    from cartogram import world_cartogram

    result = world_cartogram(
        lambda row: float(row["population"]) ** 2,
        label="pop squared",
        grid_size=(64, 128),
        render=False,
        performance_goal="quality",
    )
    assert result.world_map.gdf["value"].tolist() == [1.0, 81.0]


def test_world_cartogram_unknown_column(fake_world):
    from cartogram import world_cartogram

    with pytest.raises(KeyError, match="not_a_column"):
        world_cartogram("not_a_column", grid_size=(32, 64), render=False)
