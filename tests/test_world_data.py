"""Tests for the WorldMap loader + helpers.

Most tests are skipped when geopandas / pandas are not installed, since
the full world-map path is an optional extra. We still test the pure
CSV-parsing path (``csv_to_value_map``) unconditionally.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from cartogram.world_data import csv_to_value_map


def test_csv_to_value_map(tmp_path: Path) -> None:
    csv = tmp_path / "vals.csv"
    csv.write_text(
        'iso,x\n'
        'USA,1.5\n'
        'DEU,"12,345"\n'            # quoted thousands separator
        'FRA,\n'                    # missing -> skipped
        ',9.9\n'                    # missing iso -> skipped
        'JPN,not-a-number\n'        # bad number -> skipped
        'GBR,2.5\n',
        encoding="utf-8",
    )
    mapping = csv_to_value_map(str(csv), "iso", "x")
    assert mapping == {"USA": 1.5, "DEU": 12345.0, "GBR": 2.5}


def test_merge_by_sovereign_sums_territories() -> None:
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import Polygon

    from cartogram.world_data import WorldMap

    gdf = gpd.GeoDataFrame(
        {
            "SOV_A3": ["DN1", "DN1", "FRA"],
            "NAME": ["Denmark", "Greenland", "France"],
            "population": [5_900_000.0, 56_000.0, 67_000_000.0],
            "gdp": [400e9, 3e9, 3_000e9],
            "geometry": [
                Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
                Polygon([(2, 2), (3, 2), (3, 3), (2, 3)]),
                Polygon([(5, 0), (6, 0), (6, 1), (5, 1)]),
            ],
        },
        crs="EPSG:4326",
    )
    # Seed a gdp_per_capita column so we can verify the post-hoc recompute.
    gdf["gdp_per_capita"] = gdf["gdp"] / gdf["population"]
    wm = WorldMap(gdf=gdf, name_column="NAME", iso_column="SOV_A3")

    merged = wm.merge_by_sovereign(group_column="SOV_A3")

    assert len(merged.gdf) == 2
    dn = merged.gdf[merged.gdf["SOV_A3"] == "DN1"].iloc[0]
    assert dn["population"] == pytest.approx(5_956_000.0)
    assert dn["gdp"] == pytest.approx(403e9)
    # gdp_per_capita should be aggregated_gdp / aggregated_population,
    # NOT an average of the per-row gdp_per_capita values.
    expected_pc = 403e9 / 5_956_000.0
    assert dn["gdp_per_capita"] == pytest.approx(expected_pc, rel=1e-9)
    # The dissolved Denmark geometry must be a MultiPolygon covering
    # both the mainland box and the Greenland box.
    assert dn.geometry.geom_type == "MultiPolygon"
    assert dn.geometry.area == pytest.approx(2.0)


def test_merge_by_sovereign_unknown_column() -> None:
    gpd = pytest.importorskip("geopandas")
    from shapely.geometry import Polygon

    from cartogram.world_data import WorldMap

    gdf = gpd.GeoDataFrame(
        {"NAME": ["X"], "geometry": [Polygon([(0, 0), (1, 0), (1, 1)])]},
        crs="EPSG:4326",
    )
    wm = WorldMap(gdf=gdf, name_column="NAME")
    with pytest.raises(KeyError, match="SOV_A3"):
        wm.merge_by_sovereign()
