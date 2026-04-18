"""Optional I/O helpers. These are *only* imported when called, so the
core algorithm stays free of a GeoPandas dependency."""

from __future__ import annotations

from typing import Tuple


def read_geodataframe(path: str):
    """Load a shapefile / GeoPackage / GeoJSON into a GeoDataFrame.

    Requires ``geopandas``. Import is deferred so the core package stays
    importable without it.
    """
    import geopandas as gpd

    return gpd.read_file(path)


def polygons_and_values(gdf, value_column: str) -> Tuple[list, list]:
    """Extract a list of shapely geometries and a parallel list of numeric
    values from a GeoDataFrame column."""
    polys = list(gdf.geometry)
    values = [float(v) for v in gdf[value_column].to_numpy()]
    return polys, values
