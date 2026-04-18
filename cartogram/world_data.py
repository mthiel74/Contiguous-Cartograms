"""World-map data utilities.

Fetches the Natural Earth 1:110m admin-0 countries shapefile on demand
(cached under ``~/.cache/cartogram/``), reads it into a GeoDataFrame,
and offers a small CSV-join helper so users can bring their own values.

The Natural Earth dataset ships with ``POP_EST`` (population estimate)
and ``GDP_MD`` (GDP in millions USD) columns, which we re-expose under
the friendlier names ``population``, ``gdp`` and ``gdp_per_capita``.
"""

from __future__ import annotations

import os
import urllib.request
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence

NE_URL_110 = (
    "https://naciscdn.org/naturalearth/110m/cultural/"
    "ne_110m_admin_0_countries.zip"
)
NE_URL_50 = (
    "https://naciscdn.org/naturalearth/50m/cultural/"
    "ne_50m_admin_0_countries.zip"
)

CACHE_DIR = Path.home() / ".cache" / "cartogram"


def _cached_shapefile(resolution: str) -> Path:
    """Ensure the Natural Earth admin0 shapefile is present; return its path."""
    if resolution not in ("110m", "50m"):
        raise ValueError("resolution must be '110m' or '50m'")
    url = NE_URL_110 if resolution == "110m" else NE_URL_50
    shp_name = f"ne_{resolution}_admin_0_countries.shp"
    zip_path = CACHE_DIR / f"ne_{resolution}_admin_0_countries.zip"
    shp_path = CACHE_DIR / shp_name

    if shp_path.exists():
        return shp_path

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {url} -> {zip_path} ...", flush=True)
    urllib.request.urlretrieve(url, zip_path)
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(CACHE_DIR)
    if not shp_path.exists():
        raise RuntimeError(
            f"Expected {shp_path} in the Natural Earth zip but it was not found."
        )
    return shp_path


@dataclass
class WorldMap:
    """Wrapper around a GeoDataFrame of world countries.

    Attributes
    ----------
    gdf : geopandas.GeoDataFrame
        The loaded countries, in the given projection.
    value_column : str | None
        Name of the column whose values should drive the cartogram.
    name_column : str
        Human-readable country name (for plotting / debugging).
    iso_column : str
        ISO-3166-1 alpha-3 code column (for CSV joins).
    """

    gdf: "object"  # geopandas.GeoDataFrame
    value_column: Optional[str] = None
    name_column: str = "NAME"
    iso_column: str = "ISO_A3"

    # ------------------------------------------------------------------
    @classmethod
    def load(
        cls,
        resolution: str = "110m",
        projection: str = "EPSG:8857",
        drop_antarctica: bool = True,
    ) -> "WorldMap":
        """Load the Natural Earth admin-0 countries dataset.

        Parameters
        ----------
        resolution : '110m' or '50m'
            Source resolution. 110m (~170 kB shapefile) is fine for most
            cartograms; 50m (~5 MB) for higher-quality rendering.
        projection : str
            Any PROJ/EPSG string. Default EPSG:8857 is Equal Earth — a
            pleasing equal-area projection that keeps the world roughly
            rectangular, which is important because a cartogram is built
            on a rectangular grid.
        drop_antarctica : bool
            Antarctica's huge geographical area + negligible quantitative
            value dominates most cartograms; off by default.
        """
        import geopandas as gpd

        shp = _cached_shapefile(resolution)
        gdf = gpd.read_file(shp)
        if drop_antarctica:
            gdf = gdf[gdf["CONTINENT"] != "Antarctica"].copy()
        gdf = gdf.to_crs(projection)
        gdf = gdf.reset_index(drop=True)

        # Normalise common fields.
        if "POP_EST" in gdf.columns:
            gdf["population"] = gdf["POP_EST"].astype(float)
        if "GDP_MD" in gdf.columns:
            gdf["gdp"] = gdf["GDP_MD"].astype(float) * 1e6  # USD
        if "population" in gdf.columns and "gdp" in gdf.columns:
            pop = gdf["population"].replace(0, float("nan"))
            gdf["gdp_per_capita"] = gdf["gdp"] / pop
            gdf["gdp_per_capita"] = gdf["gdp_per_capita"].fillna(0.0)

        iso_col = "ISO_A3" if "ISO_A3" in gdf.columns else "ADM0_A3"
        name_col = "NAME" if "NAME" in gdf.columns else gdf.columns[0]
        return cls(gdf=gdf, name_column=name_col, iso_column=iso_col)

    # ------------------------------------------------------------------
    def with_values(
        self,
        values: dict | "object",
        join_key: str = "ISO_A3",
        value_key: str = "value",
        as_column: str = "value",
    ) -> "WorldMap":
        """Join an external value table onto the world map.

        Parameters
        ----------
        values : dict or pandas.DataFrame
            If a dict, keys are country codes (matching ``join_key``) and
            values are the numbers to visualise.
            If a DataFrame, it must contain ``join_key`` and ``value_key``
            columns.
        join_key : str
            Column in ``self.gdf`` to join on. Typically ``ISO_A3``.
        value_key : str
            When passing a DataFrame, the column name holding the values.
        as_column : str
            Name for the new column in the returned :class:`WorldMap`.
        """
        import pandas as pd

        gdf = self.gdf.copy()
        if isinstance(values, dict):
            gdf[as_column] = gdf[join_key].map(values).astype(float)
        else:
            merged = gdf.merge(
                values[[join_key, value_key]],
                on=join_key,
                how="left",
            )
            merged[as_column] = merged[value_key].astype(float)
            gdf = merged.drop(columns=[value_key])
        gdf[as_column] = gdf[as_column].fillna(0.0)
        return WorldMap(
            gdf=gdf,
            value_column=as_column,
            name_column=self.name_column,
            iso_column=self.iso_column,
        )

    # ------------------------------------------------------------------
    def pick(self, column: str) -> "WorldMap":
        """Declare which existing column is the cartogram value column."""
        if column not in self.gdf.columns:
            raise KeyError(f"{column!r} is not a column of the world map")
        return WorldMap(
            gdf=self.gdf,
            value_column=column,
            name_column=self.name_column,
            iso_column=self.iso_column,
        )

    # ------------------------------------------------------------------
    @property
    def bbox(self) -> tuple[float, float, float, float]:
        """Bounding box of the loaded geometries in the current CRS."""
        xmin, ymin, xmax, ymax = self.gdf.total_bounds
        return float(xmin), float(ymin), float(xmax), float(ymax)

    # ------------------------------------------------------------------
    def padded_bbox(self, pad: float = 0.03) -> tuple[float, float, float, float]:
        """Bounding box with a symmetric padding (fraction of extent)."""
        xmin, ymin, xmax, ymax = self.bbox
        dx = (xmax - xmin) * pad
        dy = (ymax - ymin) * pad
        return xmin - dx, ymin - dy, xmax + dx, ymax + dy


def csv_to_value_map(
    csv_path: str,
    iso_column: str,
    value_column: str,
    encoding: str = "utf-8",
) -> dict:
    """Parse a CSV into a ``{ISO3: float}`` mapping.

    Rows with missing / non-numeric values are skipped. This is the cheap,
    dependency-light ingestion path for the CLI; richer joins should use
    :meth:`WorldMap.with_values` directly with a pandas DataFrame.
    """
    import csv

    out: dict = {}
    with open(csv_path, newline="", encoding=encoding) as fh:
        reader = csv.DictReader(fh)
        if iso_column not in reader.fieldnames or value_column not in reader.fieldnames:
            raise KeyError(
                f"CSV must contain columns {iso_column!r} and {value_column!r}; "
                f"found {reader.fieldnames!r}"
            )
        for row in reader:
            iso = (row.get(iso_column) or "").strip()
            if not iso:
                continue
            raw = (row.get(value_column) or "").strip().replace(",", "")
            if not raw:
                continue
            try:
                out[iso] = float(raw)
            except ValueError:
                continue
    return out
