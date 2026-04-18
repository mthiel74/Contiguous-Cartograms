"""Optional demo: build a population cartogram of the contiguous US states.

This script is *optional*: it requires ``geopandas`` and downloads a
US-states shapefile (from the US Census) plus a simple population CSV
on first run, caching under ``examples/data/``. It mirrors the shape
of :mod:`examples.synthetic` but operates on real geography.

Run:

    python examples/us_states.py

Output: ``examples/output/us_states.png``.

If you don't have ``geopandas`` installed, use ``examples/synthetic.py``
instead — it has no GIS dependency.
"""

from __future__ import annotations

import os
import sys
import urllib.request
import zipfile

import matplotlib.pyplot as plt
import numpy as np

from cartogram import Cartogram, rasterize_polygons


HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")
OUT_DIR = os.path.join(HERE, "output")

STATES_URL = "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_state_20m.zip"
STATES_ZIP = os.path.join(DATA_DIR, "cb_2022_us_state_20m.zip")
STATES_SHP = os.path.join(DATA_DIR, "cb_2022_us_state_20m.shp")


def _download_if_needed() -> None:
    os.makedirs(DATA_DIR, exist_ok=True)
    if os.path.exists(STATES_SHP):
        return
    print(f"Downloading {STATES_URL} ...")
    urllib.request.urlretrieve(STATES_URL, STATES_ZIP)
    with zipfile.ZipFile(STATES_ZIP) as zf:
        zf.extractall(DATA_DIR)


def main() -> None:
    try:
        import geopandas as gpd
    except ImportError:
        print(
            "This example requires geopandas. Install with:\n"
            "    pip install geopandas",
            file=sys.stderr,
        )
        sys.exit(1)

    _download_if_needed()
    gdf = gpd.read_file(STATES_SHP)

    # Continental US only (drop Alaska, Hawaii, territories).
    keep = gdf["STUSPS"].isin(
        [s for s in gdf["STUSPS"].unique() if s not in {"AK", "HI", "PR", "VI", "GU", "MP", "AS"}]
    )
    gdf = gdf[keep].to_crs("EPSG:5070")  # Albers Equal-Area (USA).

    # A tiny hard-coded 2020 census population table (millions) — good enough
    # for a demo. Avoids an additional network dependency.
    pop = {
        "AL": 5.03, "AR": 3.01, "AZ": 7.15, "CA": 39.50, "CO": 5.77, "CT": 3.60,
        "DE": 0.99, "FL": 21.54, "GA": 10.71, "IA": 3.19, "ID": 1.84,
        "IL": 12.81, "IN": 6.79, "KS": 2.94, "KY": 4.51, "LA": 4.65,
        "MA": 7.03, "MD": 6.18, "ME": 1.36, "MI": 10.08, "MN": 5.71,
        "MO": 6.16, "MS": 2.96, "MT": 1.08, "NC": 10.44, "ND": 0.78,
        "NE": 1.96, "NH": 1.38, "NJ": 9.29, "NM": 2.12, "NV": 3.10,
        "NY": 20.20, "OH": 11.80, "OK": 3.96, "OR": 4.24, "PA": 13.00,
        "RI": 1.10, "SC": 5.12, "SD": 0.89, "TN": 6.91, "TX": 29.15,
        "UT": 3.27, "VA": 8.63, "VT": 0.64, "WA": 7.71, "WI": 5.89,
        "WV": 1.79, "WY": 0.58, "DC": 0.69,
    }
    gdf = gdf[gdf["STUSPS"].isin(pop.keys())].copy()
    gdf["pop"] = gdf["STUSPS"].map(pop)

    xmin, ymin, xmax, ymax = gdf.total_bounds
    # Add 2% padding so the cartogram has room to breathe.
    pad_x = 0.02 * (xmax - xmin)
    pad_y = 0.02 * (ymax - ymin)
    bbox = (xmin - pad_x, ymin - pad_y, xmax + pad_x, ymax + pad_y)

    ny, nx = 512, 512
    rho = rasterize_polygons(
        list(gdf.geometry),
        list(gdf["pop"]),
        bbox,
        (ny, nx),
        subpixel=2,
    )

    cart = Cartogram(rho, bbox=bbox, mean_floor=0.01)
    cart.run(tol=5e-3)

    new_geoms = cart.transform_polygons(list(gdf.geometry))
    new_gdf = gdf.copy()
    new_gdf["geometry"] = new_geoms

    os.makedirs(OUT_DIR, exist_ok=True)
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(14, 7))
    gdf.plot(ax=ax0, column="pop", cmap="viridis", edgecolor="black", linewidth=0.3)
    new_gdf.plot(ax=ax1, column="pop", cmap="viridis", edgecolor="black", linewidth=0.3)
    ax0.set_title("Lower-48 — equal-area projection")
    ax1.set_title("Lower-48 — population cartogram")
    for ax in (ax0, ax1):
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
    out_path = os.path.join(OUT_DIR, "us_states.png")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
