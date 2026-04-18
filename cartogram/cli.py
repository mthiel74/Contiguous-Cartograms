"""Command-line interface for the cartogram package.

Typical uses:

    # Built-in world map: resize countries by population.
    python -m cartogram world --value population --out pop.png

    # Same, but by GDP, at higher resolution.
    python -m cartogram world --value gdp --grid 1024 --out gdp.png

    # Join your own CSV (must have an ISO_A3 column).
    python -m cartogram world --csv my_data.csv \\
        --csv-iso-col iso3 --csv-value-col value --out mine.png

    # Add a satellite-like background image.
    python -m cartogram world --value population \\
        --background docs/earth.jpg --out pop_sat.png

    # Cartogram of an arbitrary shapefile + column.
    python -m cartogram shape path/to/states.shp --value pop2020 \\
        --out states.png
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Optional, Sequence

import numpy as np


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="python -m cartogram",
        description="Build diffusion-based contiguous cartograms.",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    w = sub.add_parser("world", help="World cartogram (Natural Earth admin-0).")
    w.add_argument(
        "--value",
        choices=["population", "gdp", "gdp_per_capita", "custom"],
        default="population",
        help="Which built-in quantity to visualise (default: population). "
        "Use 'custom' together with --csv.",
    )
    w.add_argument("--csv", help="CSV file with per-country values.")
    w.add_argument("--csv-iso-col", default="ISO_A3")
    w.add_argument("--csv-value-col", default="value")
    w.add_argument(
        "--resolution",
        choices=["110m", "50m"],
        default="110m",
        help="Natural Earth source resolution (default 110m).",
    )
    w.add_argument(
        "--grid",
        type=int,
        default=512,
        help="Density-grid side length in cells (default 512).",
    )
    w.add_argument(
        "--background",
        default=None,
        help="Optional raster image (PNG/JPG) to warp underneath the "
        "cartogram. Must cover the world bbox in the selected CRS; for "
        "equirectangular/Plate Carrée imagery use --projection EPSG:4326.",
    )
    w.add_argument(
        "--projection",
        default="EPSG:8857",
        help="Target CRS for the world map. Default is Equal Earth.",
    )
    w.add_argument("--out", required=True, help="Output PNG path.")
    w.add_argument(
        "--dpi", type=int, default=160, help="Figure DPI (default 160)."
    )
    w.add_argument(
        "--title", default=None, help="Figure title override."
    )
    w.add_argument(
        "--min-floor",
        type=float,
        default=None,
        help="Clamp per-country values from below to this fraction of the "
        "global mean. Stabilises cartograms whose value has a very heavy "
        "tail (e.g. GDP per capita). Default: no clamp.",
    )
    w.add_argument(
        "--log",
        action="store_true",
        help="Apply log1p scaling to the value column before rasterising. "
        "Produces gentler deformations for heavy-tailed quantities.",
    )
    w.add_argument(
        "--preserve-oceans",
        action="store_true",
        help="Fill ocean cells with the mean land density so that ocean "
        "area is approximately preserved; land masses are only resized "
        "relative to each other.",
    )

    s = sub.add_parser("shape", help="Cartogram from an arbitrary shapefile.")
    s.add_argument("path", help="Shapefile or other geopandas-readable file.")
    s.add_argument("--value", required=True, help="Column to use as value.")
    s.add_argument(
        "--projection",
        default=None,
        help="Optional CRS to reproject into before building the cartogram.",
    )
    s.add_argument("--grid", type=int, default=512)
    s.add_argument("--out", required=True)
    s.add_argument("--dpi", type=int, default=160)
    s.add_argument("--title", default=None)
    s.add_argument("--preserve-oceans", action="store_true")

    g = sub.add_parser("gui", help="Launch the Gradio web interface.")
    g.add_argument("--host", default="127.0.0.1")
    g.add_argument("--port", type=int, default=7860)
    g.add_argument("--share", action="store_true",
                   help="Expose via a temporary gradio.live tunnel.")

    return p


# ----------------------------------------------------------------------
def _build_cartogram(gdf, value_col: str, bbox, grid: int, preserve_oceans: bool = False):
    from .cartogram import Cartogram
    from .rasterize import rasterize_polygons

    xmin, ymin, xmax, ymax = bbox
    Lx = xmax - xmin
    Ly = ymax - ymin
    aspect = Lx / Ly
    if aspect >= 1.0:
        nx = int(grid)
        ny = max(16, int(round(grid / aspect)))
    else:
        ny = int(grid)
        nx = max(16, int(round(grid * aspect)))
    polys = list(gdf.geometry)
    values = [float(v) for v in gdf[value_col].to_numpy()]
    rho = rasterize_polygons(polys, values, bbox, (ny, nx), subpixel=2)
    cart = Cartogram(
        rho,
        bbox=bbox,
        mean_floor=0.02,
        blur_sigma=1.5,
        sea_density="auto" if preserve_oceans else None,
    )
    cart.run(tol=2e-3)
    return cart


def _plot(
    gdf,
    new_geoms,
    value_col: str,
    bbox,
    out_path: str,
    dpi: int,
    title: Optional[str],
    background_image: Optional[np.ndarray] = None,
    cartogram_label: Optional[str] = None,
):
    import matplotlib.pyplot as plt

    xmin, ymin, xmax, ymax = bbox
    aspect = (xmax - xmin) / (ymax - ymin)
    panel_width = 7.0
    fig_width = 2 * panel_width + 1.0
    fig_height = panel_width / aspect + 1.4
    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(fig_width, fig_height))
    if background_image is not None:
        # Original background on the left, cartogram-warped on the right.
        extent = (bbox[0], bbox[2], bbox[1], bbox[3])
        ax0.imshow(background_image[0], extent=extent, origin="lower", zorder=0)
        ax1.imshow(background_image[1], extent=extent, origin="lower", zorder=0)

    gdf.plot(
        ax=ax0,
        column=value_col,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.3,
        alpha=0.85 if background_image is not None else 1.0,
        legend=True,
        legend_kwds={"shrink": 0.6},
    )
    new_gdf = gdf.copy()
    new_gdf["geometry"] = new_geoms
    new_gdf.plot(
        ax=ax1,
        column=value_col,
        cmap="viridis",
        edgecolor="black",
        linewidth=0.3,
        alpha=0.85 if background_image is not None else 1.0,
        legend=True,
        legend_kwds={"shrink": 0.6},
    )
    ax0.set_title("Geographic map")
    ax1.set_title(f"Cartogram (area ∝ {cartogram_label or value_col})")
    for ax in (ax0, ax1):
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


# ----------------------------------------------------------------------
def _cmd_world(args) -> int:
    from .world_data import WorldMap, csv_to_value_map

    wm = WorldMap.load(
        resolution=args.resolution,
        projection=args.projection,
        drop_antarctica=True,
    )

    if args.value == "custom":
        if not args.csv:
            print("error: --value custom requires --csv", file=sys.stderr)
            return 2
        val_map = csv_to_value_map(args.csv, args.csv_iso_col, args.csv_value_col)
        wm = wm.with_values(val_map, join_key=wm.iso_column, as_column="value")
        value_col = "value"
    else:
        wm = wm.pick(args.value)
        value_col = args.value

    display_col = value_col
    if args.log or args.min_floor is not None:
        gdf = wm.gdf.copy()
        vals = gdf[value_col].astype(float).to_numpy()
        if args.min_floor is not None:
            floor = args.min_floor * float(np.nanmean(vals[vals > 0]))
            vals = np.maximum(vals, floor)
        if args.log:
            vals = np.log1p(np.maximum(vals, 0.0))
        shaped_col = "_cartogram_value"
        gdf[shaped_col] = vals
        wm = type(wm)(gdf=gdf, name_column=wm.name_column, iso_column=wm.iso_column).pick(shaped_col)
        value_col = shaped_col

    bbox = wm.padded_bbox(pad=0.02)
    cart = _build_cartogram(
        wm.gdf,
        value_col,
        bbox,
        args.grid,
        preserve_oceans=args.preserve_oceans,
    )
    new_geoms = cart.transform_polygons(list(wm.gdf.geometry))

    background = None
    if args.background:
        from PIL import Image
        from .warp import warp_image

        img = np.asarray(Image.open(args.background).convert("RGB"))
        # Flip so pixel (0,0) is at ymin (matplotlib extent convention).
        img = img[::-1]
        warped = warp_image(cart, img, image_bbox=bbox, out_bbox=bbox)
        background = (img, warped)

    title = args.title or f"World cartogram — area ∝ {display_col}"
    plot_col = display_col if display_col in wm.gdf.columns else value_col
    _plot(
        wm.gdf,
        new_geoms,
        plot_col,
        bbox,
        args.out,
        args.dpi,
        title,
        background_image=background,
        cartogram_label=display_col,
    )
    print(f"wrote {args.out}")
    return 0


def _cmd_shape(args) -> int:
    import geopandas as gpd

    gdf = gpd.read_file(args.path)
    if args.projection:
        gdf = gdf.to_crs(args.projection)
    if args.value not in gdf.columns:
        print(
            f"error: column {args.value!r} not in shapefile; "
            f"available: {list(gdf.columns)}",
            file=sys.stderr,
        )
        return 2

    xmin, ymin, xmax, ymax = gdf.total_bounds
    dx = (xmax - xmin) * 0.02
    dy = (ymax - ymin) * 0.02
    bbox = (xmin - dx, ymin - dy, xmax + dx, ymax + dy)

    cart = _build_cartogram(
        gdf, args.value, bbox, args.grid,
        preserve_oceans=args.preserve_oceans,
    )
    new_geoms = cart.transform_polygons(list(gdf.geometry))

    title = args.title or f"Cartogram — area ∝ {args.value}"
    _plot(gdf, new_geoms, args.value, bbox, args.out, args.dpi, title)
    print(f"wrote {args.out}")
    return 0


def _cmd_gui(args) -> int:
    from .gui import launch

    launch(host=args.host, port=args.port, share=args.share)
    return 0


# ----------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _build_parser().parse_args(argv)
    if args.cmd == "world":
        return _cmd_world(args)
    if args.cmd == "shape":
        return _cmd_shape(args)
    if args.cmd == "gui":
        return _cmd_gui(args)
    return 2
