"""Gradio web GUI for building world cartograms from any indicator.

Run with either

    python -m cartogram gui

or directly

    python -m cartogram.gui

The UI exposes four data paths — built-in (Natural Earth attributes),
World Bank API, WHO GHO API, and CSV upload — wraps the core algorithm
with sane defaults, and offers a single toggle for ocean-area
preservation as well as optional shaded-relief backgrounds.
"""

from __future__ import annotations

import io
import tempfile
from pathlib import Path
from typing import Optional, Tuple

import numpy as np

from .cartogram import Cartogram
from .data_sources import INDICATOR_CATALOGUE, catalogue_labels, fetch, indicator_by_label
from .rasterize import rasterize_polygons
from .world_data import WorldMap


DEFAULT_PROJECTION = "EPSG:8857"   # Equal Earth
BUILTIN_OPTIONS = ["population", "gdp", "gdp_per_capita"]


# ----------------------------------------------------------------------
def _apply_transforms(
    gdf, column: str, log_scale: bool, min_floor: float
):
    import pandas as pd

    values = gdf[column].astype(float).to_numpy()
    if min_floor and min_floor > 0:
        land = values[values > 0]
        if land.size:
            floor = float(min_floor) * float(land.mean())
            values = np.maximum(values, floor)
    if log_scale:
        values = np.log1p(np.maximum(values, 0.0))
    out = gdf.copy()
    out["_cart_value"] = values
    return out


def _compute_ocean_share(gdf_orig, gdf_new, bbox) -> tuple[float, float]:
    """Return (orig_ocean_share, new_ocean_share) as fractions of the bbox."""
    from shapely.geometry import box
    from shapely.ops import unary_union

    B = box(*bbox)
    land0 = unary_union(list(gdf_orig.geometry)).intersection(B)
    land1 = unary_union(list(gdf_new.geometry)).intersection(B)
    total = B.area
    return (1.0 - land0.area / total, 1.0 - land1.area / total)


# ----------------------------------------------------------------------
def _render(
    gdf, new_geoms, label: str, bbox, background_image=None
):
    import matplotlib.pyplot as plt

    xmin, ymin, xmax, ymax = bbox
    aspect = (xmax - xmin) / (ymax - ymin)
    panel = 6.5
    fig, (a, b) = plt.subplots(
        1, 2, figsize=(2 * panel + 1.0, panel / aspect + 1.2)
    )

    if background_image is not None:
        extent = (xmin, xmax, ymin, ymax)
        a.imshow(background_image[0], extent=extent, origin="lower")
        b.imshow(background_image[1], extent=extent, origin="lower")
        alpha = 0.35
    else:
        alpha = 0.9

    gdf.plot(
        ax=a, column="_cart_value" if "_cart_value" in gdf.columns else None,
        cmap="viridis", edgecolor="black", linewidth=0.3, alpha=alpha,
    )
    gnew = gdf.copy()
    gnew["geometry"] = new_geoms
    gnew.plot(
        ax=b, column="_cart_value" if "_cart_value" in gnew.columns else None,
        cmap="viridis", edgecolor="black", linewidth=0.3, alpha=alpha,
    )
    a.set_title("Geographic map")
    b.set_title(f"Cartogram — area ∝ {label}")
    for ax in (a, b):
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=140)
    plt.close(fig)
    buf.seek(0)
    from PIL import Image
    return Image.open(buf).convert("RGB")


# ----------------------------------------------------------------------
def _load_background():
    """Return ``(rgb, (lon_min, lat_min, lon_max, lat_max))`` or ``None``.

    Reuses the NE1_50M_SR_W shaded-relief raster if examples/world_satellite
    has already downloaded it; otherwise returns None."""
    from pathlib import Path

    tif = Path.home() / ".cache" / "cartogram" / "NE1_50M_SR_W.tif"
    if not tif.exists():
        # Try to fetch it quietly.
        try:
            from examples.world_satellite import _ensure_raster
            tif = _ensure_raster()
        except Exception:
            return None
    try:
        import rasterio
        with rasterio.open(tif) as ds:
            img = ds.read()
            if img.shape[0] >= 3:
                rgb = np.stack([img[0], img[1], img[2]], axis=-1)
            else:
                g = img[0]
                rgb = np.stack([g, g, g], axis=-1)
            l, bt, r, t = ds.bounds
            return rgb[::-1], (float(l), float(bt), float(r), float(t))
    except Exception:
        return None


# ----------------------------------------------------------------------
def generate(
    source: str,
    builtin_column: str,
    indicator_label: str,
    year: int,
    csv_file,
    csv_iso_column: str,
    csv_value_column: str,
    grid: int,
    log_scale: bool,
    min_floor: float,
    preserve_oceans: bool,
    background: str,
    projection: str,
) -> Tuple["object", str]:
    """Main callback wired to the "Generate cartogram" button."""
    use_wgs84 = background == "Shaded relief (Natural Earth)"
    crs = "EPSG:4326" if use_wgs84 else (projection or DEFAULT_PROJECTION)

    wm = WorldMap.load(projection=crs)

    # --- source resolution ---
    if source == "Built-in (Natural Earth)":
        if builtin_column not in BUILTIN_OPTIONS:
            return None, f"error: unknown built-in column {builtin_column!r}"
        wm = wm.pick(builtin_column)
        display_label = builtin_column
    elif source == "World Bank / WHO (catalogue)":
        ind = indicator_by_label(indicator_label)
        values = fetch(ind.source, ind.code, year=year or None)
        wm = wm.with_values(values, join_key=wm.iso_column, as_column="value")
        wm = wm.pick("value")
        display_label = f"{ind.label} ({ind.unit or '-'})"
    elif source == "Upload CSV":
        if csv_file is None:
            return None, "Please upload a CSV file first."
        import pandas as pd
        df = pd.read_csv(csv_file.name if hasattr(csv_file, "name") else csv_file)
        if csv_iso_column not in df.columns or csv_value_column not in df.columns:
            return None, (
                f"CSV is missing columns {csv_iso_column!r}/{csv_value_column!r}; "
                f"found: {list(df.columns)}"
            )
        val_map = {
            str(row[csv_iso_column]).strip(): float(row[csv_value_column])
            for _, row in df.iterrows()
            if str(row.get(csv_iso_column, "")).strip()
            and not pd.isna(row.get(csv_value_column))
        }
        wm = wm.with_values(val_map, join_key=wm.iso_column, as_column="value")
        wm = wm.pick("value")
        display_label = f"{csv_value_column} (uploaded)"
    else:
        return None, f"Unknown source {source!r}"

    gdf = _apply_transforms(wm.gdf, wm.value_column, log_scale, min_floor)

    # --- cartogram build ---
    bbox = wm.padded_bbox(pad=0.02)
    background_image = None
    if use_wgs84:
        bg = _load_background()
        if bg is not None:
            rgb, raster_bbox = bg
            xmin = max(bbox[0], raster_bbox[0]); ymin = max(bbox[1], raster_bbox[1])
            xmax = min(bbox[2], raster_bbox[2]); ymax = min(bbox[3], raster_bbox[3])
            bbox = (xmin, ymin, xmax, ymax)

    xmin, ymin, xmax, ymax = bbox
    aspect = (xmax - xmin) / (ymax - ymin)
    nx = int(grid)
    ny = max(16, int(round(grid / aspect)))

    rho = rasterize_polygons(
        list(gdf.geometry), list(gdf["_cart_value"]),
        bbox, (ny, nx), subpixel=2,
    )
    cart = Cartogram(
        rho, bbox=bbox, mean_floor=0.02, blur_sigma=1.5,
        sea_density="auto" if preserve_oceans else None,
    )
    cart.run(tol=2e-3)
    new_geoms = cart.transform_polygons(list(gdf.geometry))

    if use_wgs84 and background_image is None:
        bg = _load_background()
        if bg is not None:
            from .warp import warp_image
            rgb, raster_bbox = bg
            warped = warp_image(
                cart, rgb, image_bbox=raster_bbox, out_bbox=bbox,
                out_shape=(rgb.shape[0] // 2, rgb.shape[1] // 2),
            )
            # Crop the original to the working bbox.
            H, W = rgb.shape[:2]
            sx = W / (raster_bbox[2] - raster_bbox[0])
            sy = H / (raster_bbox[3] - raster_bbox[1])
            j0 = max(0, int(round((xmin - raster_bbox[0]) * sx)))
            j1 = min(W, int(round((xmax - raster_bbox[0]) * sx)))
            i0 = max(0, int(round((ymin - raster_bbox[1]) * sy)))
            i1 = min(H, int(round((ymax - raster_bbox[1]) * sy)))
            orig_crop = rgb[i0:i1, j0:j1]
            background_image = (orig_crop, warped)

    # --- ocean-share diagnostic ---
    try:
        gnew = gdf.copy(); gnew["geometry"] = new_geoms
        ob0, ob1 = _compute_ocean_share(gdf, gnew, bbox)
        ocean_msg = (
            f"Ocean area share — original: **{ob0:.1%}**, cartogram: **{ob1:.1%}**"
        )
    except Exception as e:
        ocean_msg = f"(ocean share computation failed: {e})"

    img = _render(gdf, new_geoms, display_label, bbox, background_image=background_image)

    info = (
        f"**Source:** {source}\n\n"
        f"**Indicator:** {display_label}\n\n"
        f"**Countries with data:** {int((gdf['_cart_value'] > 0).sum())} / {len(gdf)}\n\n"
        f"**Grid:** {ny} × {nx}   **log:** {log_scale}   "
        f"**min-floor:** {min_floor}   **preserve-oceans:** {preserve_oceans}\n\n"
        f"{ocean_msg}\n"
    )
    return img, info


# ----------------------------------------------------------------------
def build_interface():
    import gradio as gr

    with gr.Blocks(title="Contiguous Cartograms") as demo:
        gr.Markdown("# Contiguous Cartograms")
        gr.Markdown(
            "Resize every country by any quantity — built-in Natural Earth "
            "attributes, live World Bank / WHO indicators, or your own CSV. "
            "The algorithm is the classic Gastner–Newman diffusion cartogram "
            "(PNAS 2004). See the README for the math."
        )
        with gr.Row():
            with gr.Column(scale=1):
                source = gr.Radio(
                    ["Built-in (Natural Earth)",
                     "World Bank / WHO (catalogue)",
                     "Upload CSV"],
                    value="Built-in (Natural Earth)",
                    label="Data source",
                )
                builtin_column = gr.Radio(
                    BUILTIN_OPTIONS, value="population",
                    label="Built-in column", visible=True,
                )
                indicator_label = gr.Dropdown(
                    choices=catalogue_labels(),
                    value=catalogue_labels()[0],
                    label="Indicator (World Bank / WHO)",
                    visible=False,
                )
                year = gr.Number(
                    value=2022, label="Year (0 = latest available)",
                    precision=0, visible=False,
                )
                csv_file = gr.File(label="CSV file", visible=False,
                                   file_types=[".csv"])
                csv_iso_column = gr.Textbox(value="ISO_A3",
                                            label="CSV: ISO3 column",
                                            visible=False)
                csv_value_column = gr.Textbox(value="value",
                                              label="CSV: value column",
                                              visible=False)

                def _on_source_change(s):
                    import gradio as g
                    return (
                        g.update(visible=(s == "Built-in (Natural Earth)")),
                        g.update(visible=(s == "World Bank / WHO (catalogue)")),
                        g.update(visible=(s == "World Bank / WHO (catalogue)")),
                        g.update(visible=(s == "Upload CSV")),
                        g.update(visible=(s == "Upload CSV")),
                        g.update(visible=(s == "Upload CSV")),
                    )
                source.change(
                    _on_source_change, inputs=source,
                    outputs=[builtin_column, indicator_label, year,
                             csv_file, csv_iso_column, csv_value_column],
                )

                gr.Markdown("### Cartogram settings")
                grid = gr.Slider(128, 1024, value=384, step=32,
                                 label="Density grid side")
                log_scale = gr.Checkbox(label="Apply log1p to values",
                                        value=False)
                min_floor = gr.Slider(0.0, 1.0, value=0.0, step=0.01,
                                      label="Minimum value as fraction of mean")
                preserve_oceans = gr.Checkbox(
                    label="Preserve ocean area (recommended)", value=True,
                )
                background = gr.Radio(
                    ["None", "Shaded relief (Natural Earth)"],
                    value="None", label="Background",
                )
                projection = gr.Textbox(
                    value=DEFAULT_PROJECTION,
                    label="Projection CRS (ignored with shaded relief)",
                )

                generate_btn = gr.Button("Generate cartogram",
                                         variant="primary")

            with gr.Column(scale=2):
                out_image = gr.Image(label="Cartogram", type="pil",
                                     height=520)
                out_text = gr.Markdown()

        generate_btn.click(
            fn=generate,
            inputs=[source, builtin_column, indicator_label, year,
                    csv_file, csv_iso_column, csv_value_column,
                    grid, log_scale, min_floor, preserve_oceans,
                    background, projection],
            outputs=[out_image, out_text],
        )

    return demo


def launch(host: str = "127.0.0.1", port: int = 7860, share: bool = False):
    """Launch the Gradio GUI."""
    demo = build_interface()
    demo.launch(server_name=host, server_port=port, share=share)


if __name__ == "__main__":
    launch()
