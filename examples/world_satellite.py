"""World cartogram with a Natural Earth shaded-relief background.

Downloads a small (~2 MB) public-domain Natural Earth raster on first
run (cached under ``~/.cache/cartogram/``), builds a population or GDP
cartogram in lat/lon coordinates, and warps the raster through the
same deformation so the satellite-like imagery deforms consistently
with the country polygons.

Run:

    python examples/world_satellite.py                   # population
    python examples/world_satellite.py --value gdp       # GDP
    python examples/world_satellite.py --value gdp_per_capita --log

Output: ``docs/images/world_satellite_<value>.png``.
"""

from __future__ import annotations

import argparse
import os
import urllib.request
import zipfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from cartogram import Cartogram, WorldMap, rasterize_polygons, warp_image


CACHE = Path.home() / ".cache" / "cartogram"
RASTER_URL = "https://naciscdn.org/naturalearth/50m/raster/NE1_50M_SR_W.zip"
RASTER_NAME = "NE1_50M_SR_W.tif"
OUT_DIR = Path(__file__).resolve().parent.parent / "docs" / "images"


def _ensure_raster() -> Path:
    CACHE.mkdir(parents=True, exist_ok=True)
    tif = CACHE / RASTER_NAME
    if tif.exists():
        return tif
    zpath = CACHE / "NE1_50M_SR_W.zip"
    print(f"Downloading {RASTER_URL} ...", flush=True)
    urllib.request.urlretrieve(RASTER_URL, zpath)
    with zipfile.ZipFile(zpath) as zf:
        zf.extractall(CACHE)
    if not tif.exists():
        # Some Natural Earth archives nest the tif under a subdirectory.
        matches = list(CACHE.rglob(RASTER_NAME))
        if matches:
            matches[0].replace(tif)
    if not tif.exists():
        raise RuntimeError(f"Could not find {RASTER_NAME} in the downloaded archive")
    return tif


def _load_raster_as_wgs84_rgb(tif: Path):
    """Return ``(rgb, bbox_xyxy)`` where ``rgb`` is (H, W, 3) and
    ``bbox_xyxy`` is ``(lon_min, lat_min, lon_max, lat_max)``.

    Tries rasterio if available (for correct geo-referencing). Falls
    back to assuming the raster is a global Plate Carrée image covering
    (-180, -90, 180, 90)."""
    try:
        import rasterio

        with rasterio.open(tif) as ds:
            img = ds.read()  # (bands, H, W)
            if img.shape[0] >= 3:
                rgb = np.stack([img[0], img[1], img[2]], axis=-1)
            else:
                gray = img[0]
                rgb = np.stack([gray, gray, gray], axis=-1)
            left, bottom, right, top = ds.bounds
            return rgb, (float(left), float(bottom), float(right), float(top))
    except Exception:
        pass

    arr = np.asarray(Image.open(tif).convert("RGB"))
    return arr, (-180.0, -90.0, 180.0, 90.0)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--value", choices=["population", "gdp", "gdp_per_capita"], default="population"
    )
    parser.add_argument("--grid", type=int, default=512)
    parser.add_argument("--log", action="store_true")
    parser.add_argument("--min-floor", type=float, default=None)
    parser.add_argument("--out", default=None)
    args = parser.parse_args()

    tif = _ensure_raster()
    rgb, raster_bbox = _load_raster_as_wgs84_rgb(tif)
    # Raster pixel (0, 0) is in the top-left (lat = +90). Flip so that
    # imshow with origin="lower" sees low-y at the bottom.
    rgb = rgb[::-1]

    # Load countries in WGS84 so the raster and polygons share a CRS.
    wm = WorldMap.load(projection="EPSG:4326").pick(args.value)
    vals = wm.gdf[args.value].astype(float).to_numpy()
    if args.min_floor is not None:
        floor = args.min_floor * float(np.nanmean(vals[vals > 0]))
        vals = np.maximum(vals, floor)
    if args.log:
        vals = np.log1p(np.maximum(vals, 0.0))
    gdf = wm.gdf.copy()
    gdf["_cart"] = vals

    # Clip the working bbox to the raster's extent.
    b = wm.padded_bbox(pad=0.01)
    xmin = max(b[0], raster_bbox[0])
    ymin = max(b[1], raster_bbox[1])
    xmax = min(b[2], raster_bbox[2])
    ymax = min(b[3], raster_bbox[3])
    bbox = (xmin, ymin, xmax, ymax)

    aspect = (xmax - xmin) / (ymax - ymin)
    nx = int(args.grid)
    ny = max(16, int(round(args.grid / aspect)))
    rho = rasterize_polygons(list(gdf.geometry), list(gdf["_cart"]), bbox, (ny, nx), subpixel=2)
    cart = Cartogram(rho, bbox=bbox, mean_floor=0.02, blur_sigma=1.5)
    print(f"advecting {ny*nx} grid points ...", flush=True)
    cart.run(tol=3e-3)

    new_geoms = cart.transform_polygons(list(gdf.geometry))
    print("warping raster through deformation ...", flush=True)
    warped = warp_image(cart, rgb, image_bbox=raster_bbox, out_bbox=bbox,
                        out_shape=(rgb.shape[0] // 2, rgb.shape[1] // 2))
    # Crop the original raster to the working bbox for side-by-side.
    from cartogram.warp import _inverse_transform  # noqa: F401  (not needed here)
    ry = rgb.shape[0]; rx = rgb.shape[1]
    def _crop(img, src_bbox, dst_bbox):
        sxmin, symin, sxmax, symax = src_bbox
        dxmin, dymin, dxmax, dymax = dst_bbox
        H, W = img.shape[:2]
        sx = W / (sxmax - sxmin)
        sy = H / (symax - symin)
        j0 = max(0, int(round((dxmin - sxmin) * sx)))
        j1 = min(W, int(round((dxmax - sxmin) * sx)))
        i0 = max(0, int(round((dymin - symin) * sy)))
        i1 = min(H, int(round((dymax - symin) * sy)))
        return img[i0:i1, j0:j1]

    orig_crop = _crop(rgb, raster_bbox, bbox)

    # -------- plot --------
    panel_w = 8.0
    fig_w = 2 * panel_w + 1
    fig_h = panel_w / aspect + 1.2
    fig, (a, b2) = plt.subplots(1, 2, figsize=(fig_w, fig_h))
    extent = (bbox[0], bbox[2], bbox[1], bbox[3])
    a.imshow(orig_crop, extent=extent, origin="lower")
    b2.imshow(warped, extent=extent, origin="lower")
    gdf.plot(ax=a, facecolor="none", edgecolor="black", linewidth=0.3)
    gnew = gdf.copy(); gnew["geometry"] = new_geoms
    gnew.plot(ax=b2, facecolor="none", edgecolor="black", linewidth=0.3)
    a.set_title(f"Geographic map + Natural Earth relief")
    b2.set_title(f"Cartogram (area ∝ {args.value})")
    for ax in (a, b2):
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out = args.out or str(OUT_DIR / f"world_satellite_{args.value}.png")
    fig.suptitle(f"World — {args.value} with shaded-relief background")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
