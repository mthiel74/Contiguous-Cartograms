"""Runnable demo: a 4x4 grid of square regions with a non-uniform
population, turned into a contiguous cartogram.

Output: ``examples/output/synthetic.png`` — two panels (original /
cartogram) filled and colour-coded by population.
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
from shapely.geometry import Polygon

from cartogram import Cartogram, rasterize_polygons


# ------------------------------------------------------------------
# Build a 4x4 grid of unit squares with a made-up population field.
def build_regions(ncols: int = 4, nrows: int = 4):
    polys = []
    pops = []
    rng = np.random.default_rng(0)
    # Pick a single "hot" cell that is ~10x denser than the rest.
    hot = (nrows // 2, ncols // 2)
    for i in range(nrows):
        for j in range(ncols):
            polys.append(
                Polygon(
                    [
                        (j, i),
                        (j + 1, i),
                        (j + 1, i + 1),
                        (j, i + 1),
                    ]
                )
            )
            if (i, j) == hot:
                pops.append(12.0)
            else:
                pops.append(1.0 + 0.3 * rng.random())
    return polys, np.array(pops)


def main() -> None:
    polys, pops = build_regions()
    ncols = 4
    nrows = 4
    bbox = (0.0, 0.0, float(ncols), float(nrows))

    # Rasterise to a 256 x 256 density grid.
    ny, nx = 256, 256
    rho = rasterize_polygons(polys, pops, bbox, (ny, nx), subpixel=3)

    cart = Cartogram(rho, bbox=bbox, mean_floor=0.01)
    cart.run(tol=1e-3)

    new_polys = cart.transform_polygons(polys)

    # --- plot ---
    out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
    os.makedirs(out_dir, exist_ok=True)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(10, 5))
    cmap = plt.get_cmap("viridis")
    vmin, vmax = pops.min(), pops.max()
    for ax, polygons, title in (
        (ax0, polys, "Original (equal-area)"),
        (ax1, new_polys, "Cartogram (area ∝ population)"),
    ):
        for poly, val in zip(polygons, pops):
            colour = cmap((val - vmin) / (vmax - vmin + 1e-12))
            xs, ys = poly.exterior.xy
            ax.fill(xs, ys, color=colour, edgecolor="black", linewidth=0.5)
        ax.set_aspect("equal")
        ax.set_title(title)
        ax.set_xlim(-0.2, ncols + 0.2)
        ax.set_ylim(-0.2, nrows + 0.2)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle("Diffusion cartogram — 4x4 grid with one dense cell")
    fig.tight_layout()

    out_path = os.path.join(out_dir, "synthetic.png")
    fig.savefig(out_path, dpi=150)
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
