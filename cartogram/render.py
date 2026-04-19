"""Publication-ready matplotlib figures for cartograms.

The single useful function in here is :func:`side_by_side`, which
produces the two-panel geographic + cartogram figure used in the
notebook and demo images. It matches the Wolfram port's
``WorldCartogram`` output layout (same colour scale across both
panels, equal aspect, ``log`` colour-norm when the data is
heavy-tailed) and handles the grey-out path for the
``missing_countries="grey"`` mode.
"""

from __future__ import annotations

from typing import Any, Iterable, Optional

import numpy as np


GREY = "#d9d9d9"


def side_by_side(
    gdf,
    deformed: Iterable,
    label: str = "value",
    *,
    cmap: str = "plasma",
    value_column: str = "value",
    has_value_mask: Optional[np.ndarray] = None,
    figsize: tuple = (18, 6),
    edgecolor: str = "black",
    linewidth: float = 0.2,
    show_titles: bool = True,
):
    """Render a two-panel geographic + cartogram figure.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Source geometries with a numeric ``value_column``.
    deformed : list of shapely geometries
        Output of :meth:`Cartogram.transform_polygons` on
        ``list(gdf.geometry)`` (same order, same length).
    label : str
        Human-readable metric name used in the cartogram panel title.
    cmap : str
        Matplotlib colormap.
    value_column : str
        Column on ``gdf`` holding the metric values. Defaults to
        ``"value"`` to match what :func:`world_cartogram` writes.
    has_value_mask : np.ndarray[bool], optional
        When given, rows where the mask is False are drawn in neutral
        grey on both panels and omitted from the colour-scale
        calibration (matches the ``missing_countries="grey"`` mode).
    figsize, edgecolor, linewidth, show_titles :
        Forwarded to matplotlib. Defaults match the demo images.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object. Not shown; save with ``fig.savefig(...)`` or
        display it in whatever environment you prefer.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize
    import geopandas as gpd

    vals = gdf[value_column].to_numpy(dtype=float)
    if has_value_mask is None:
        mask = np.ones(len(vals), dtype=bool)
    else:
        mask = np.asarray(has_value_mask, dtype=bool)
        if mask.shape != (len(vals),):
            raise ValueError(
                f"has_value_mask shape {mask.shape} does not match "
                f"gdf length {len(vals)}"
            )

    finite = vals[mask & np.isfinite(vals) & (vals > 0)]
    if finite.size == 0:
        norm = Normalize(vmin=0, vmax=1)
    else:
        vmin = float(finite.min())
        vmax = float(finite.max())
        if vmax / max(vmin, 1e-12) > 50:
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)

    valued_gdf = gdf.loc[mask]
    missing_gdf = gdf.loc[~mask]

    fig, (ax_geo, ax_cart) = plt.subplots(
        1, 2, figsize=figsize, constrained_layout=True
    )

    if not missing_gdf.empty:
        missing_gdf.plot(
            ax=ax_geo, color=GREY,
            edgecolor=edgecolor, linewidth=linewidth,
        )
    if not valued_gdf.empty:
        valued_gdf.plot(
            column=value_column, ax=ax_geo, cmap=cmap, norm=norm,
            edgecolor=edgecolor, linewidth=linewidth,
        )
    if show_titles:
        ax_geo.set_title("Geographic")
    ax_geo.set_aspect("equal")
    ax_geo.set_xticks([]); ax_geo.set_yticks([])

    deformed_list = list(deformed)
    if len(deformed_list) != len(gdf):
        raise ValueError(
            f"deformed has {len(deformed_list)} geometries but gdf has "
            f"{len(gdf)} rows"
        )
    missing_idx = np.where(~mask)[0]
    valued_idx = np.where(mask)[0]
    if missing_idx.size:
        missing_def = gpd.GeoDataFrame(
            {value_column: vals[~mask]},
            geometry=[deformed_list[i] for i in missing_idx],
            crs=gdf.crs,
        )
        missing_def.plot(
            ax=ax_cart, color=GREY,
            edgecolor=edgecolor, linewidth=linewidth,
        )
    if valued_idx.size:
        valued_def = gpd.GeoDataFrame(
            {value_column: vals[mask]},
            geometry=[deformed_list[i] for i in valued_idx],
            crs=gdf.crs,
        )
        valued_def.plot(
            column=value_column, ax=ax_cart, cmap=cmap, norm=norm,
            edgecolor=edgecolor, linewidth=linewidth,
        )
    if show_titles:
        ax_cart.set_title(
            f"Cartogram (area \u221d {label}, oceans preserved)"
        )
    ax_cart.set_aspect("equal")
    ax_cart.set_xticks([]); ax_cart.set_yticks([])

    return fig
