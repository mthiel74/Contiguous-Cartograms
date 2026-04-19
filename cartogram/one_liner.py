"""A single-call driver over the Natural Earth world map.

Mirrors the Wolfram ``WorldCartogram[metric]`` helper: one function,
sensible defaults, returns a cartogram object and optionally a
side-by-side matplotlib figure. Built on top of the existing
:class:`cartogram.Cartogram`, :class:`cartogram.WorldMap`, and
:func:`cartogram.rasterize_polygons` primitives.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Literal, Mapping, Optional, Sequence, Tuple, Union

import numpy as np

from .cartogram import Cartogram
from .rasterize import rasterize_polygons


Metric = Union[str, Mapping[str, float], Callable[[Any], float]]
BBox = Tuple[float, float, float, float]


@dataclass
class WorldCartogramResult:
    """Bundle returned by :func:`world_cartogram`.

    Attributes
    ----------
    cartogram : Cartogram
        The driven cartogram, with ``.run()`` already called so it is
        ready for ``.transform`` / ``.transform_polygons`` calls.
    world_map : WorldMap
        The underlying world map with the metric attached as a column
        named ``"value"``.
    deformed_geometries : list
        Each country's exterior + interior rings after the cartogram
        deformation, as a list of shapely geometries in the same order
        as ``world_map.gdf``.
    figure : matplotlib.figure.Figure or None
        The side-by-side geographic + cartogram figure, if
        ``render=True`` was passed.
    """

    cartogram: Cartogram
    world_map: Any  # WorldMap — late import to avoid GIS stack at import time
    deformed_geometries: list
    figure: Any = None


def world_cartogram(
    metric: Metric,
    *,
    label: Optional[str] = None,
    resolution: str = "110m",
    projection: str = "EPSG:8857",
    drop_antarctica: bool = True,
    merge_sovereign: bool = False,
    grid_size: Tuple[int, int] = (256, 512),
    bbox: Optional[BBox] = None,
    mean_floor: float = 0.02,
    blur_sigma: float = 1.5,
    sea_density: Union[str, float, None] = "auto",
    min_floor_fraction: Optional[float] = None,
    rho_ceil_multiplier: Optional[float] = None,
    performance_goal: Literal["speed", "quality"] = "speed",
    snapshots: int = 60,
    tol: float = 1e-3,
    render: bool = True,
    cmap: str = "plasma",
) -> WorldCartogramResult:
    """Build a world cartogram in one call.

    Parameters
    ----------
    metric : str, dict, or callable
        What to drive the cartogram with.

        * **str** — a column name on the Natural Earth admin-0 layer
          (``"population"``, ``"gdp"``, ``"gdp_per_capita"``) or on the
          raw shapefile (``"POP_EST"``, ``"GDP_MD"``, …).
        * **dict** — ``{ISO3 or sovereign code: value}``. Countries not
          present take the floor value.
        * **callable** — ``f(row) -> float`` applied to each GeoDataFrame
          row; useful for derived quantities.
    label : str, optional
        Display label for the cartogram panel. Defaults to the metric
        name (for ``str``) or ``"custom metric"`` (otherwise).
    resolution, projection, drop_antarctica : str, str, bool
        Passed through to :meth:`WorldMap.load`.
    merge_sovereign : bool
        If True, call :meth:`WorldMap.merge_by_sovereign` so that
        overseas territories (Greenland under Denmark, French Guiana
        under France, Alaska/Hawaii-style behaviour for multi-part
        sovereigns) are unioned into a single row per sovereign state.
    grid_size : (ny, nx)
        Rasterisation grid. Higher is smoother, slower.
    bbox : (xmin, ymin, xmax, ymax), optional
        Physical extent to rasterise. Defaults to the padded map bbox.
    mean_floor, blur_sigma, sea_density :
        Passed to :class:`Cartogram`. Default ``sea_density="auto"``
        preserves ocean area — the overwhelmingly most useful mode.
    min_floor_fraction, rho_ceil_multiplier : float, optional
        Tame heavy-tailed distributions. If set, clamp per-country
        values below ``fraction * mean(values)`` and cap rasterised
        cell density at ``multiplier * mean(positive cells)`` before
        handing the grid to the solver.
    performance_goal, snapshots, tol : passed to :meth:`Cartogram.run`.
    render : bool
        If True, also build a matplotlib figure with the geographic map
        on the left and the cartogram on the right, and return it on
        the result.
    cmap : str
        Colormap for both panels when ``render=True``.

    Returns
    -------
    WorldCartogramResult
        The cartogram, the annotated world map, the deformed geometries,
        and (if requested) the matplotlib figure.

    Examples
    --------
    >>> from cartogram import world_cartogram
    >>> r = world_cartogram("population")
    >>> r.figure.savefig("world_population.png", dpi=150)  # doctest: +SKIP

    >>> r = world_cartogram({"USA": 1.0, "CHN": 1.2}, label="demo")  # doctest: +SKIP
    """
    from .world_data import WorldMap  # lazy — GIS stack optional

    wm = WorldMap.load(
        resolution=resolution,
        projection=projection,
        drop_antarctica=drop_antarctica,
    )
    if merge_sovereign:
        wm = wm.merge_by_sovereign()

    value_label = label or (metric if isinstance(metric, str) else "custom metric")

    # Resolve metric -> gdf["value"]
    gdf = wm.gdf
    values: np.ndarray
    if isinstance(metric, str):
        if metric not in gdf.columns:
            raise KeyError(
                f"metric {metric!r} is not a column of the world map; "
                f"available columns include: "
                f"{[c for c in gdf.columns if c != gdf.geometry.name]}"
            )
        values = gdf[metric].astype(float).to_numpy()
    elif isinstance(metric, Mapping):
        iso_col = wm.iso_column if wm.iso_column in gdf.columns else "SOV_A3"
        mapped = gdf[iso_col].map(lambda k: float(metric.get(k, 0.0)))
        values = mapped.to_numpy(dtype=float)
    elif callable(metric):
        values = np.array([float(metric(row)) for _, row in gdf.iterrows()])
    else:
        raise TypeError(
            "metric must be a column name, a dict, or a callable"
        )
    gdf = gdf.copy()
    gdf["value"] = values

    # Heavy-tail clipping — mirrors the Wolfram defaults for per-capita
    # metrics. Opt-in because absolute quantities rarely need it.
    if min_floor_fraction is not None:
        positive = values[values > 0]
        if positive.size:
            floor = float(min_floor_fraction) * float(positive.mean())
            gdf["value"] = np.maximum(gdf["value"].to_numpy(), floor)
            values = gdf["value"].to_numpy()

    # Choose a bbox. Default: the map's padded total bounds.
    if bbox is None:
        bbox = wm.padded_bbox(pad=0.02)

    # Rasterise polygons weighted by `value`.
    ny, nx = grid_size
    rho = rasterize_polygons(
        list(gdf.geometry),
        list(gdf["value"].astype(float)),
        bbox,
        (ny, nx),
    )

    if rho_ceil_multiplier is not None:
        pos = rho[rho > 0]
        if pos.size:
            ceiling = float(rho_ceil_multiplier) * float(pos.mean())
            rho = np.minimum(rho, ceiling)

    cart = Cartogram(
        rho, bbox=bbox,
        mean_floor=mean_floor,
        blur_sigma=blur_sigma,
        sea_density=sea_density,
    )
    cart.run(
        tol=tol,
        performance_goal=performance_goal,
        snapshots=snapshots,
    )
    deformed = cart.transform_polygons(list(gdf.geometry))

    fig = None
    if render:
        fig = _render_side_by_side(gdf, deformed, value_label, cmap=cmap)

    wm_out = WorldMap(
        gdf=gdf,
        value_column="value",
        name_column=wm.name_column,
        iso_column=wm.iso_column,
    )
    return WorldCartogramResult(
        cartogram=cart,
        world_map=wm_out,
        deformed_geometries=deformed,
        figure=fig,
    )


def _render_side_by_side(gdf, deformed, label: str, cmap: str = "plasma"):
    """Two-panel geographic + cartogram figure.

    Same log1p colour scale on both panels so the eye compares directly.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    vals = gdf["value"].to_numpy(dtype=float)
    finite = vals[np.isfinite(vals) & (vals > 0)]
    if finite.size == 0:
        norm = Normalize(vmin=0, vmax=1)
        color_vals = vals
    else:
        vmin = float(finite.min())
        vmax = float(finite.max())
        if vmax / max(vmin, 1e-12) > 50:
            norm = LogNorm(vmin=vmin, vmax=vmax)
        else:
            norm = Normalize(vmin=vmin, vmax=vmax)
        color_vals = np.where(vals > 0, vals, vmin)

    fig, (ax_geo, ax_cart) = plt.subplots(
        1, 2, figsize=(18, 6), constrained_layout=True
    )
    gdf.plot(
        column="value",
        ax=ax_geo, cmap=cmap, norm=norm,
        edgecolor="black", linewidth=0.2,
    )
    ax_geo.set_title("Geographic")
    ax_geo.set_aspect("equal")
    ax_geo.set_xticks([]); ax_geo.set_yticks([])

    # Build a temporary GeoDataFrame of the deformed geometries so we
    # can piggyback on geopandas' plot() + colourmap.
    import geopandas as gpd
    deformed_gdf = gpd.GeoDataFrame(
        {"value": vals},
        geometry=list(deformed),
        crs=gdf.crs,
    )
    deformed_gdf.plot(
        column="value",
        ax=ax_cart, cmap=cmap, norm=norm,
        edgecolor="black", linewidth=0.2,
    )
    ax_cart.set_title(f"Cartogram (area ∝ {label}, oceans preserved)")
    ax_cart.set_aspect("equal")
    ax_cart.set_xticks([]); ax_cart.set_yticks([])

    return fig
