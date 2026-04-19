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
    missing_countries: Literal["hide", "grey"] = "hide",
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
    missing_countries : {"hide", "grey"}
        How countries whose metric is missing / zero are handled on the
        rendered figure. ``"hide"`` (default) matches the Wolfram
        behaviour: those countries are dropped from both panels and
        contribute nothing to the density field. ``"grey"`` keeps them
        on the map in a neutral grey — they still do not deform the
        cartogram (their weight stays at the floor), but the reader
        can see the full world and tell "no data" apart from "zero".

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
        from .country_lookup import canonical, normalise_mapping

        canonical_map = normalise_mapping(metric)
        iso_col = wm.iso_column if wm.iso_column in gdf.columns else "SOV_A3"

        def _lookup(row_key: str) -> float:
            c = canonical(row_key)
            if c is None:
                c = row_key.strip().upper() if isinstance(row_key, str) else ""
            return float(canonical_map.get(c, 0.0))

        mapped = gdf[iso_col].map(_lookup)
        values = mapped.to_numpy(dtype=float)
    elif callable(metric):
        values = np.array([float(metric(row)) for _, row in gdf.iterrows()])
    else:
        raise TypeError(
            "metric must be a column name, a dict, or a callable"
        )
    gdf = gdf.copy()
    gdf["value"] = values

    # Track which rows are "missing" (NaN or <= 0 on a positive metric)
    # so we can honour missing_countries on the render side. We keep the
    # mask on the gdf so downstream consumers (plotting, joins) can see
    # it too.
    has_value_mask = np.isfinite(values) & (values > 0)
    gdf["__has_value"] = has_value_mask

    # Heavy-tail clipping — mirrors the Wolfram defaults for per-capita
    # metrics. Opt-in because absolute quantities rarely need it. Only
    # applied to rows that actually have a value; missing rows stay at
    # zero so they don't inflate the cartogram.
    if min_floor_fraction is not None:
        positive = values[has_value_mask]
        if positive.size:
            floor = float(min_floor_fraction) * float(positive.mean())
            clamped = np.where(
                has_value_mask,
                np.maximum(gdf["value"].to_numpy(), floor),
                gdf["value"].to_numpy(),
            )
            gdf["value"] = clamped
            values = clamped

    if missing_countries == "hide":
        # Drop missing rows from the gdf entirely; rasterisation /
        # polygon warp / plot all become naturally consistent.
        gdf = gdf.loc[has_value_mask].reset_index(drop=True)
        values = gdf["value"].to_numpy(dtype=float)
        has_value_mask = np.ones(len(gdf), dtype=bool)

    # Choose a bbox. Default: the map's padded total bounds.
    if bbox is None:
        bbox = wm.padded_bbox(pad=0.02)

    # Rasterise polygons weighted by `value`. Only rows that carry a
    # value contribute to the density field -- "grey" rows are drawn
    # on the plot but must not distort the cartogram.
    ny, nx = grid_size
    raster_geoms = [
        geom for geom, keep in zip(gdf.geometry, has_value_mask) if keep
    ]
    raster_values = [
        float(v) for v, keep in zip(gdf["value"], has_value_mask) if keep
    ]
    if not raster_geoms:
        raise ValueError(
            "no country carries a value for this metric; nothing to rasterise"
        )
    rho = rasterize_polygons(
        raster_geoms,
        raster_values,
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
        fig = _render_side_by_side(
            gdf, deformed, value_label, cmap=cmap,
            has_value_mask=has_value_mask,
        )

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


def _render_side_by_side(
    gdf, deformed, label: str,
    cmap: str = "plasma",
    has_value_mask=None,
):
    """Two-panel geographic + cartogram figure.

    Same log1p colour scale on both panels so the eye compares
    directly. If ``has_value_mask`` is supplied, rows where the mask is
    False are rendered in neutral grey on both panels (see the
    ``missing_countries="grey"`` option).
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, Normalize

    vals = gdf["value"].to_numpy(dtype=float)
    if has_value_mask is None:
        mask = np.ones(len(vals), dtype=bool)
    else:
        mask = np.asarray(has_value_mask, dtype=bool)

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

    grey = "#d9d9d9"
    valued_gdf = gdf.loc[mask]
    missing_gdf = gdf.loc[~mask]

    fig, (ax_geo, ax_cart) = plt.subplots(
        1, 2, figsize=(18, 6), constrained_layout=True
    )

    if not missing_gdf.empty:
        missing_gdf.plot(
            ax=ax_geo, color=grey,
            edgecolor="black", linewidth=0.2,
        )
    if not valued_gdf.empty:
        valued_gdf.plot(
            column="value", ax=ax_geo, cmap=cmap, norm=norm,
            edgecolor="black", linewidth=0.2,
        )
    ax_geo.set_title("Geographic")
    ax_geo.set_aspect("equal")
    ax_geo.set_xticks([]); ax_geo.set_yticks([])

    import geopandas as gpd
    deformed_list = list(deformed)
    if missing_gdf.empty:
        missing_deformed = gpd.GeoDataFrame(
            {"value": []}, geometry=[], crs=gdf.crs
        )
    else:
        missing_idx = np.where(~mask)[0]
        missing_deformed = gpd.GeoDataFrame(
            {"value": vals[~mask]},
            geometry=[deformed_list[i] for i in missing_idx],
            crs=gdf.crs,
        )
    valued_idx = np.where(mask)[0]
    valued_deformed = gpd.GeoDataFrame(
        {"value": vals[mask]},
        geometry=[deformed_list[i] for i in valued_idx],
        crs=gdf.crs,
    )

    if not missing_deformed.empty:
        missing_deformed.plot(
            ax=ax_cart, color=grey,
            edgecolor="black", linewidth=0.2,
        )
    if not valued_deformed.empty:
        valued_deformed.plot(
            column="value", ax=ax_cart, cmap=cmap, norm=norm,
            edgecolor="black", linewidth=0.2,
        )
    ax_cart.set_title(f"Cartogram (area ∝ {label}, oceans preserved)")
    ax_cart.set_aspect("equal")
    ax_cart.set_xticks([]); ax_cart.set_yticks([])

    return fig
