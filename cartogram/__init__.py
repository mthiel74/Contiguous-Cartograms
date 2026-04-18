"""Diffusion-based contiguous cartograms (Gastner & Newman 2004).

Public surface:

    DiffusionSolver    — analytic diffusion on a rectangle with Neumann BCs
    Cartogram          — high-level driver tying diffusion to point advection
    rasterize_polygons — polygon list -> density grid
    advect_points      — integrate points through the diffusion-induced flow
"""

from .diffusion import DiffusionSolver
from .advect import advect_points
from .cartogram import Cartogram
from .rasterize import rasterize_polygons

__all__ = [
    "DiffusionSolver",
    "Cartogram",
    "advect_points",
    "rasterize_polygons",
    "WorldMap",
    "warp_image",
]


def __getattr__(name: str):
    # Lazy-load the optional GeoPandas-dependent surface so the core
    # package stays importable in environments without the GIS stack.
    if name == "WorldMap":
        from .world_data import WorldMap
        return WorldMap
    if name == "warp_image":
        from .warp import warp_image
        return warp_image
    raise AttributeError(f"module 'cartogram' has no attribute {name!r}")

__version__ = "0.1.0"
