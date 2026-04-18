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
]

__version__ = "0.1.0"
