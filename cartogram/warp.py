"""Warp an arbitrary raster image through a :class:`Cartogram` deformation.

The deformation :math:`\\Phi` is known forward — the cartogram stores the
post-advection position of every original grid point. To warp an image
whose pixels live in original coordinates into one whose pixels live in
cartogram coordinates, we need the *inverse* deformation
:math:`\\Phi^{-1}`. We build it by scattered-data interpolation:
``(Φ(grid_i), grid_i)`` pairs become the training set, queried at the
target pixel positions.

This is exact for piecewise-linear images and good-to-subpixel-accuracy
for the deformations produced by the Gastner–Newman flow, which is
smooth away from the boundary.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np

from .cartogram import Cartogram


def warp_image(
    cart: Cartogram,
    image: np.ndarray,
    image_bbox: Tuple[float, float, float, float] | None = None,
    out_bbox: Tuple[float, float, float, float] | None = None,
    out_shape: Tuple[int, int] | None = None,
    fill: float | Tuple[float, ...] = 0.0,
) -> np.ndarray:
    """Warp ``image`` through the cartogram deformation.

    Parameters
    ----------
    cart : Cartogram
        Must have had :meth:`Cartogram.run` called on it.
    image : np.ndarray, shape (H, W) or (H, W, C)
        Source image, in the coordinate frame described by ``image_bbox``.
        Pixels are assumed to be equally spaced over ``image_bbox``; the
        pixel at ``image[i, j]`` corresponds to physical coordinate
        ``(xmin + (j+0.5) dx, ymin + (i+0.5) dy)``.
    image_bbox : (xmin, ymin, xmax, ymax), optional
        Extent of ``image`` in the same coordinate system as the
        cartogram. Defaults to the cartogram's own bbox.
    out_bbox : same, optional
        Extent of the output image. Defaults to the cartogram's bbox.
    out_shape : (H_out, W_out), optional
        Shape of the output image. Defaults to the input image's shape.
    fill : float or per-channel tuple
        Value used for output pixels whose pre-image is outside the
        source image (this can happen near the cartogram boundary).

    Returns
    -------
    np.ndarray of same dtype as ``image``
        The warped image, shape ``(H_out, W_out[, C])``.
    """
    if cart._final_grid_points is None:
        raise RuntimeError("cartogram has not been run() yet")

    image = np.asarray(image)
    if image.ndim not in (2, 3):
        raise ValueError("image must be 2-D or 3-D (H, W, C)")

    if image_bbox is None:
        image_bbox = cart.bbox
    if out_bbox is None:
        out_bbox = cart.bbox
    if out_shape is None:
        out_shape = image.shape[:2]

    H_out, W_out = out_shape
    oxmin, oymin, oxmax, oymax = out_bbox
    odx = (oxmax - oxmin) / W_out
    ody = (oymax - oymin) / H_out

    xs = oxmin + (np.arange(W_out) + 0.5) * odx
    ys = oymin + (np.arange(H_out) + 0.5) * ody
    X, Y = np.meshgrid(xs, ys)
    query = np.column_stack([X.ravel(), Y.ravel()])

    # Inverse deformation — query the scattered-point map
    # (cartogram_pt -> original_pt) at our target cartogram pixels.
    pre = _inverse_transform(cart, query)

    # Bilinear sampling from the source image.
    ixmin, iymin, ixmax, iymax = image_bbox
    H_in, W_in = image.shape[:2]
    idx = (pre[:, 0] - ixmin) / (ixmax - ixmin) * W_in - 0.5
    idy = (pre[:, 1] - iymin) / (iymax - iymin) * H_in - 0.5

    inside = (
        (idx >= 0) & (idx <= W_in - 1) & (idy >= 0) & (idy <= H_in - 1)
    )
    idx_c = np.clip(idx, 0, W_in - 1)
    idy_c = np.clip(idy, 0, H_in - 1)
    j0 = np.floor(idx_c).astype(np.intp)
    i0 = np.floor(idy_c).astype(np.intp)
    j1 = np.minimum(j0 + 1, W_in - 1)
    i1 = np.minimum(i0 + 1, H_in - 1)
    tx = idx_c - j0
    ty = idy_c - i0

    def _sample(channel: np.ndarray) -> np.ndarray:
        f00 = channel[i0, j0]
        f01 = channel[i0, j1]
        f10 = channel[i1, j0]
        f11 = channel[i1, j1]
        return (
            (1 - ty) * ((1 - tx) * f00 + tx * f01)
            + ty * ((1 - tx) * f10 + tx * f11)
        )

    if image.ndim == 2:
        vals = _sample(image.astype(np.float64))
        out = np.full(H_out * W_out, fill if np.ndim(fill) == 0 else fill[0], dtype=np.float64)
        out[inside] = vals[inside]
        return out.reshape(H_out, W_out).astype(image.dtype)
    else:
        C = image.shape[2]
        out = np.empty((H_out * W_out, C), dtype=np.float64)
        fill_vec = np.broadcast_to(np.asarray(fill, dtype=np.float64), (C,))
        out[:] = fill_vec
        for c in range(C):
            vals = _sample(image[..., c].astype(np.float64))
            out[inside, c] = vals[inside]
        return out.reshape(H_out, W_out, C).astype(image.dtype)


def _inverse_transform(cart: Cartogram, query: np.ndarray) -> np.ndarray:
    """Find original coordinates ``p`` such that ``Φ(p) ≈ query``.

    Uses ``scipy.interpolate.LinearNDInterpolator`` on the scattered
    ``(Φ(grid_i), grid_i)`` pairs. We fall back to the nearest-neighbour
    value for queries that land outside the deformed grid's convex hull.
    """
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator

    ny, nx = cart.solver.shape
    xs, ys = cart.solver.grid_coords()
    gx, gy = np.meshgrid(xs, ys)
    orig = np.column_stack([gx.ravel(), gy.ravel()])
    deformed = cart._final_grid_points  # shape (ny*nx, 2)

    lin_x = LinearNDInterpolator(deformed, orig[:, 0])
    lin_y = LinearNDInterpolator(deformed, orig[:, 1])
    near_x = NearestNDInterpolator(deformed, orig[:, 0])
    near_y = NearestNDInterpolator(deformed, orig[:, 1])

    px = lin_x(query)
    py = lin_y(query)
    bad = np.isnan(px) | np.isnan(py)
    if np.any(bad):
        px[bad] = near_x(query[bad])
        py[bad] = near_y(query[bad])
    return np.column_stack([px, py])
