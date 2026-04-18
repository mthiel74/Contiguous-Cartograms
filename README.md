# Contiguous Cartograms

A Python implementation of the **Gastner–Newman diffusion-based cartogram algorithm**
for producing *contiguous* density-equalizing cartograms.

A cartogram is a map in which the area of each region is rescaled to be
proportional to some quantity of interest (population, GDP, number of cases,
electoral votes, …). A *contiguous* cartogram keeps neighbouring regions
attached — the topology is preserved — while the shapes are continuously
deformed.

## Mathematical background

We are given a density field `ρ(x, y)` on the plane (e.g. population per
unit area). We want a smooth, invertible deformation `Φ : R² → R²` such
that the pulled-back density `(ρ ∘ Φ⁻¹) · |det DΦ⁻¹|` is uniform.

Gastner & Newman (PNAS, 2004) observed that one can construct such a
deformation by *diffusing* the density until it is uniform, and letting
the material of the map be carried along with the induced flow. Concretely,
they solve the linear diffusion equation

```
∂ρ/∂t = ∇²ρ ,    ρ(x, 0) = ρ₀(x) ,
```

on a rectangle with Neumann (no-flux) boundary conditions. The associated
velocity field is the normalised flux

```
v(x, t) = − ∇ρ(x, t) / ρ(x, t) .
```

Every map point is then advected by this velocity field,

```
dr/dt = v(r(t), t) ,    r(0) = r₀ ,
```

and we take the limit `t → ∞`, in which `ρ` has relaxed to its spatial
mean and the flow stops. The resulting positions `r(∞)` are the cartogram
coordinates of the original points `r₀`.

Because the deformation is the flow of a smooth velocity field, it is a
homeomorphism: regions that were contiguous remain contiguous, and
(non-degenerate) boundaries stay non-self-intersecting.

The implementation here follows the spectral / FFT approach of the
original paper: the diffusion equation is solved exactly in Fourier space
on a padded periodic grid, and coordinates are integrated with an
adaptive Runge–Kutta scheme.

## Repository layout

```
cartogram/              core Python package
    __init__.py
    diffusion.py        FFT-based diffusion solver
    advect.py           adaptive RK integration of map points
    cartogram.py        high-level Cartogram class
    rasterize.py        polygon → density-grid rasteriser
    io_utils.py         optional GeoPandas I/O helpers
examples/
    synthetic.py        runs the algorithm on a toy 4×4 grid of cells
    us_states.py        (optional) runs on real US-state shapefiles
tests/
    test_diffusion.py   conservation + steady-state tests
    test_advect.py      flow-preservation tests
README.md
CLAUDE.md
requirements.txt
```

## Installation

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

`geopandas` is only required for the real-world example; the core
algorithm and synthetic demo depend solely on `numpy`, `scipy`,
`shapely` and `matplotlib`.

## Quick start

```python
import numpy as np
from cartogram import Cartogram

# A 256 × 256 density grid (e.g. rasterised population).
rho = np.ones((256, 256))
rho[64:192, 64:192] = 5.0   # a dense central block

cart = Cartogram(rho, bbox=(0, 0, 1, 1))
cart.run()

# Deform an arbitrary set of (x, y) points:
pts  = np.array([[0.25, 0.25], [0.75, 0.75]])
pts2 = cart.transform(pts)
```

See `examples/synthetic.py` for a complete runnable demo that produces a
before/after PNG of a 4×4 grid of square regions.

## References

* M. T. Gastner & M. E. J. Newman, *Diffusion-based method for producing
  density-equalizing maps*, PNAS **101** (20), 7499–7504 (2004).
* M. T. Gastner, V. Seguy & P. More, *Fast flow-based algorithm for
  creating density-equalizing map projections*, PNAS **115** (10),
  E2156–E2164 (2018).

## License

MIT — see `LICENSE`.
