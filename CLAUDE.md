# CLAUDE.md — project notes for Claude Code

This file tells Claude Code how to work on this repository.

## What this project is

A from-scratch Python implementation of the Gastner–Newman diffusion
cartogram algorithm. The scientific reference is:

> M. T. Gastner & M. E. J. Newman, "Diffusion-based method for producing
> density-equalizing maps", *PNAS* **101**(20), 7499–7504 (2004).

The goal is a small, auditable, pure-Python/NumPy implementation — not a
wrapper around the C `cart` or `go-cart` binaries.

## Directory conventions

* `cartogram/` — the importable package. Keep it dependency-light:
  `numpy`, `scipy`, `shapely`, `matplotlib`. Do **not** import
  `geopandas` at the top level of any core module — put it behind a
  local import inside `io_utils.py` so the core algorithm stays usable
  without the heavy GIS stack.
* `examples/` — runnable scripts. They should import from `cartogram`
  and write output (PNGs, new shapefiles) into `examples/output/`,
  which is git-ignored.
* `tests/` — `pytest`-style tests. They must run without network access
  and without any optional dependency.

## Algorithmic conventions

* Grids are stored as `numpy` arrays with shape `(ny, nx)` — row index
  `i` corresponds to `y`, column index `j` to `x`. Every function that
  takes a grid must document this.
* The physical bounding box is carried around as `bbox = (xmin, ymin,
  xmax, ymax)`. Conversions between physical coordinates and grid
  indices live in `cartogram.diffusion._grid_coords`.
* Diffusion is solved in Fourier space on a **padded, reflected** grid
  so that Neumann (no-flux) boundary conditions are enforced exactly.
  The padding factor defaults to 2. Do not change this default without
  updating the tests in `test_diffusion.py`.
* Time stepping for coordinate advection uses the Cash–Karp embedded
  RK4(5) pair from `scipy.integrate.solve_ivp` with
  `method="RK45"`, `rtol=1e-6`, `atol=1e-9`. The velocity field is
  evaluated by bilinear interpolation on the diffused density grid.

## Style rules

* Type-annotate all public functions.
* Prefer vectorised NumPy over Python loops; the one place a loop is
  acceptable is the outer time-stepping driver in `Cartogram.run`.
* No comments that just restate the code. Comments are reserved for
  *why* (non-obvious invariants, references to equation numbers in the
  paper).

## Testing

Run the full suite with

```
pytest -q
```

Tests must pass before any commit is pushed. The synthetic example

```
python examples/synthetic.py
```

should finish in well under 30 s on a laptop and produce
`examples/output/synthetic.png`.

## Git workflow

* One logical change per commit. Commit messages use the imperative
  mood ("Add diffusion solver", not "Added…").
* After each meaningful change, commit **and push** to `origin/main`.
* Never force-push; never rewrite published history.

## Things to avoid

* Do not add a GUI, a web UI, or a CLI wrapper unless explicitly asked.
* Do not vendor large example datasets into the repository. Real-world
  examples should fetch data on demand and cache under
  `examples/data/` (git-ignored).
* Do not introduce a new top-level dependency without updating both
  `requirements.txt` and this file.
