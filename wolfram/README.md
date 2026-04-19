# CartogramWL — Wolfram-Language port

A from-scratch port of the Gastner–Newman diffusion cartogram
algorithm to the Wolfram Language. Mirrors the Python package in the
parent directory: same math, same grid convention, same behaviour.

## Demos

### Synthetic — 4×4 grid with one dense cell

![synthetic](docs/images/synthetic.png)

### World — population cartogram (oceans preserved)

Built from Wolfram's bundled `CountryData`.

![world](docs/images/world.png)

### Other built-in metrics

Driven by a single parameterised script (`examples/world_metric.wls`)
and the reusable `WorldCartogram[]` function (see below). Every image
on this page is produced from the same pipeline — only the metric,
the `"ColorFunction"` option, and the `"Background"` option change.

| Metric | Image |
|---|---|
| `GDP` (USD / year) | ![gdp](docs/images/world_gdp.png) |
| `GDPPerCapita` | ![gdppc](docs/images/world_gdp_per_capita.png) |
| `PopulationDensity` | ![popdens](docs/images/world_population_density.png) |

### Custom colour schemes

Pass any Wolfram `ColorData[…]` gradient as `"ColorFunction"`. Example
using `TemperatureMap` on the population-density cartogram:

![popdens-temp](docs/images/world_population_density_temperaturemap.png)

### `GeoBackground -> "Satellite"`

Set `"Background" -> "Satellite"` and each country is textured with
the corresponding patch of satellite imagery — on the cartogram side
the imagery is deformed to match the density-equalised geometry.
Example, the GDP cartogram with the satellite surface of the Earth:

![gdp-sat](docs/images/world_gdp_satellite.png)

## Requirements

* Wolfram Engine / Mathematica ≥ 13.0 (tested on WolframKernel 15.0).
* `wolframscript` on `$PATH` for the command-line runner.

Nothing else — the package uses only built-in symbols.

## Layout

```
wolfram/
    CartogramWL/
        CartogramWL.wl     # the package (algorithm + WorldCartogram)
    examples/
        synthetic.wls      # 4x4 synthetic demo
        world.wls          # world population cartogram
        world_metric.wls   # CLI wrapper over WorldCartogram[]
    tests/
        test_basic.wls     # regression tests
    docs/images/           # generated PNGs committed with the repo
```

## Running

```bash
cd wolfram

# Sanity tests (a few seconds):
wolframscript -file tests/test_basic.wls

# Synthetic demo (~5 s):
wolframscript -file examples/synthetic.wls

# World demo (~1 min, uses Wolfram's CountryData so no download):
wolframscript -file examples/world.wls

# Parameterised driver: <metric> [<colorScheme>] [<bg>]
wolframscript -file examples/world_metric.wls gdp
wolframscript -file examples/world_metric.wls gdp_per_capita
wolframscript -file examples/world_metric.wls population_density
wolframscript -file examples/world_metric.wls gdp SunsetColors satellite
wolframscript -file examples/world_metric.wls population_density TemperatureMap
```

### Calling from a notebook

```wolfram
Needs["CartogramWL`"];

(* One-liner — Gastner-Newman cartogram of every country's population. *)
fig = WorldCartogram["Population"];

(* Change the colour scheme. *)
WorldCartogram["GDP", "ColorFunction" -> ColorData["TemperatureMap"]]

(* Use satellite imagery instead of solid colours — each country is
   textured with the corresponding patch of GeoBackground -> "Satellite",
   deformed to match the cartogram. *)
WorldCartogram["GDP", "Background" -> "Satellite"]

(* Any CountryData property works too — e.g. life expectancy. *)
WorldCartogram["LifeExpectancy",
   "ColorFunction" -> ColorData["ThermometerColors"]]

(* Or supply a custom {label, fn} pair for derived quantities. *)
WorldCartogram[{"CO2 per capita",
   Function[c, Quantity[QuantityMagnitude[CountryData[c, "CO2Emissions"]] /
                        QuantityMagnitude[CountryData[c, "Population"]],
                        "Tonnes/Person/Year"]]}]

(* Feed your own data in as an Association — keys may be CountryData
   identifiers ("UnitedStates") or common names ("United States");
   countries missing from the Association are skipped. *)
myData = <|"Germany" -> 83.0, "France" -> 67.5, "Italy" -> 59.1,
           "Spain" -> 47.6, "Poland" -> 38.0|>;
WorldCartogram[{"my metric", myData},
   "BoundingBox" -> {-15., 30., 45., 70.}]

(* Missing-data handling: "Hide" (default) drops countries with no
   value; "ShowGrey" draws them in grey with no density contribution. *)
WorldCartogram[{"G20 population", g20Data},
   "MissingCountries" -> "ShowGrey"]
```

`?WorldCartogram` inside a notebook prints the full option list.

## API cheat sheet

```wolfram
Needs["CartogramWL`"];

(* 1. Build a cartogram from a density grid. *)
solver = DiffusionSolver[rho0, {xmin, ymin, xmax, ymax}];
rhoAtT = DensityAt[solver, t];

(* 2. Flow arbitrary points through the velocity field. *)
pts    = AdvectPoints[solver, {{x1, y1}, {x2, y2}, ...}];

(* 3. High-level driver. *)
cart = Cartogram[rho, bbox,
    "MeanFloor"   -> 0.02,
    "BlurSigma"   -> 1.5,
    "SeaDensity"  -> "auto"];     (* "auto" = preserve ocean area *)
cart = CartogramRun[cart, "Tol" -> 2.*^-3];

newPts   = CartogramTransform[cart, pts];
newPoly  = CartogramTransformPolygon[cart, polyCoords];

(* 4. Rasterise polygons + values to a density grid. *)
rho = RasterizePolygons[polyList, values, bbox, {ny, nx}, subpixel];
```

Conventions:

* Density arrays are `rho[[i, j]]` with row index `i ↦ y`, column
  index `j ↦ x` (same as NumPy).
* `bbox = {xmin, ymin, xmax, ymax}`.
* Polygons are lists of `{x, y}` pairs (single exterior ring per
  polygon — pass one ring at a time to `CartogramTransformPolygon`).

## Algorithm notes (same as the Python side)

1. Diffusion is solved analytically on the Neumann box via a 2-D
   DCT-II cosine expansion. Wolfram's `FourierDCT[·, 2]` /
   `FourierDCT[·, 3]` are already N-dimensional and orthonormal, so
   the forward/inverse pair is a single call each. The density at any
   time `t` is then

   ```
   rho(·, t) = idct2[ coeffs * Exp[-lambda * t] ]
   ```

   where `lambda[i, j] = pi^2 ((j/Lx)^2 + (i/Ly)^2)` — the standard
   Neumann eigenvalues of -Laplacian.

2. The velocity field `v = -grad(rho)/rho` is evaluated on the grid
   by central differences and interpolated bilinearly at advected
   points. Points are integrated through `v` with `NDSolveValue`
   using an explicit RK5 scheme.

3. The input density is pre-smoothed with `GaussianFilter` (default
   σ = 1 cell) — mathematically equivalent to starting diffusion at
   `t = σ² dx dy / 2` — so the polygon-rasterisation step-functions
   don't make the initial velocity field stiff.

4. `SeaDensity -> "auto"` fills ocean cells (cells with zero input
   density) with the mean of the land cells, so the spatial mean is
   uniform across the domain and only the *relative* sizes of
   countries change. Oceans stay put.

## Relationship to the Python package

The Wolfram port stays intentionally thin:

| Concern | Python | Wolfram |
|---|---|---|
| Diffusion solver | `cartogram/diffusion.py`       | `DiffusionSolver[…]`            |
| Advection        | `cartogram/advect.py` (SciPy)  | `AdvectPoints[…]` (`NDSolveValue`) |
| High-level class | `Cartogram`                    | `Cartogram[…]` association       |
| Polygon warp     | `transform_polygons`           | `CartogramTransformPolygon`      |
| Rasteriser       | `rasterize_polygons`           | `RasterizePolygons[…]`           |
| Ocean preservation | `sea_density="auto"`         | `"SeaDensity" -> "auto"`         |

The Python package is the more feature-complete of the two — it has
the Gradio GUI, the live World Bank/WHO fetchers, the Natural Earth
shaded-relief warper, and a CLI. The Wolfram port is designed for
users who'd rather stay inside Mathematica and reach the core
algorithm + built-in `CountryData` in a few lines.

## License

MIT — same as the parent repository.
