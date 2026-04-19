# experimental/ — improvement ideas for CartogramWL

This directory is a sandbox. Nothing in the public `wolfram/CartogramWL`
package is modified. We benchmark the existing implementation,
prototype alternatives, and measure the wins (or losses) honestly.

## Baseline profile — where the time actually goes

From `profile_baseline.wls`, one full world-GDP cartogram at the
default `(192, 384)` grid:

| stage                  | seconds | share |
|---|---:|---:|
| `CountryData` pull     |   5.6 |   2% |
| `RasterizePolygons`    |   9.6 |   3% |
| `Cartogram` build      |   0.1 |  <1% |
| **`CartogramRun`**     | **320.0** | **94%** |
| polygon warp (233×)    |   9.5 |   3% |
| **total**              | **345** | |

So `CartogramRun` is the only stage that matters for wall-clock.
Everything else could be made free and the user wouldn't notice.

## The win: snapshot-cached advection RHS

`CartogramRun` integrates the advection ODE

    dr/dt = -∇ρ(r, t) / ρ(r, t)

with every grid cell centre (192 × 384 ≈ 74 k points) packed into a
single NDSolveValue state vector. The expensive part is evaluating
the RHS — at every RK45 sub-step `DensityAndGradientAt[solver, t]`
does a fresh inverse DCT of the full 192 × 384 spectrum and a
finite-difference gradient on the result.

`CartogramWLFast.wl` keeps the NDSolveValue adaptive stepping (a
fixed-step RK4 blows up in stiff high-density regions — we tried
and the result diverged by 99 % from the baseline) but pre-computes
a table of `{rho(t_k), ∂x ρ(t_k), ∂y ρ(t_k)}` at N log-spaced time
samples, and linearly interpolates in t during advection. Each RHS
call now costs three array blends + three bilinear samples — no
inverse DCTs. The gradient uses array shifts instead of `Do` loops.

### Measured speed-up (world GDP, macOS Apple Silicon, Wolfram 15)

| nSnaps | `CartogramRun` time | speed-up | max grid err | mean | median |
|---:|---:|---:|---:|---:|---:|
| baseline  | 316 s | 1.0×  |  — | — | — |
| 25   |  55 s | 5.8×  | 21.8° | 2.46° | 1.40° |
| 50   |  31 s | 10.3× | 16.7° | 1.35° | 0.55° |
| 100  |  29 s | 10.9× | 13.1° | 0.66° | 0.16° |

Errors are in degrees of longitude/latitude on a 360 × 145° frame.
At rendering resolution (900 px wide) one pixel is ~0.4°, so median
error at 50 snapshots is already subpixel.

Going from 50 to 100 snapshots barely changes wall-clock (the
bottleneck shifts to NDSolveValue's own stepping logic) but halves
the worst-point error. The fast defaults ship with `"Snapshots" -> 60`.

### Pipeline totals

| stage | baseline | fast (50 snaps) |
|---|---:|---:|
| rasterise      |  10 s |  10 s |
| cartogram run  | 316 s |  31 s |
| polygon warp   |  10 s |   6 s |
| **total**      | **336 s** | **47 s** |
| speed-up       | 1.0× | **7.2×** |

## Things we tried that didn't pay off

### `FastRasterizePolygons`
Replacing the `RegionMember`-based subpixel loop with a
`Rasterize[Graphics[Polygon[…]]]` pipeline sounds obvious, but:

* One `Rasterize` call per country has enough fixed overhead to be
  roughly as slow as the baseline's subpixel loop, and
* The anti-aliased edge pixels disagree with the baseline's hard
  subpixel integration by ~20 % on cells the polygon partially
  covers. It's "different wrong", not "better".

Rasterisation is 3 % of wall-clock. Not worth fighting.

### Custom fixed-step RK4 over the grid
Ditching NDSolveValue gave a 127× speed-up on `CartogramRun` … and
disagreed with the baseline by 357° out of 360°. The velocity field
near dense countries is genuinely stiff; without adaptive step
control the integrator overshoots and points leave the domain. The
adaptive Cash–Karp step control inside NDSolveValue is doing real
work; we need to either keep it or implement our own PI controller.

### `FastCartogramTransformPolygon`
Vectorising the edge-densification `AppendTo` loop into a single
`Join` + `Table` is a genuine but small win: 1.5× on the 10 s
polygon-warp stage, i.e. ~3 s saved end-to-end.

## How to use

```wolfram
Get["/path/to/wolfram/CartogramWL/CartogramWL.wl"];
Get["/path/to/wolfram/experimental/CartogramWLFast.wl"];

fig = CartogramWLFast`FastWorldCartogram["GDP"];
Export["world_gdp.png", fig]
```

Same options as the baseline `WorldCartogram` plus `"Snapshots"` on
`FastCartogramRun` for tuning accuracy/speed.

## Ideas left on the table

Speed (biggest → smallest expected win):

* **Advect a coarser grid, bilinearly interpolate back.** The
  cartogram deformation is smooth; 96 × 192 points (4× fewer)
  recover the full field under interpolation. Could take the 31 s
  fast run down to ~10 s, at which point we're I/O-bound.
* **Compile the bilinear sampler.** `Compile[…, CompilationTarget ->
  "C"]` on the inner point-sampling loop. `NDSolveValue` currently
  calls a top-level Wolfram function 100+ times; compiling it would
  drop the per-call overhead.
* **Interpolate the DCT coefficients in log-space** between
  snapshots instead of rho itself — mathematically exact for
  separable modes, since `e^{-λt_1} ^(1-α) · e^{-λt_2}^α = e^{-λt}`.
  Eliminates the 16° worst-case error, at the cost of an
  element-wise multiplication + exp per snapshot lookup (still much
  cheaper than a fresh inverse DCT).
* **Parallel polygon warp.** `ParallelMap` over the 233 countries
  works today; the startup cost on `LaunchKernels[]` eats most of
  the win on a single invocation but is free in a batch pipeline
  (e.g. the Gradio demo UI on the Python side).
* **GPU `FourierDCT`** via `CUDAFourier` / `OpenCLLink` — only
  relevant if we bump the grid resolution past 512 × 1024.

Correctness / quality:

* **Antimeridian splitting.** Russia, Fiji, Chukotka wrap around
  ±180°. The baseline clips; a proper split before rasterising would
  make the cartogram geometry exactly correct near the dateline.
* **Adaptive snapshot spacing** driven by `||ρ(t_k) − ρ(t_{k−1})||`
  rather than fixed log-spacing — could drop to 20 snapshots with
  no accuracy loss.
* **Sparse rasteriser.** For continent-scale cartograms most cells
  stay zero. Storing rho sparse would save memory (though not time,
  since the diffusion step forces it dense).

Ergonomics:

* **Native-math notebook builder.** The published
  `community/cartograms.nb` uses `$…$` placeholders that render as
  MathJax on Wolfram Community but show as literal text in the
  notebook's own PDF. A version of `build_notebook.wls` that emits
  equations as Wolfram box expressions (`FractionBox`, `SubscriptBox`,
  `FormBox`) would render everywhere. Straightforward, just
  tedious.
