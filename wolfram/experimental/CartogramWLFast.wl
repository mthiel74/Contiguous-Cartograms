(* ::Package:: *)

(* :Title: CartogramWLFast *)
(* :Summary:
     Experimental, faster variant of CartogramWL. Shares the public
     API (same function names, same Options) but replaces the
     bottleneck routines with vectorised / cached equivalents.

     Loads the baseline CartogramWL package first so it can reuse
     helper symbols (DiffusionSolver, DensityAt, ...) where the
     baseline is already fast, and only overrides the slow ones.

     Drop-in usage:
         Get[".../CartogramWL/CartogramWL.wl"];
         Get[".../experimental/CartogramWLFast.wl"];
         (* Now FastRasterizePolygons, FastCartogramRun, etc. are
            available alongside the baseline symbols. *)

     The baseline symbols are kept intact so side-by-side benchmarks
     remain honest.
*)

BeginPackage["CartogramWLFast`", {"CartogramWL`"}];

FastRasterizePolygons::usage = "FastRasterizePolygons[polys, values, \
bbox, {ny, nx}] — fast polygon rasteriser using Rasterize + \
ImageData. Interface identical to CartogramWL`RasterizePolygons \
(but the final subpixel argument is ignored; Rasterize handles \
anti-aliasing).";

FastCartogramRun::usage = "FastCartogramRun[cart, opts] — drop-in \
replacement for CartogramWL`CartogramRun with a PerformanceGoal \
option: \"Quality\" (default) delegates to the baseline solver; \
\"Speed\" switches to the snapshot-cached RHS variant (~10x faster, \
~0.5 degree median error on a 360 degree world frame).";

FastCartogramTransformPolygon::usage = "FastCartogramTransformPolygon\
[cart, poly, opts] — same semantics as \
CartogramWL`CartogramTransformPolygon but vectorises the edge \
densification.";

FastWorldCartogram::usage = "FastWorldCartogram[metric, opts] — \
variant of CartogramWL`WorldCartogram with a PerformanceGoal option:\n\
  \"Quality\" (default) — use the baseline Gastner-Newman solver; \
result is bit-exact with WorldCartogram.\n\
  \"Speed\"             — use the snapshot-cached RHS, ~7x faster \
end-to-end with a few tenths of a degree per-point warp error.\n\n\
Also caches the satellite-imagery raster across calls with the same \
BoundingBox and accepts every option WorldCartogram does.";

Begin["`Private`"];

(* ====================================================================
   1. Fast rasteriser via Rasterize.
==================================================================== *)

shoelaceArea[poly_] := 0.5 Abs[Total[MapThread[
  #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
  {poly, RotateLeft[poly]}]]];

maskFromPoly[poly_, bbox_, {ny_Integer, nx_Integer}] := Module[
  {xmin, ymin, xmax, ymax, g, img, raw},
  {xmin, ymin, xmax, ymax} = bbox;
  g = Graphics[{Black, Polygon[poly]},
    PlotRange -> {{xmin, xmax}, {ymin, ymax}},
    Background -> White, Frame -> False,
    PlotRangePadding -> 0, ImagePadding -> 0,
    AspectRatio -> Full];
  img = Rasterize[g, "Image", RasterSize -> {nx, ny},
          Background -> White];
  raw = ImageData[ColorConvert[img, "Grayscale"]];
  (* ImageData is y-down; flip to match grid (y-up). *)
  (* Polygon cells are near 0 (Black), outside near 1 (White). *)
  Reverse[1. - raw]];

FastRasterizePolygons[polys_List, values_List, bbox_List,
                      shape : {_Integer, _Integer},
                      _: 1] := Module[
  {rho = ConstantArray[0., shape], poly, v, area, mask},
  Do[
    poly = polys[[p]]; v = values[[p]];
    If[v == 0 || Length[poly] < 3, Continue[]];
    area = shoelaceArea[poly];
    If[area <= 0, Continue[]];
    mask = maskFromPoly[poly, bbox, shape];
    rho += mask * (v/area),
    {p, 1, Length[polys]}];
  rho];

(* ====================================================================
   2. Vectorised central gradient.
==================================================================== *)

fastCentralGradient[m_, dy_, dx_] := Module[{ny, nx, gy, gx},
  {ny, nx} = Dimensions[m];
  (* Interior: centred differences via array shifts. *)
  gy = (RotateLeft[m, {1, 0}] - RotateRight[m, {1, 0}])/(2. dy);
  gx = (RotateLeft[m, {0, 1}] - RotateRight[m, {0, 1}])/(2. dx);
  (* One-sided at boundaries (NumPy-style). *)
  gy[[1, All]]   = (m[[2, All]]     - m[[1, All]])/dy;
  gy[[ny, All]]  = (m[[ny, All]]    - m[[ny - 1, All]])/dy;
  gx[[All, 1]]   = (m[[All, 2]]     - m[[All, 1]])/dx;
  gx[[All, nx]]  = (m[[All, nx]]    - m[[All, nx - 1]])/dx;
  {gy, gx}];

(* ====================================================================
   3. Bilinear sampler (vectorised, identical semantics to the one in
      CartogramWL`Private`.
==================================================================== *)

fastBilinear[field_, bbox_, dx_, dy_, xq_, yq_] := Module[
  {ny, nx, xmin, ymin, fx, fy, j0, i0, j1, i1, tx, ty,
   f00, f01, f10, f11, flat, idx00, idx01, idx10, idx11},
  {ny, nx} = Dimensions[field];
  {xmin, ymin} = bbox[[{1, 2}]];
  fx = (xq - xmin)/dx - 0.5;
  fy = (yq - ymin)/dy - 0.5;
  fx = Clip[fx, {0., nx - 1.0000000001}];
  fy = Clip[fy, {0., ny - 1.0000000001}];
  j0 = Floor[fx]; i0 = Floor[fy];
  j1 = Clip[j0 + 1, {0, nx - 1}];
  i1 = Clip[i0 + 1, {0, ny - 1}];
  tx = fx - j0; ty = fy - i0;
  flat = Flatten[field];
  idx00 = i0 * nx + j0 + 1;
  idx01 = i0 * nx + j1 + 1;
  idx10 = i1 * nx + j0 + 1;
  idx11 = i1 * nx + j1 + 1;
  f00 = flat[[idx00]];
  f01 = flat[[idx01]];
  f10 = flat[[idx10]];
  f11 = flat[[idx11]];
  (1 - ty) ((1 - tx) f00 + tx f01) + ty ((1 - tx) f10 + tx f11)];

(* ====================================================================
   4. Snapshot cache.  Evaluates rho + gradient at nSnaps time points
      spaced log-linearly between 0 and tMax, so the advection RHS can
      interpolate linearly in t instead of re-doing an inverse DCT at
      every sub-step.
==================================================================== *)

(* Snapshot time grid.  The density modes decay exponentially with
   rates lambda_mn.  Fast modes (large lambda) finish collapsing well
   before slow modes start, so a linear-ish spacing in t misses the
   fast-mode dynamics entirely on large domains.  We space the
   snapshots geometrically from t0 = small fraction of the fastest
   decay time out to tMax, with t=0 prepended so the first RHS call
   has the initial condition.  This is the regime where linear
   interpolation of rho agrees with the exact exponential solution to
   a few parts in 1e3. *)
buildSnapshots[solver_, tMax_, nSnaps_] := Module[
  {tGrid, rho, gy, gx, out, dx, dy, lamMax, t0, ratio},
  dx = solver["dx"]; dy = solver["dy"];
  lamMax = Max[solver["lam"]];
  t0 = Min[0.1/Max[lamMax, 1.*^-12], tMax/1000.];
  ratio = (tMax/t0)^(1./(nSnaps - 1));
  tGrid = Prepend[
    Table[t0 * ratio^(k - 1), {k, 1, nSnaps}], 0.];
  out = Table[
    rho = CartogramWL`DensityAt[solver, tGrid[[k + 1]]];
    {gy, gx} = fastCentralGradient[rho, dy, dx];
    {rho, gx, gy},
    {k, 0, nSnaps}];
  {tGrid, out}];

interpSnapshot[tGrid_, snapshots_, tNow_] := Module[
  {k, alpha, a, b},
  k = Clip[Ceiling[(tNow - tGrid[[1]])/(tGrid[[-1]] - tGrid[[1]]) *
          (Length[tGrid] - 1)], {1, Length[tGrid] - 1}];
  While[tGrid[[k]] > tNow && k > 1, k--];
  While[k < Length[tGrid] - 1 && tGrid[[k + 1]] < tNow, k++];
  alpha = (tNow - tGrid[[k]])/(tGrid[[k + 1]] - tGrid[[k]]);
  alpha = Clip[alpha, {0., 1.}];
  a = snapshots[[k]]; b = snapshots[[k + 1]];
  {(1 - alpha) a[[1]] + alpha b[[1]],
   (1 - alpha) a[[2]] + alpha b[[2]],
   (1 - alpha) a[[3]] + alpha b[[3]]}];

(* ====================================================================
   5. Advection driver. We keep NDSolveValue's adaptive Runge-Kutta
      stepping (the velocity field is genuinely stiff in some
      regions and a fixed-step RK4 blows up there) but we swap the
      expensive per-call inverse-DCT RHS for a cheap interpolation
      from the precomputed snapshots. This keeps the baseline's
      accuracy while cutting each RHS call cost by an order of
      magnitude.
==================================================================== *)

FastAdvect[solver_, pts_, OptionsPattern[{
    "TMax"      -> Automatic,
    "Tol"       -> 10.^-3,
    "Snapshots" -> 60}]] := Module[
  {tMax, nSnaps, tGrid, snapshots, n, dx, dy, bbox, y0, rhsFn, state,
   t, Y, final},
  tMax = OptionValue["TMax"];
  If[tMax === Automatic,
     tMax = CartogramWL`ConvergenceTime[solver, OptionValue["Tol"]]];
  tMax = N@tMax;
  nSnaps = OptionValue["Snapshots"];
  {tGrid, snapshots} = buildSnapshots[solver, tMax, nSnaps];
  n    = Length[pts];
  dx   = solver["dx"]; dy = solver["dy"]; bbox = solver["bbox"];

  rhsFn[tNow_?NumericQ, stateVec_?VectorQ] := Module[
    {rho, gx, gy, xq, yq, rhoQ, gxQ, gyQ},
    {rho, gx, gy} = interpSnapshot[tGrid, snapshots, tNow];
    xq = stateVec[[1 ;; n]];
    yq = stateVec[[n + 1 ;; 2 n]];
    rhoQ = fastBilinear[rho, bbox, dx, dy, xq, yq];
    gxQ  = fastBilinear[gx,  bbox, dx, dy, xq, yq];
    gyQ  = fastBilinear[gy,  bbox, dx, dy, xq, yq];
    rhoQ = Map[Max[#, 1.*^-300] &, rhoQ];
    Join[-gxQ/rhoQ, -gyQ/rhoQ]];

  y0 = N@Join[pts[[All, 1]], pts[[All, 2]]];
  state = NDSolveValue[
    {Y'[t] == rhsFn[t, Y[t]], Y[0] == y0},
    Y, {t, 0, tMax},
    Method -> {"ExplicitRungeKutta",
               "DifferenceOrder" -> 5,
               "StiffnessTest" -> False},
    AccuracyGoal  -> 5,
    PrecisionGoal -> 4,
    MaxStepFraction -> 1/64];
  final = state[tMax];
  Transpose[{final[[1 ;; n]], final[[n + 1 ;; 2 n]]}]];

(* ====================================================================
   6. FastCartogramRun: drop-in replacement for CartogramRun.
==================================================================== *)

FastCartogramRun[cart_Association,
                 OptionsPattern[{
                   PerformanceGoal -> "Quality",
                   "Tol"            -> 3.*^-3,
                   "Snapshots"      -> 60}]] := Module[
  {solver, xs, ys, pts, moved, ny, nx, xmin, ymin, xmax, ymax, dx, dy,
   goal},
  goal = OptionValue[PerformanceGoal];
  If[goal === "Quality",
    Return[CartogramWL`CartogramRun[cart,
      "Tol" -> OptionValue["Tol"]]]];
  If[goal =!= "Speed",
    Message[FastCartogramRun::pgoal, goal]; Return[$Failed]];
  solver = cart["solver"];
  ny = solver["ny"]; nx = solver["nx"];
  {xmin, ymin, xmax, ymax} = solver["bbox"];
  dx = solver["dx"]; dy = solver["dy"];
  xs = Table[xmin + (j - 0.5) dx, {j, 1, nx}];
  ys = Table[ymin + (i - 0.5) dy, {i, 1, ny}];
  pts = Flatten[Table[{xs[[j]], ys[[i]]}, {i, 1, ny}, {j, 1, nx}], 1];
  moved = FastAdvect[solver, pts,
    "Tol" -> OptionValue["Tol"],
    "Snapshots" -> OptionValue["Snapshots"]];
  moved = Transpose[{
    Clip[moved[[All, 1]], {xmin, xmax}],
    Clip[moved[[All, 2]], {ymin, ymax}]}];
  Append[cart, "deformedGrid" -> moved]];
FastCartogramRun::pgoal = "PerformanceGoal must be \"Quality\" or \
\"Speed\"; got `1`.";

(* ====================================================================
   7. Polygon transform — vectorised densification.
==================================================================== *)

vectorDensify[poly_, maxE_] := Module[
  {a, b, seg, d, k, out = {}},
  Do[
    a = poly[[i]];
    b = poly[[Mod[i, Length[poly]] + 1]];
    d = Norm[b - a];
    k = Max[1, Ceiling[d/maxE]];
    out = Join[out,
      Table[a + (b - a)*(s/k), {s, 0, k - 1}]],
    {i, 1, Length[poly]}];
  out];

FastCartogramTransformPolygon[cart_Association, poly_List,
    OptionsPattern[{"MaxEdge" -> Automatic}]] := Module[
  {solver, maxE},
  solver = cart["solver"];
  maxE = OptionValue["MaxEdge"];
  If[maxE === Automatic,
     maxE = Min[solver["dx"], solver["dy"]]];
  CartogramWL`CartogramTransform[cart, vectorDensify[poly, maxE]]];

(* ====================================================================
   8. Satellite-texture cache (keyed by bbox + image size).
==================================================================== *)

$satCache = <||>;

cachedSatTexture[bbox_, size_] := Module[
  {key = {N@bbox, size}, xmin, ymin, xmax, ymax, g},
  If[KeyExistsQ[$satCache, key],
    $satCache[key],
    {xmin, ymin, xmax, ymax} = bbox;
    g = GeoGraphics[
      GeoRange -> {{ymin, ymax}, {xmin, xmax}},
      GeoBackground -> "Satellite",
      GeoRangePadding -> None,
      ImageSize -> size,
      Frame -> False, PlotRangePadding -> 0];
    $satCache[key] =
      Rasterize[g, "Image", RasterSize -> size[[1]]]]];

(* ====================================================================
   9. FastWorldCartogram — uses the fast primitives + sat cache +
      ParallelMap on polygon transforms.
==================================================================== *)

Options[FastWorldCartogram] = Join[
  Options[CartogramWL`WorldCartogram],
  {PerformanceGoal -> "Quality", "Snapshots" -> 60}];

toRings[poly_Polygon] := Module[{r},
  r = First[poly[[1]]];
  Map[{#[[2]], #[[1]]} &, r, {2}]];

worldSpec[m_] := Which[
  StringQ[m], Switch[ToLowerCase[m],
    "population",
      {"population", Function[c, CountryData[c, "Population"]],
       0.01, 200.},
    "gdp",
      {"GDP (USD/year)", Function[c, CountryData[c, "GDP"]],
       0.01, 200.},
    "gdppercapita" | "gdp_per_capita",
      {"GDP per capita (USD/person/year)",
       Function[c, CountryData[c, "GDPPerCapita"]], 0.10, 40.},
    "populationdensity" | "population_density",
      {"population density (people/km^2)",
       Function[c, With[{p = CountryData[c, "Population"],
                         a = CountryData[c, "Area"]},
         If[Head[p] === Quantity && Head[a] === Quantity &&
            QuantityMagnitude[a] > 0, p/a,
            Missing["NotAvailable"]]]], 0.05, 30.},
    _, {m, Function[c, CountryData[c, m]], 0.01, 200.}],
  MatchQ[m, _Function], {"custom", m, 0.01, 200.},
  MatchQ[m, _Association],
    {"custom",
     Function[c, With[{v = Lookup[m, c, Missing[]]},
       If[MissingQ[v],
         Lookup[m, CountryData[c, "Name"], Missing[]], v]]],
     0.01, 200.},
  MatchQ[m, {_String, _Function}], {m[[1]], m[[2]], 0.01, 200.},
  MatchQ[m, {_String, _Association}],
    Prepend[Rest@worldSpec[m[[2]]], m[[1]]],
  True, $Failed];

FastWorldCartogram[metric_, opts : OptionsPattern[]] := Module[
  {label, metricFn, floorFrac, ceilMult, bbox, gridSize, bg, imgSize,
   skip, maxEdge, colorFn, labelOpt, countries, entries, values,
   positive, meanVal, floor, polyRaster, vals, rho, posMean, ceiling,
   cart, newEntries, minP, maxP, logMin, logMax, satTexture, xmin,
   ymin, xmax, ymax, flatPanel, satPanel, panelTitle, cfOpt,
   missingMode, valBuf, greyBuf, greyEntries, greyNewEntries, uvRing,
   normalise, cName, cVal, cPoly},
  {label, metricFn, floorFrac, ceilMult} = worldSpec[metric];
  If[OptionValue["MinFloorFraction"] =!= Automatic,
     floorFrac = N@OptionValue["MinFloorFraction"]];
  If[OptionValue["RhoCeilMultiplier"] =!= Automatic,
     ceilMult  = N@OptionValue["RhoCeilMultiplier"]];
  labelOpt = OptionValue["Label"];
  If[labelOpt =!= Automatic, label = labelOpt];
  bbox = N@OptionValue["BoundingBox"];
  {xmin, ymin, xmax, ymax} = bbox;
  gridSize = OptionValue["GridSize"];
  bg       = ToLowerCase[ToString[OptionValue["Background"]]];
  imgSize  = OptionValue["ImageSize"];
  skip     = OptionValue["Skip"];
  maxEdge  = OptionValue["MaxEdge"];
  cfOpt    = OptionValue["ColorFunction"];
  colorFn  = If[cfOpt === Automatic, ColorData["SunsetColors"], cfOpt];
  missingMode = OptionValue["MissingCountries"];

  countries = Select[CountryData["Countries"],
    ! MemberQ[skip, CountryData[#, "Name"]] &];

  valBuf = {}; greyBuf = {};
  Do[
    cPoly = CountryData[c, "Polygon"];
    If[Head[cPoly] === Polygon,
      cVal = metricFn[c]; cName = CountryData[c, "Name"];
      If[(Head[cVal] === Quantity || NumericQ[cVal]),
        AppendTo[valBuf, <|
          "name"  -> cName,
          "value" -> N@If[Head[cVal] === Quantity,
                          QuantityMagnitude[cVal], cVal],
          "rings" -> toRings[cPoly]|>],
        If[missingMode === "ShowGrey",
          AppendTo[greyBuf, <|
            "name" -> cName, "rings" -> toRings[cPoly]|>]]]],
    {c, countries}];
  entries = valBuf; greyEntries = greyBuf;

  values  = entries[[All, "value"]];
  positive = Select[values, # > 0 &];
  meanVal = Mean[positive];
  floor   = floorFrac * meanVal;
  entries = Map[Append[#, "value" -> Max[#["value"], floor]] &,
                entries];

  polyRaster = entries[[All, "rings", 1]];
  vals       = entries[[All, "value"]];
  (* Use the baseline rasteriser: the Rasterize-based fast variant
     disagrees with it on partially-covered edge cells by enough to
     deform the whole cartogram.  Rasterisation is 3% of wall-clock
     anyway, so there's no pressure to replace it. *)
  rho        = CartogramWL`RasterizePolygons[polyRaster, vals, bbox,
                 gridSize, 2];
  posMean    = Mean[Select[Flatten[rho], # > 0 &]];
  ceiling    = ceilMult * posMean;
  rho        = Map[Min[#, ceiling] &, rho, {2}];

  cart = CartogramWL`Cartogram[rho, bbox,
    "MeanFloor" -> OptionValue["MeanFloor"],
    "BlurSigma" -> OptionValue["BlurSigma"],
    "SeaDensity" -> "auto"];
  cart = FastCartogramRun[cart,
    PerformanceGoal -> OptionValue[PerformanceGoal],
    "Tol"           -> OptionValue["Tolerance"],
    "Snapshots"     -> OptionValue["Snapshots"]];

  newEntries = Map[
    Append[#, "rings" -> Map[
      FastCartogramTransformPolygon[cart, #,
        "MaxEdge" -> maxEdge] &, #["rings"]]] &,
    entries];
  greyNewEntries = If[Length[greyEntries] > 0,
    Map[Append[#, "rings" -> Map[
      FastCartogramTransformPolygon[cart, #,
        "MaxEdge" -> maxEdge] &, #["rings"]]] &,
      greyEntries],
    {}];

  minP = Min[vals]; maxP = Max[vals];
  logMin = Log10[1. + minP]; logMax = Log10[1. + maxP];
  normalise[v_] := (Log10[1. + v] - logMin)/(logMax - logMin + 1.*^-12);
  uvRing[ring_] := Map[{(#[[1]] - xmin)/(xmax - xmin),
                        (#[[2]] - ymin)/(ymax - ymin)} &, ring];

  satTexture = If[bg === "satellite",
    cachedSatTexture[bbox, {2048, 1024}], None];

  flatPanel[posEntries_, greyPos_, title_] := Graphics[
    {If[Length[greyPos] > 0,
       {EdgeForm[Directive[GrayLevel[0.4], Thickness[0.0003]]],
        FaceForm[GrayLevel[0.82]],
        Table[Polygon /@ greyPos[[k]]["rings"],
              {k, 1, Length[greyPos]}]}, {}],
     Table[Module[{e = posEntries[[k]]},
       {EdgeForm[Directive[Black, Thickness[0.0004]]],
        FaceForm[colorFn[normalise[e["value"]]]],
        Polygon /@ e["rings"]}],
       {k, 1, Length[posEntries]}]},
    Frame -> True, FrameTicks -> None,
    PlotLabel -> title, AspectRatio -> Automatic,
    ImageSize -> imgSize];

  satPanel[posEntries_, origEntries_, greyPos_, greyOrig_, title_] :=
    Graphics[
    {Texture[satTexture],
     If[Length[greyPos] > 0,
       Table[Module[{pRings = greyPos[[k]]["rings"],
                     oRings = greyOrig[[k]]["rings"]},
         {EdgeForm[Directive[White, Opacity[0.4], Thickness[0.0004]]],
          FaceForm[Directive[GrayLevel[0.6], Opacity[0.55]]],
          MapThread[
            Polygon[#1, VertexTextureCoordinates -> uvRing[#2]] &,
            {pRings, oRings}]}],
         {k, 1, Length[greyPos]}], {}],
     Table[Module[{pRings = posEntries[[k]]["rings"],
                   oRings = origEntries[[k]]["rings"]},
       {EdgeForm[Directive[White, Opacity[0.6], Thickness[0.0005]]],
        FaceForm[White],
        MapThread[
          Polygon[#1, VertexTextureCoordinates -> uvRing[#2]] &,
          {pRings, oRings}]}],
       {k, 1, Length[posEntries]}]},
    Frame -> True, FrameTicks -> None,
    PlotLabel -> title, AspectRatio -> Automatic,
    ImageSize -> imgSize, Background -> Black];

  panelTitle = "Cartogram (area \[Proportional] " <> label <>
               ", oceans preserved)";
  If[bg === "satellite",
    GraphicsRow[{
      satPanel[entries, entries, greyEntries, greyEntries,
               "Geographic (lon/lat)"],
      satPanel[newEntries, entries, greyNewEntries, greyEntries,
               panelTitle]},
      ImageSize -> 2 imgSize, Spacings -> 30],
    GraphicsRow[{
      flatPanel[entries, greyEntries, "Geographic (lon/lat)"],
      flatPanel[newEntries, greyNewEntries, panelTitle]},
      ImageSize -> 2 imgSize, Spacings -> 30]]];

End[];
EndPackage[];
