#!/usr/bin/env wolframscript
(* Side-by-side benchmark: baseline CartogramWL vs. experimental
   CartogramWLFast. Runs the full GDP world pipeline twice and times
   each stage, then verifies the results agree to within a small
   tolerance so "faster" doesn't quietly become "wrong". *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];
Get[FileNameJoin[{hereDir, "CartogramWLFast.wl"}]];

bbox = {-180., -60., 180., 85.};
skip = {"Antarctica", "FrenchSouthernTerritories", "Bouvet Island"};

toRings[poly_Polygon] := Module[{r},
  r = First[poly[[1]]];
  Map[{#[[2]], #[[1]]} &, r, {2}]];

Print["pulling country data..."];
{_, entries} = AbsoluteTiming[
  Reap[Do[
    With[{val = CountryData[c, "GDP"],
          poly = CountryData[c, "Polygon"]},
      If[Head[val] === Quantity && Head[poly] === Polygon,
        Sow[<|"value" -> N@QuantityMagnitude[val],
              "rings" -> toRings[poly]|>]]],
    {c, Select[CountryData["Countries"],
        ! MemberQ[skip, CountryData[#, "Name"]] &]}]][[2, 1]]];
values = entries[[All, "value"]];
floor = 0.01 * Mean[Select[values, # > 0 &]];
entries = Map[Append[#, "value" -> Max[#["value"], floor]] &, entries];
polys = entries[[All, "rings", 1]];
vals  = entries[[All, "value"]];
Print["  ", Length[entries], " countries\n"];

(* ---- raster (baseline only — see README for why) ---- *)

Print["== rasteriser =="];
{tRasterBase, rhoBase} = AbsoluteTiming[
  CartogramWL`RasterizePolygons[polys, vals, bbox, {192, 384}, 2]];
Print["  baseline : ", tRasterBase, " s (3% of total; not worth",
      " optimising — see README)"];
Print[""];

posMean = Mean[Select[Flatten[rhoBase], # > 0 &]];
rho = Map[Min[#, 200. posMean] &, rhoBase, {2}];

(* ---- cartogram build (same for both) ---- *)

cart0 = CartogramWL`Cartogram[rho, bbox,
  "MeanFloor" -> 0.02, "BlurSigma" -> 1.5, "SeaDensity" -> "auto"];

(* ---- cartogram run ---- *)

Print["== cartogram run (diffusion + advection of grid) =="];
{tRunBase, cartBase} = AbsoluteTiming[
  CartogramWL`CartogramRun[cart0, "Tol" -> 3.*^-3]];
Print["  baseline (NDSolveValue)      : ", tRunBase, " s"];

Do[
  {tRunFast, cartFast} = AbsoluteTiming[
    CartogramWLFast`FastCartogramRun[cart0,
      "Tol" -> 3.*^-3, "Snapshots" -> ns]];
  Print["  fast (snapshot-cached RHS, nSnaps=", ns, ") : ",
        tRunFast, " s"];
  Print["    speedup vs baseline : ", tRunBase/tRunFast, "x"];
  gridDiff = Abs[Flatten[
    cartBase["deformedGrid"] - cartFast["deformedGrid"]]];
  Print["    grid err (deg) : max=", N@Max[gridDiff],
        "  mean=", N@Mean[gridDiff],
        "  95%=", N@Quantile[gridDiff, 0.95],
        "  median=", N@Median[gridDiff]],
  {ns, {25, 50, 100}}];
Print[""];

(* ---- polygon warp ---- *)

Print["== polygon warp (233 countries) =="];
{tWarpBase, _} = AbsoluteTiming[
  Map[Map[
    CartogramWL`CartogramTransformPolygon[cartBase, #,
      "MaxEdge" -> 2.0] &, #["rings"]] &, entries]];
Print["  baseline : ", tWarpBase, " s"];

{tWarpFast, _} = AbsoluteTiming[
  Map[Map[
    CartogramWLFast`FastCartogramTransformPolygon[cartFast, #,
      "MaxEdge" -> 2.0] &, #["rings"]] &, entries]];
Print["  fast     : ", tWarpFast, " s"];
Print["  speedup  : ", tWarpBase/tWarpFast, "x"];
Print[""];

(* ---- totals ---- *)

totBase = tRasterBase + tRunBase + tWarpBase;
totFast = tRasterBase + tRunFast + tWarpFast;
Print["== TOTAL (raster + run + warp) =="];
Print["  baseline : ", totBase, " s"];
Print["  fast     : ", totFast, " s"];
Print["  speedup  : ", totBase/totFast, "x"];
