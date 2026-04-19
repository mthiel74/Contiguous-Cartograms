#!/usr/bin/env wolframscript
(* Profile the three heavy stages of the baseline world cartogram
   pipeline: RasterizePolygons, CartogramRun (diffusion + advection),
   and the per-country polygon warp. Prints absolute timings so we
   know what's worth optimising. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

bbox = {-180., -60., 180., 85.};
skip = {"Antarctica", "FrenchSouthernTerritories", "Bouvet Island"};

toRings[poly_Polygon] := Module[{r},
  r = First[poly[[1]]];
  Map[{#[[2]], #[[1]]} &, r, {2}]];

Print["[t=0] collecting country polygons + GDP..."];
{tPull, entries} = AbsoluteTiming[
  Reap[Do[
    With[{val = CountryData[c, "GDP"],
          poly = CountryData[c, "Polygon"]},
      If[Head[val] === Quantity && Head[poly] === Polygon,
        Sow[<|"value" -> N@QuantityMagnitude[val],
              "rings" -> toRings[poly]|>]]],
    {c, Select[CountryData["Countries"],
        ! MemberQ[skip, CountryData[#, "Name"]] &]}]][[2, 1]]];
Print["  CountryData pull   : ", tPull, " s, ", Length[entries],
      " countries"];

values  = entries[[All, "value"]];
floor   = 0.01 * Mean[Select[values, # > 0 &]];
entries = Map[Append[#, "value" -> Max[#["value"], floor]] &, entries];

polys = entries[[All, "rings", 1]];
vals  = entries[[All, "value"]];

{tRaster, rho} = AbsoluteTiming[
  CartogramWL`RasterizePolygons[polys, vals, bbox, {192, 384}, 2]];
Print["  RasterizePolygons  : ", tRaster, " s"];

posMean = Mean[Select[Flatten[rho], # > 0 &]];
rho = Map[Min[#, 200. posMean] &, rho, {2}];

{tCart, cart} = AbsoluteTiming[
  CartogramWL`Cartogram[rho, bbox,
    "MeanFloor" -> 0.02, "BlurSigma" -> 1.5,
    "SeaDensity" -> "auto"]];
Print["  Cartogram build    : ", tCart, " s"];

{tRun, cart2} = AbsoluteTiming[
  CartogramWL`CartogramRun[cart, "Tol" -> 3.*^-3]];
Print["  CartogramRun       : ", tRun, " s"];

{tWarp, newEntries} = AbsoluteTiming[
  Map[Append[#, "rings" -> Map[
    CartogramWL`CartogramTransformPolygon[cart2, #,
      "MaxEdge" -> 2.0] &, #["rings"]]] &,
    entries]];
Print["  Polygon warp       : ", tWarp, " s"];

Print[""];
Print["TOTAL (ex CountryData pull) = ",
      tRaster + tCart + tRun + tWarp, " s"];
