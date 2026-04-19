#!/usr/bin/env wolframscript
(* World-population cartogram using Wolfram's built-in CountryData.

   Output: wolfram/docs/images/world.png — side-by-side plain-geography
   and population-cartogram images.

   Runs in ~1 min on a modern laptop. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

(* ---- Convert a CountryData polygon into a list of {lon, lat} rings.
   The structure is Polygon[GeoPosition[{ring1, ring2, ...}]] where
   each ring is a list of {lat, lon} pairs. *)
toRings[poly_Polygon] := Module[{inner, rings},
  inner = poly[[1]];
  rings = First[inner];
  (* Swap to {lon, lat} for x-east, y-north. *)
  Map[{#[[2]], #[[1]]} &, rings, {2}]
];

(* ---- Pull every country that has both a polygon and a population.
   Skip Antarctica and micro-countries to avoid degeneracies. *)
skip = {"Antarctica", "FrenchSouthernTerritories", "Bouvet Island"};
countries = Select[CountryData["Countries"],
  ! MemberQ[skip, CountryData[#, "Name"]] &];

Print["pulling ", Length[countries], " countries..."];
entries = Reap[Do[
  With[{name = CountryData[c, "Name"],
        pop  = CountryData[c, "Population"],
        poly = CountryData[c, "Polygon"]},
    If[Head[pop] === Quantity && Head[poly] === Polygon,
      Sow[<|"name" -> name,
            "pop"  -> N@QuantityMagnitude[pop],
            "rings" -> toRings[poly]|>]
    ]
  ],
  {c, countries}]][[2, 1]];
Print["  kept ", Length[entries], " with polygon+population"];

(* ---- Cut bbox to the inhabited world (drop tiny far-Pacific slivers). *)
allPts = Flatten[entries[[All, "rings"]], 2];
{xmin, xmax} = {-180., 180.};
{ymin, ymax} = {-60., 85.};
bbox = {xmin, ymin, xmax, ymax};

(* ---- Rasterise: each country contributes population/area to its
   cells. "area" here is 2-D area in {lon, lat} space, which is not
   equal-area, but for the diffusion step it is enough. *)
Print["rasterising..."];
polyForRaster = entries[[All, "rings", 1]];   (* outer ring only *)
pops = entries[[All, "pop"]];
rho = CartogramWL`RasterizePolygons[polyForRaster, pops, bbox,
         {192, 384}, 2];
Print["  rho max=", Max[Flatten[rho]],
      " mean=", Mean[Flatten[rho]]];

Print["building cartogram (ocean preservation on)..."];
cart = CartogramWL`Cartogram[rho, bbox,
   "MeanFloor" -> 0.02, "BlurSigma" -> 1.5,
   "SeaDensity" -> "auto"];
cart = CartogramWL`CartogramRun[cart, "Tol" -> 3.*^-3];

Print["transforming polygons..."];
transformRing[ring_] := CartogramWL`CartogramTransformPolygon[cart,
                           ring, "MaxEdge" -> 2.0];
newEntries = Map[
  MapAt[Map[transformRing, #] &, #, "rings"] &, entries];

Print["rendering..."];
minP = Min[pops]; maxP = Max[pops];
colorFn[v_] := ColorData["SunsetColors"][
  (Log10[1. + v] - Log10[1. + minP])/
  (Log10[1. + maxP] - Log10[1. + minP])];

makeGraphics[entryList_, title_] := Graphics[
  Table[
    Module[{e = entryList[[k]], ringList, col},
      ringList = e["rings"];
      col = colorFn[e["pop"]];
      {EdgeForm[Directive[Black, Thickness[0.0005]]],
       FaceForm[col],
       Polygon /@ ringList}
    ],
    {k, 1, Length[entryList]}],
  Frame -> True, FrameTicks -> None,
  PlotLabel -> title, AspectRatio -> Automatic,
  ImageSize -> 900
];

fig = GraphicsRow[
  {makeGraphics[entries,    "Geographic (lon/lat)"],
   makeGraphics[newEntries,
    "Cartogram (area \[Proportional] population, oceans preserved)"]},
  ImageSize -> 1600, Spacings -> 30];

outDir = FileNameJoin[{hereDir, "..", "docs", "images"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir]];
outPath = FileNameJoin[{outDir, "world.png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
