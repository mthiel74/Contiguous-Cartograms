#!/usr/bin/env wolframscript
(* Visualise a GDP-driven cartogram as a grid deformation.

   Left panel  : countries (geographic) with a regular 15 degree lat/lon
                 grid overlaid.
   Right panel : the same countries and the same grid, both warped
                 through the cartogram.

   Output: wolfram/docs/images/distortion_grid.png *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

bbox = {-180., -60., 180., 85.};
{xmin, ymin, xmax, ymax} = bbox;

toRings[poly_Polygon] := Module[{inner, rings},
  inner = poly[[1]]; rings = First[inner];
  Map[{#[[2]], #[[1]]} &, rings, {2}]];

Print["collecting GDP for every country..."];
skip = {"Antarctica", "FrenchSouthernTerritories", "Bouvet Island"};
countries = Select[CountryData["Countries"],
  ! MemberQ[skip, CountryData[#, "Name"]] &];

entries = Reap[Do[
  With[{val  = CountryData[c, "GDP"],
        poly = CountryData[c, "Polygon"]},
    If[Head[val] === Quantity && Head[poly] === Polygon,
      Sow[<|"name"  -> CountryData[c, "Name"],
            "value" -> N@QuantityMagnitude[val],
            "rings" -> toRings[poly]|>]]],
  {c, countries}]][[2, 1]];
Print["  kept ", Length[entries], " countries"];

values  = entries[[All, "value"]];
floor   = 0.01 * Mean[Select[values, # > 0 &]];
entries = Map[Append[#, "value" -> Max[#["value"], floor]] &, entries];

(* Split every country into its constituent rings and distribute
   the GDP across rings by area so Alaska/Hawaii contribute to the
   density field alongside the mainland US. *)
{polyRaster, vals} = Transpose[Flatten[Map[
  Function[e, Module[{rs, areas, tot},
    rs    = e["rings"];
    areas = Map[0.5 Abs[Total[MapThread[
              #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
              {#, RotateLeft[#]}]]] &, rs];
    tot   = Total[areas];
    If[tot > 0,
      MapThread[{#1, e["value"] * #2/tot} &, {rs, areas}],
      {}]]],
  entries], 1]];
Print["rasterising..."];
rho     = CartogramWL`RasterizePolygons[polyRaster, vals, bbox,
            {192, 384}, 2];
posMean = Mean[Select[Flatten[rho], # > 0 &]];
rho     = Map[Min[#, 200. * posMean] &, rho, {2}];

Print["running cartogram..."];
cart = CartogramWL`Cartogram[rho, bbox,
  "MeanFloor" -> 0.02, "BlurSigma" -> 1.5, "SeaDensity" -> "auto"];
cart = CartogramWL`CartogramRun[cart, "Tol" -> 3.*^-3];

Print["warping polygons..."];
transformRing[ring_] :=
  CartogramWL`CartogramTransformPolygon[cart, ring, "MaxEdge" -> 2.0];
newEntries = Map[
  Append[#, "rings" -> Map[transformRing, #["rings"]]] &, entries];

Print["building distortion grid..."];
regularLines[] := Module[
  {nlines = 13, ncols = 25, xs, ys},
  xs = Table[xmin + j (xmax - xmin)/(ncols - 1), {j, 0, ncols - 1}];
  ys = Table[ymin + i (ymax - ymin)/(nlines - 1), {i, 0, nlines - 1}];
  Join[
    Table[Line[Table[{x, y}, {x, xs}]], {y, ys}],
    Table[Line[Table[{x, y}, {y, ys}]], {x, xs}]]];

gridLines     = regularLines[];
cartGridLines = CartogramWL`CartogramDistortionGrid[cart,
                  "Lines" -> 13, "Columns" -> 25];

Print["rendering..."];
countryStyleGeo  = {EdgeForm[Directive[Black, Thickness[0.0004]]],
                    FaceForm[GrayLevel[0.88]]};
countryStyleCart = {EdgeForm[Directive[Black, Thickness[0.0004]]],
                    FaceForm[Directive[RGBColor[1, 0.73, 0.35],
                                       Opacity[0.8]]]};
gridStyle = Directive[RGBColor[0.85, 0.18, 0.20], Thickness[0.0009]];

makeCountries[entryList_, style_] := Table[
  Module[{e = entryList[[k]]},
    {style, Polygon /@ e["rings"]}],
  {k, 1, Length[entryList]}];

leftPanel = Graphics[{
  makeCountries[entries, countryStyleGeo],
  {gridStyle, gridLines}},
  Frame -> True, FrameTicks -> None,
  PlotLabel -> "Geographic (lon/lat) + regular 15\[Degree] grid",
  AspectRatio -> Automatic, ImageSize -> 900];

rightPanel = Graphics[{
  makeCountries[newEntries, countryStyleCart],
  {gridStyle, cartGridLines}},
  Frame -> True, FrameTicks -> None,
  PlotLabel -> "Cartogram (area \[Proportional] GDP) + same grid, deformed",
  AspectRatio -> Automatic, ImageSize -> 900];

fig = GraphicsRow[{leftPanel, rightPanel},
  ImageSize -> 1800, Spacings -> 30,
  PlotLabel -> Style[
    "Cartogram distortion grid \[LongDash] world GDP",
    FontSize -> 14, Bold]];

outDir  = FileNameJoin[{hereDir, "..", "docs", "images"}];
outPath = FileNameJoin[{outDir, "distortion_grid.png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
