#!/usr/bin/env wolframscript
(* Synthetic 4x4-grid demo, analogous to examples/synthetic.py.

   Produces wolfram/docs/images/synthetic.png — a side-by-side
   "original / cartogram" plot of a 4x4 grid of unit squares where
   one central cell has ~10x the population of the rest. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
PrependTo[$Path, FileNameJoin[{hereDir, "..", "CartogramWL"}]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

ncols = 4; nrows = 4;
bbox = {0., 0., N@ncols, N@nrows};

(* Hot cell at (i,j) = (2,2) 0-indexed. *)
hot = {2, 2};
SeedRandom[0];
polys = Flatten[Table[
  {{{j - 1., i - 1.}, {j, i - 1.}, {j, i}, {j - 1., i}}},
  {i, 1, nrows}, {j, 1, ncols}], 2];
pops = Flatten@Table[
  If[{i, j} - 1 === hot, 12., 1. + 0.3 RandomReal[]],
  {i, 1, nrows}, {j, 1, ncols}];

Print["rasterizing..."];
rho = CartogramWL`RasterizePolygons[polys, pops, bbox, {192, 192}, 3];
Print["rho mean: ", Mean[Flatten[rho]],
      ", rho max: ", Max[Flatten[rho]]];

Print["building cartogram..."];
cart = CartogramWL`Cartogram[rho, bbox,
   "MeanFloor" -> 0.02, "BlurSigma" -> 1.5];
cart = CartogramWL`CartogramRun[cart, "Tol" -> 2.*^-3];

Print["transforming polygons..."];
newPolys = CartogramWL`CartogramTransformPolygon[cart, #] & /@ polys;

Print["rendering..."];
minPop = Min[pops]; maxPop = Max[pops];
colorFn[v_] := ColorData["SunsetColors"][(v - minPop)/(maxPop - minPop + 0.001)];
makeGraphics[ps_, title_] := Graphics[
  MapThread[{EdgeForm[Black], FaceForm[colorFn[#2]], Polygon[#1]} &,
            {ps, pops}],
  Frame -> True, FrameTicks -> None,
  PlotLabel -> title, AspectRatio -> Automatic,
  PlotRangePadding -> Scaled[0.02],
  PlotRange -> {{-0.3, ncols + 0.3}, {-0.3, nrows + 0.3}}
];
fig = GraphicsRow[
  {makeGraphics[polys, "Original"],
   makeGraphics[newPolys, "Cartogram (area \[Proportional] population)"]},
  Spacings -> 30,
  PlotLabel -> "Diffusion cartogram \[LongDash] 4x4 grid with one dense cell",
  ImageSize -> 900
];

outDir = FileNameJoin[{hereDir, "..", "docs", "images"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir]];
outPath = FileNameJoin[{outDir, "synthetic.png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
