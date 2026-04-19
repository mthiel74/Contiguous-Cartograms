#!/usr/bin/env wolframscript
(* Just render the US polygon set (mainland + Alaska + Hawaii) as
   orange filled polygons on a North-America-sized frame so we can
   see unambiguously that Alaska is being drawn. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

rings = CartogramWL`CountryPolygonRings["UnitedStates"];
Print["US rings: ", Length[rings]];

fig = Graphics[{
  EdgeForm[Directive[Black, Thickness[0.0008]]],
  FaceForm[RGBColor[1, 0.55, 0.20]],
  Polygon /@ rings},
  Frame -> True,
  PlotRange -> {{-180., -60.}, {15., 75.}},
  AspectRatio -> Automatic,
  ImageSize -> 1200,
  Background -> RGBColor[0.92, 0.96, 1.0]];

outPath = FileNameJoin[{hereDir, "alaska_draw.png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
