#!/usr/bin/env wolframscript
(* Render a zoomed-in North America GDP cartogram before and after
   the ring-distribution fix, so we can see whether Alaska actually
   grows with the mainland now. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

(* Use the real WorldCartogram but restrict the bbox to North
   America so Alaska is clearly visible. *)
fig = CartogramWL`WorldCartogram["GDP",
   "BoundingBox" -> {-170., 10., -50., 75.},
   "Skip"        -> {"Antarctica", "FrenchSouthernTerritories",
                     "Bouvet Island"},
   "ImageSize"   -> 900];

outPath = FileNameJoin[{hereDir, "alaska_check.png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
