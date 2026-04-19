#!/usr/bin/env wolframscript
(* Minimal test: use the synthetic 4x4 dense-block example that the
   regression suite pins, compare the deformed block area reached by
   Quality vs Speed. If the fast path is fundamentally broken, we'll
   see it here with no world-map confounds. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];
Get[FileNameJoin[{hereDir, "CartogramWLFast.wl"}]];

polyArea[p_] := 0.5 Abs[Total[MapThread[
  #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
  {p, RotateLeft[p]}]]];

rho = ConstantArray[1., {48, 48}];
rho[[18 ;; 30, 18 ;; 30]] = 9.0;
bbox = {0., 0., 1., 1.};

cart0 = CartogramWL`Cartogram[rho, bbox,
  "MeanFloor" -> 0.02, "BlurSigma" -> 1.0];

cartQ = CartogramWLFast`FastCartogramRun[cart0,
  PerformanceGoal -> "Quality", "Tol" -> 3.*^-3];
cartS = CartogramWLFast`FastCartogramRun[cart0,
  PerformanceGoal -> "Speed",   "Tol" -> 3.*^-3,
  "Snapshots" -> 60];

block = {{18/48., 18/48.}, {30/48., 18/48.},
         {30/48., 30/48.}, {18/48., 30/48.}};
movedQ = CartogramWL`CartogramTransform[cartQ, block];
movedS = CartogramWL`CartogramTransform[cartS, block];

a0 = polyArea[block];
aQ = polyArea[movedQ];
aS = polyArea[movedS];

Print["input block area        : ", a0];
Print["  Quality-mode block    : ", aQ, "  (x", aQ/a0, ")"];
Print["  Speed-mode  block    : ", aS, "  (x", aS/a0, ")"];
Print["  relative underflow    : ", 1 - aS/aQ, " (0 = no drift)"];

(* Also measure the full grid disagreement. *)
gridDiff = Abs[Flatten[
  cartQ["deformedGrid"] - cartS["deformedGrid"]]];
Print[""];
Print["grid point disagreement (unit square, extent 1):"];
Print["  max    = ", N@Max[gridDiff]];
Print["  mean   = ", N@Mean[gridDiff]];
Print["  median = ", N@Median[gridDiff]];
Print["  95%ile = ", N@Quantile[gridDiff, 0.95]];
