#!/usr/bin/env wolframscript
(* Basic regression tests for CartogramWL. Exits 0 on success, 1 on
   failure. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

$failures = 0;
check[name_, cond_] := If[TrueQ[cond],
  Print["PASS: ", name],
  (Print["FAIL: ", name]; $failures++)];

(* 1. Round-trip: DensityAt[solver, 0] reproduces the input density. *)
Module[{rho, s, r0},
  SeedRandom[11];
  rho = Table[1. + 0.3 RandomReal[], {16}, {24}] + 2.0;
  s = CartogramWL`DiffusionSolver[rho, {0., 0., 2.0, 1.0}];
  r0 = CartogramWL`DensityAt[s, 0.0];
  check["DensityAt[t=0] == rho0", Max[Abs[rho - r0]] < 1.*^-10]];

(* 2. Mean conservation under diffusion. *)
Module[{rho, s, m0, m1, m2},
  SeedRandom[12];
  rho = Table[1. + 0.5 RandomReal[], {20}, {20}] + 2.0;
  s = CartogramWL`DiffusionSolver[rho, {0., 0., 1.0, 1.0}];
  m0 = Mean[Flatten[rho]];
  m1 = Mean[Flatten@CartogramWL`DensityAt[s, 0.1]];
  m2 = Mean[Flatten@CartogramWL`DensityAt[s, 10.0]];
  check["mean conserved at t=0.1", Abs[m1 - m0] < 1.*^-10];
  check["mean conserved at t=10",  Abs[m2 - m0] < 1.*^-10]];

(* 3. Density converges to uniform for large t. *)
Module[{rho, s, rhoInf},
  SeedRandom[13];
  rho = Table[1. + 0.5 Sin[0.5 i] Cos[0.3 j], {i, 1, 24}, {j, 1, 24}] + 2.0;
  s = CartogramWL`DiffusionSolver[rho, {0., 0., 1.0, 1.0}];
  rhoInf = CartogramWL`DensityAt[s, CartogramWL`ConvergenceTime[s, 1.*^-8]];
  check["density relaxes to uniform for large t",
    StandardDeviation[Flatten[rhoInf]] < 1.*^-4]];

(* 4. Uniform density yields identity flow. *)
Module[{rho, s, pts, out, err},
  rho = ConstantArray[1., {16, 20}];
  s = CartogramWL`DiffusionSolver[rho, {0., 0., 2.0, 1.0}];
  pts = {{0.5, 0.5}, {1.0, 0.5}, {0.1, 0.9}};
  out = CartogramWL`AdvectPoints[s, pts, "TMax" -> 3.0];
  err = Max[Abs[Flatten[pts - out]]];
  check["uniform density => identity flow (err < 1e-6)", err < 1.*^-6]];

(* 5. Dense block expands under Cartogram. *)
Module[{rho, cart, block, moved, area0, area1},
  rho = ConstantArray[1., {48, 48}];
  rho[[18 ;; 30, 18 ;; 30]] = 9.0;
  cart = CartogramWL`Cartogram[rho, {0., 0., 1.0, 1.0},
     "MeanFloor" -> 0.02, "BlurSigma" -> 1.0];
  cart = CartogramWL`CartogramRun[cart, "Tol" -> 3.*^-3];
  block = {{18/48., 18/48.}, {30/48., 18/48.},
           {30/48., 30/48.}, {18/48., 30/48.}};
  moved = CartogramWL`CartogramTransform[cart, block];
  area0 = 0.5 Abs[Total[MapThread[
    #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
    {block, RotateLeft[block]}]]];
  area1 = 0.5 Abs[Total[MapThread[
    #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
    {moved, RotateLeft[moved]}]]];
  check["dense block at least doubles in area", area1 > 2 area0]];

(* 6. Corners stay on corners under the cartogram deformation. *)
Module[{rho, cart, corners, moved, err},
  rho = ConstantArray[1., {32, 32}];
  rho[[10 ;; 22, 10 ;; 22]] = 5.0;
  cart = CartogramWL`Cartogram[rho, {0., 0., 1.0, 1.0},
    "MeanFloor" -> 0.02, "BlurSigma" -> 1.0];
  cart = CartogramWL`CartogramRun[cart, "Tol" -> 3.*^-3];
  corners = {{0., 0.}, {1., 0.}, {1., 1.}, {0., 1.}};
  moved = CartogramWL`CartogramTransform[cart, corners];
  err = Max[Abs[Flatten[corners - moved]]];
  check["corners stay within 0.03 of original", err < 0.03]];

Print["=== ", If[$failures == 0, "ALL PASS", $failures, " FAILURE(S)"]];
Quit[If[$failures == 0, 0, 1]]
