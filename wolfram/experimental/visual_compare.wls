#!/usr/bin/env wolframscript
(* Render the GDP cartogram twice — once with PerformanceGoal
   "Quality" (baseline solver) and once with "Speed" (snapshot-
   cached RHS) — and stack them so we can eyeball the difference. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];
Get[FileNameJoin[{hereDir, "CartogramWLFast.wl"}]];

Print["PerformanceGoal -> \"Quality\" ..."];
{tQ, figQ} = AbsoluteTiming[
  CartogramWLFast`FastWorldCartogram["GDP",
    PerformanceGoal -> "Quality"]];
Print["  ", tQ, " s"];

Print["PerformanceGoal -> \"Speed\"   ..."];
{tS, figS} = AbsoluteTiming[
  CartogramWLFast`FastWorldCartogram["GDP",
    PerformanceGoal -> "Speed"]];
Print["  ", tS, " s"];
Print["speedup: ", tQ/tS, "x"];

out = Column[{
  Labeled[figQ, Style[
    "PerformanceGoal -> \"Quality\"   (" <>
    ToString[NumberForm[tQ, 3]] <> " s)",
    14, Bold, Darker@Red], Top],
  Labeled[figS, Style[
    "PerformanceGoal -> \"Speed\"     (" <>
    ToString[NumberForm[tS, 3]] <> " s, " <>
    ToString[NumberForm[tQ/tS, 3]] <> "x)",
    14, Bold, Darker@Green], Top]},
  Spacings -> 1];

outPath = FileNameJoin[{hereDir, "visual_compare_gdp.png"}];
Export[outPath, out, ImageResolution -> 120];
Print["wrote ", outPath];
