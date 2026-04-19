#!/usr/bin/env wolframscript
(* Parameterised world cartogram — thin CLI wrapper over
   WorldCartogram[].

   Usage:
       wolframscript -file world_metric.wls <metric> [<colorScheme>] [<bg>]

   metric      : population | gdp | gdp_per_capita | population_density
                 (or any CountryData property name).
   colorScheme : any Wolfram ColorData gradient name (default SunsetColors)
                 — e.g. TemperatureMap, ThermometerColors, AvocadoColors.
   bg          : flat (default) | satellite
                 — "satellite" textures each country with the corresponding
                 patch of GeoBackground -> "Satellite" imagery instead of a
                 solid colour fill.

   Output: wolfram/docs/images/world_<metric>[_<scheme>][_<bg>].png *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

metric = If[Length[$ScriptCommandLine] >= 2,
            $ScriptCommandLine[[2]], "population"];
colorScheme = If[Length[$ScriptCommandLine] >= 3,
                 $ScriptCommandLine[[3]], "SunsetColors"];
bgMode = If[Length[$ScriptCommandLine] >= 4,
            ToLowerCase[$ScriptCommandLine[[4]]], "flat"];
If[! MemberQ[{"flat", "satellite"}, bgMode],
   Print["unknown bg mode ", bgMode, " (expected flat|satellite)"];
   Quit[2]];

Print["metric: ", metric, "  scheme: ", colorScheme, "  bg: ", bgMode];

fig = CartogramWL`WorldCartogram[metric,
   "ColorFunction" -> ColorData[colorScheme],
   "Background"    -> bgMode];
If[fig === $Failed, Quit[3]];

outDir = FileNameJoin[{hereDir, "..", "docs", "images"}];
If[! DirectoryQ[outDir], CreateDirectory[outDir]];
outStem = "world_" <> ToLowerCase[metric];
If[colorScheme =!= "SunsetColors",
   outStem = outStem <> "_" <> ToLowerCase[colorScheme]];
If[bgMode =!= "flat",
   outStem = outStem <> "_" <> bgMode];
outPath = FileNameJoin[{outDir, outStem <> ".png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
