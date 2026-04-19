#!/usr/bin/env wolframscript
(* Build wolfram/community/cartograms_update.nb — a short follow-up
   "comment" notebook to the main cartograms.nb post. Summarises the
   recent changes (Alaska/Hawaii handling, Speed default, Association
   input, missing-country display, satellite background, distortion
   grid), links back to the GitHub repository and the Wolfram section
   of the main README, and shows how to load the package and call it
   in two lines. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
imgDir  = FileNameJoin[{hereDir, "..", "docs", "images"}];

imgWidth = 680;

loadImg[name_] := With[{p = FileNameJoin[{imgDir, name}]},
  If[FileExistsQ[p],
    Image[Import[p], ImageSize -> imgWidth],
    Style["[missing: " <> name <> "]", Red]]];

imgCell[name_, caption_String: ""] := Cell[
  BoxData[ToBoxes[
    If[caption === "",
      loadImg[name],
      Column[{
        loadImg[name],
        Pane[Style[caption, Italic, FontSize -> 10,
                   FontFamily -> "Helvetica"],
             imgWidth, Alignment -> Center]},
        Alignment -> Center, Spacings -> 0.6]]]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{0, 0}, {10, 14}}];

(* URL constants — single source of truth for the links we reference. *)
repoUrl    = "https://github.com/mthiel74/Contiguous-Cartograms";
readmeWL   = repoUrl <> "#wolfram-language-implementation";
wolframDir = repoUrl <> "/tree/main/wolfram";
wlReadme   = repoUrl <> "/blob/main/wolfram/README.md";
wlFileView = repoUrl <> "/blob/main/wolfram/CartogramWL/CartogramWL.wl";
wlFileRaw  = "https://raw.githubusercontent.com/mthiel74/Contiguous-\
Cartograms/main/wolfram/CartogramWL/CartogramWL.wl";
mainNb     = repoUrl <>
             "/blob/main/wolfram/community/cartograms.nb";

link[text_, url_] := ButtonBox[text,
  BaseStyle -> "Hyperlink",
  ButtonData -> {URL[url], None}];

cells = {
  Cell["CartogramWL update \[LongDash] Alaska, Hawaii, and a faster \
default", "Title"],

  Cell["A short follow-up note on the diffusion-cartogram package \
introduced in the main post: what the recent commits changed, and \
where to download the current version.",
       "Subtitle"],

  Cell[TextData[{"Marco Thiel, April 2026 \[LongDash] ",
                 link["github.com/mthiel74/Contiguous-Cartograms",
                      repoUrl]}],
       "Text", TextAlignment -> Center, FontSlant -> Italic],

  (* ---- 1. TL;DR ---- *)

  Cell["1. What changed", "Section"],

  Cell[TextData[{"Since the original post the package ",
    StyleBox["CartogramWL", "Input"],
    " has picked up several improvements that are relevant if you \
tried the earlier version. The biggest single visible change is \
that ",
    StyleBox["WorldCartogram", "Input"],
    " now correctly handles multi-part countries \[LongDash] in \
particular, the United States is represented as the contiguous 48 \
states plus Alaska and Hawaii rather than just the lower 48."}],
       "Text"],

  Cell[TextData[{
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    "Alaska, Hawaii and other \"multi-part\" countries are now \
rendered and deformed correctly. ",
    StyleBox["CountryData[\"UnitedStates\", \"Polygon\"]", "Input"],
    " only returns the contiguous 48 states; a new helper ",
    StyleBox["CountryPolygonRings[c]", "Input"],
    " splices the Alaska and Hawaii ",
    StyleBox["AdministrativeDivision", "Input"],
    " polygons back in before rasterisation. The same hook can be \
used to bring in France's overseas departments, Portugal's \
Azores/Madeira, Denmark's Greenland, and so on.\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    "A new ",
    StyleBox["PerformanceGoal -> \"Speed\"", "Input"],
    " default pre-computes a handful of log-spaced density/gradient \
snapshots and linearly interpolates between them while the \
advection integrator runs. That cuts a full world cartogram from \
about 5 minutes to under a minute with sub-pixel agreement. ",
    StyleBox["PerformanceGoal -> \"Quality\"", "Input"],
    " restores the bit-exact baseline.\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    "Any ",
    StyleBox["ColorData[...]", "Input"],
    " gradient can be supplied as ",
    StyleBox["\"ColorFunction\"", "Input"],
    "; ",
    StyleBox["\"Background\" -> \"Satellite\"", "Input"],
    " textures each country with the corresponding patch of Wolfram's \
satellite imagery and deforms the imagery along with the cartogram.\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    "The metric argument is flexible: a recognised string, any \
CountryData property name, a custom ",
    StyleBox["Function[c, ...]", "Input"],
    ", or an ",
    StyleBox["Association", "Input"],
    " like ",
    StyleBox["<|\"Germany\" -> 83., \"France\" -> 67.5, ...|>", "Input"],
    " with missing countries optionally drawn in grey via ",
    StyleBox["\"MissingCountries\" -> \"ShowGrey\"", "Input"],
    ".\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    "A ",
    StyleBox["CartogramDistortionGrid[cart]", "Input"],
    " helper returns the warped-grid lines so you can visualise the \
deformation field directly."}],
       "Text"],

  (* ---- 2. Alaska story ---- *)

  Cell["2. Alaska and Hawaii: what this actually looks like",
       "Section"],

  Cell[TextData[{"The figure below shows the GDP cartogram of North \
America rendered twice with identical parameters. In the left \
column the US polygon comes straight from ",
    StyleBox["CountryData", "Input"],
    "; in the right column it comes from the new ",
    StyleBox["CountryPolygonRings", "Input"],
    ". Top row is the geographic input (US highlighted in orange), \
bottom row is the cartogram. With the fix the United States ends up \
with 28 mainland rings + 77 Alaska rings (mainland + Aleutians) + 8 \
Hawaiian rings, all sharing a single country value, and Alaska and \
Hawaii inflate alongside the mainland instead of being silently \
dropped."}], "Text"],

  imgCell["alaska_before_after.png",
    "Figure. Top row: geographic inputs. Bottom row: resulting \
cartogram. Left column uses bare CountryData[\"UnitedStates\", \
\"Polygon\"]; right column uses CartogramWL`CountryPolygonRings, \
which splices in the Alaska and Hawaii AdministrativeDivision \
polygons. The United States is highlighted in orange. Every ring of \
the same country shares the country's per-area density."],

  Cell[TextData[{"For reference, here is the full US ring set drawn \
in one go \[LongDash] 28 mainland rings, 77 Alaska rings, and 8 \
Hawaiian rings \[LongDash] so you can see what actually feeds the \
algorithm:"}], "Text"],

  imgCell["united_states_polygon.png",
    "Figure. CartogramWL`CountryPolygonRings[\"UnitedStates\"] returns \
113 exterior rings in a single list."],

  (* ---- 3. How to use it ---- *)

  Cell["3. Getting the current version", "Section"],

  Cell[TextData[{"The package is a single ",
    StyleBox[".wl", "Input"],
    " file with no external dependencies \[LongDash] pure Wolfram \
Language. The simplest way to try it right now is to download the \
file directly from GitHub and ",
    StyleBox["Get", "Input"],
    " it:"}], "Text"],

  Cell[BoxData[
"url  = \"" <> wlFileRaw <> "\";\n\
path = URLDownload[url];\n\
Get[path];\n\n\
(* Built-in metric *)\n\
WorldCartogram[\"Population\"]\n\n\
(* Custom colour scheme + satellite background *)\n\
WorldCartogram[\"GDP\",\n\
  \"ColorFunction\" -> ColorData[\"TemperatureMap\"],\n\
  \"Background\"    -> \"Satellite\"]\n\n\
(* Your own data as an Association *)\n\
WorldCartogram[\n\
  {\"my metric\", <|\"Germany\" -> 83., \"France\" -> 67.5, ...|>},\n\
  \"MissingCountries\" -> \"ShowGrey\"]"],
       "Input"],

  Cell[TextData[{"That is one call to ",
    StyleBox["URLDownload", "Input"],
    " followed by one call to ",
    StyleBox["Get", "Input"],
    "; no paclet registration, no git clone, nothing to install. \
After that all the public symbols (",
    StyleBox["WorldCartogram", "Input"],
    ", ",
    StyleBox["CartogramDistortionGrid", "Input"],
    ", ",
    StyleBox["CountryPolygonRings", "Input"],
    ", and the low-level primitives ",
    StyleBox["DiffusionSolver", "Input"],
    ", ",
    StyleBox["AdvectPoints", "Input"],
    ", ",
    StyleBox["Cartogram", "Input"],
    ", ",
    StyleBox["CartogramTransform", "Input"],
    ", ",
    StyleBox["RasterizePolygons", "Input"],
    ") live in the ",
    StyleBox["CartogramWL`", "Input"],
    " context. ",
    StyleBox["?WorldCartogram", "Input"],
    " prints the full option list."}], "Text"],

  Cell[TextData[{"Power users who want the examples, the test suite, \
the experimental benchmark sandbox, or the Python twin can clone the \
repository: ",
    link[repoUrl, repoUrl], "."}], "Text"],

  (* ---- 4. Where to read more ---- *)

  Cell["4. Documentation links", "Section"],

  Cell[TextData[{
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    link["Wolfram section of the repository README",
         readmeWL],
    " \[LongDash] headline figure, one-line call, full option table \
(ColorFunction, Background, MissingCountries, PerformanceGoal, \[Ellipsis]), \
and inline examples. Written for a general Wolfram user and probably \
the best single entry point if you have just arrived here.\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    link["wolfram/README.md", wlReadme],
    " \[LongDash] API cheat sheet, algorithm notes, grid conventions, \
and the mapping between the Wolfram and Python entry points.\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    link["wolfram/CartogramWL/CartogramWL.wl", wlFileView],
    " \[LongDash] the package source, ~700 lines of commented \
Wolfram Language. The direct raw-file link used above is: ",
    link[wlFileRaw, wlFileRaw], ".\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    link["wolfram/community/cartograms.nb", mainNb],
    " \[LongDash] the main walkthrough notebook that this note \
comments on.\n",
    StyleBox["\[Bullet] ", FontWeight -> Bold],
    link["the wolfram/ subtree on GitHub", wolframDir],
    " \[LongDash] examples (",
    StyleBox["world.wls", "Input"], ", ",
    StyleBox["world_metric.wls", "Input"], ", ",
    StyleBox["distortion_grid.wls", "Input"], ", ",
    StyleBox["alaska_before_after.wls", "Input"],
    "), the 7-case regression suite, and the experimental benchmark \
sandbox that records what was measured and what did or did not pay \
off."}], "Text"],

  Cell[TextData[{"Feedback, bug reports and PRs all welcome on \
GitHub. If you spot a country whose ",
    StyleBox["CountryData", "Input"],
    " polygon is incomplete in the same way the US was, the fix is \
one line in ",
    StyleBox["extraCountryRings", "Input"],
    " \[LongDash] happy to merge."}],
       "Text"]
};

nb = Notebook[cells,
  WindowSize -> {1100, 800},
  WindowMargins -> {{Automatic, 0}, {Automatic, 0}},
  FrontEndVersion -> "15.0 for Mac OS X ARM (64-bit) (January 2026)",
  StyleDefinitions -> "Default.nb",
  Saveable -> True];

outNb  = FileNameJoin[{hereDir, "cartograms_update.nb"}];
outPdf = FileNameJoin[{hereDir, "cartograms_update.pdf"}];

Print["writing notebook -> ", outNb];
Export[outNb, nb];

Print["writing pdf      -> ", outPdf];
Export[outPdf, nb];

Print["done; nb=", FileByteCount[outNb],
      "  pdf=", FileByteCount[outPdf]];
