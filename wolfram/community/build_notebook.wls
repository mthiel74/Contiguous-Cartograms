#!/usr/bin/env wolframscript
(* Build wolfram/community/cartograms.nb — a Wolfram Notebook suitable
   for posting to the Wolfram Community, with native math typesetting
   (no literal $...$ placeholders). Also exports the same content as
   a PDF for offline layout review. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
imgDir  = FileNameJoin[{hereDir, "..", "docs", "images"}];

(* Target image width (in Wolfram display units). Chosen so that the
   rasterised content bleeds to the edges of the printable area of
   the default PDF page (Letter, standard margins) without
   overflowing. *)
imgWidth = 680;

loadImg[name_] := With[
  {p = FileNameJoin[{imgDir, name}]},
  If[FileExistsQ[p],
    Image[Import[p], ImageSize -> imgWidth],
    Style["[missing image: " <> name <> "]", Red]]];

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

(* -----------------------------------------------------------------
   Math typesetting helpers. Every equation is built as a Box
   expression and wrapped in FormBox[..., TraditionalForm] so Wolfram
   draws it with real math fonts (fractions, subscripts, Greek,
   partials, etc.) both in the notebook and in the exported PDF.
----------------------------------------------------------------- *)

boldR := StyleBox["r", FontWeight -> Bold];
boldV := StyleBox["v", FontWeight -> Bold];
boldN := StyleBox["n", FontWeight -> Bold];

(* inline math embedded inside a Text cell *)
mathI[box_] := Cell[BoxData[FormBox[box, TraditionalForm]],
  "InlineFormula", FontFamily -> "Times"];

(* centred display equation — its own cell *)
mathD[box_] := Cell[BoxData[FormBox[box, TraditionalForm]],
  "DisplayFormula", TextAlignment -> Center,
  FontFamily -> "Times"];

(* convenience: wrap Wolfram expressions in HoldForm and convert *)
tbx[expr_] := ToBoxes[HoldForm[expr], TraditionalForm];

(* ------- inline fragments ------- *)

vField      = RowBox[{boldV, "=",
                RowBox[{"-",
                  FractionBox[RowBox[{"\[Del]", "\[Rho]"}], "\[Rho]"]}]}];
rhoOfR      = RowBox[{"\[Rho]", "(", boldR, ")"}];
rhoOfRt     = RowBox[{"\[Rho]", "(", boldR, ",", " ", "t", ")"}];
Mc          = SubscriptBox["M", "c"];
Ac          = SubscriptBox["A", "c"];
rho0OfR     = RowBox[{SubscriptBox["\[Rho]", "0"], "(", boldR, ")"}];
lambdaMN    = SubscriptBox["\[Lambda]", "m n"];
rhoHatMN    = SubscriptBox[OverscriptBox["\[Rho]", "^"], "m n"];
tStar       = SuperscriptBox["t", "\[Star]"];
oNLogN      = RowBox[{"O", "(", "N", " ", Log, " ", "N", ")"}];
sigmaSym    = "\[Sigma]";
tZero       = SubscriptBox["t", "0"];
sigmaSqDxDy = RowBox[{SuperscriptBox["\[Sigma]", "2"], " ",
                      "\[CapitalDelta]", "x", " ",
                      "\[CapitalDelta]", "y", "/", "2"}];

(* ------- display equations ------- *)

heatEq = RowBox[{
  RowBox[{SubscriptBox["\[PartialD]", "t"], " ", rhoOfRt}],
  "=",
  RowBox[{SuperscriptBox["\[Del]", "2"], " ", rhoOfRt}],
  ",", "\[ThickSpace]\[ThickSpace]\[ThickSpace]",
  RowBox[{boldN, "\[CenterDot]", "\[Del]", "\[Rho]"}],
  "=", "0", ",", "\[ThickSpace]\[ThickSpace]",
  "on", "\[ThickSpace]",
  RowBox[{"\[PartialD]", "\[CapitalOmega]"}]}];

velEq = RowBox[{
  FractionBox[RowBox[{"d", " ", boldR}], "d t"], "=",
  RowBox[{boldV, "(", boldR, ",", " ", "t", ")"}], "=",
  RowBox[{"-",
    FractionBox[
      RowBox[{"\[Del]", " ", rhoOfRt}],
      rhoOfRt]}]}];

densityEq = RowBox[{rho0OfR, "=", FractionBox[Mc, Ac],
  ",", "\[ThickSpace]\[ThickSpace]",
  "for", "\[ThickSpace]",
  RowBox[{boldR, "\[Element]", "c"}]}];

dctExpansion = RowBox[{
  rhoOfRt, "=",
  RowBox[{
    UnderoverscriptBox["\[Sum]", RowBox[{"m", ",", "n"}], ""],
    rhoHatMN, " ",
    SuperscriptBox["e",
      RowBox[{"-", lambdaMN, " ", "t"}]], " ",
    RowBox[{"cos", "(",
      FractionBox[RowBox[{"\[Pi]", " ", "m", " ", "x"}],
                  SubscriptBox["L", "x"]], ")"}], " ",
    RowBox[{"cos", "(",
      FractionBox[RowBox[{"\[Pi]", " ", "n", " ", "y"}],
                  SubscriptBox["L", "y"]], ")"}]}]}];

lambdaEq = RowBox[{lambdaMN, "=",
  SuperscriptBox["\[Pi]", "2"], " ",
  RowBox[{"(",
    RowBox[{
      FractionBox[SuperscriptBox["m", "2"],
                  SubsuperscriptBox["L", "x", "2"]], "+",
      FractionBox[SuperscriptBox["n", "2"],
                  SubsuperscriptBox["L", "y", "2"]]}], ")"}]}];

cells = {

 Cell["Diffusion Cartograms in the Wolfram Language",
      "Title"],

 Cell["A self-contained pure-Wolfram implementation of the \
Gastner\[Dash]Newman density-equalising map algorithm.",
      "Subtitle"],

 Cell[TextData[{"Marco Thiel, April 2026 \[LongDash] ",
                ButtonBox["github.com/mthiel74/Contiguous-Cartograms",
                  BaseStyle -> "Hyperlink",
                  ButtonData -> {URL["https://github.com/mthiel74/\
Contiguous-Cartograms"], None}]}],
      "Text", TextAlignment -> Center, FontSlant -> Italic],

 (* ---- Abstract ---- *)

 Cell["Abstract", "Section"],

 Cell[TextData[{
   "Cartograms are maps in which the area of each region is deformed \
to be proportional to a statistical quantity \[LongDash] population, \
GDP, number of COVID cases, votes cast, and so on \[LongDash] while \
the topology (who is next to whom) is preserved. Gastner and Newman \
(2004) showed how to produce such maps by interpreting the quantity \
as a density and letting it diffuse to uniformity; every point of \
the map is advected by the resulting velocity field ",
   mathI[vField],
   ". This notebook presents a compact, pure-Wolfram-Language port \
of that algorithm. It ships as a paclet-style package ",
   StyleBox["CartogramWL", "Input"],
   " exposing both low-level primitives (",
   StyleBox["DiffusionSolver", "Input"], ", ",
   StyleBox["Cartogram", "Input"], ", ",
   StyleBox["CartogramTransform", "Input"],
   ") and a one-liner ",
   StyleBox["WorldCartogram[metric]", "Input"],
   " that pulls country polygons and statistical values from the \
built-in ",
   StyleBox["CountryData", "Input"],
   " and returns a ready-to-publish figure. It supports any ",
   StyleBox["ColorData", "Input"],
   " gradient, a ",
   StyleBox["GeoBackground -> \"Satellite\"", "Input"],
   " mode that textures each country with deformed satellite \
imagery, user-supplied ",
   StyleBox["Association", "Input"],
   " input, visualisation of the underlying space deformation via a \
distortion grid, and a ",
   StyleBox["PerformanceGoal -> \"Speed\"", "Input"],
   " default that runs the full world pipeline in under a minute."}],
      "Text"],

 imgCell["world_gdp_satellite.png",
         "Figure 1. Left: countries on the geographic lon/lat grid, \
each textured with satellite imagery. Right: the same countries \
with area proportional to GDP, the satellite imagery deformed \
along with the geometry. Ocean area is preserved."],

 (* ---- Section 1: Introduction ---- *)

 Cell["1. What is a cartogram?", "Section"],

 Cell[TextData[{"A geographic map of the world devotes enormous \
screen area to Russia, Canada, and the Sahara, and a few pixels \
each to Bangladesh, South Korea, or the Netherlands. If the \
question being asked is \"where do people live?\" or \"where is \
the economy?\", this is the opposite of what we want. \
Density-equalising cartograms solve this by reshaping each \
region so that its area tells the story of a variable of \
interest."}],
      "Text"],

 Cell[TextData[{"Gastner & Newman's elegant insight (PNAS 2004) is \
that a map in which density ",
   mathI[rhoOfR],
   " is uniform is exactly what we want, and diffusion is the \
natural process that takes any non-uniform density to a uniform \
one. We therefore simulate the heat equation"}],
      "Text"],

 mathD[heatEq],

 Cell[TextData[{"on a rectangular domain with Neumann boundary \
conditions (no flux through the box edge, so mass is conserved), \
and flow each geographic point through the instantaneous \
velocity"}],
      "Text"],

 mathD[velEq],

 Cell[TextData[{"As ",
   mathI[RowBox[{"t", "\[Rule]", "\[Infinity]"}]],
   " the density relaxes to its spatial mean and the advected \
points settle into their cartogram positions. Regions of high \
initial density draw space inwards and therefore expand; regions \
of low density are pushed outwards."}],
      "Text"],

 (* ---- Section 2: Method ---- *)

 Cell["2. The algorithm, step by step", "Section"],

 Cell["2.1 Rasterise polygons to a density grid", "Subsection"],

 Cell[TextData[{"Given a collection of country polygons and a \
per-country metric ",
   mathI[Mc],
   " (population, GDP, \[Ellipsis]), we set the density inside \
polygon ",
   mathI["c"],
   " to"}],
      "Text"],

 mathD[densityEq],

 Cell[TextData[{"where ",
   mathI[Ac],
   " is the geographic area of ",
   mathI["c"],
   ". Integrated over ",
   mathI["c"],
   " we recover the metric ",
   mathI[Mc],
   ". Ocean cells are filled with the mean land density so the \
spatial mean is uniform across the domain and oceans retain their \
size; only the relative sizes of countries change."}],
      "Text"],

 Cell[TextData[{"The rasterisation uses a ",
   StyleBox["RegionMember", "Input"],
   "-based subpixel test for anti-aliasing on thin/small \
countries."}],
      "Text"],

 Cell["2.2 Solve the diffusion equation", "Subsection"],

 Cell[TextData[{"The Laplacian on a rectangle with Neumann boundary \
conditions is diagonalised by the cosine basis. Wolfram's built-in ",
   StyleBox["FourierDCT[\[CenterEllipsis], 2]", "Input"],
   " and ",
   StyleBox["FourierDCT[\[CenterEllipsis], 3]", "Input"],
   " implement the DCT-II and its inverse in any number of \
dimensions with an orthonormal convention, so the full 2D \
spectrum is one call each. Expanding"}],
      "Text"],

 mathD[dctExpansion],

 Cell[TextData[{"with eigenvalues ",
   mathI[lambdaEq],
   ", the density at any time ",
   mathI["t"],
   " is evaluated in ",
   mathI[oNLogN],
   " from the cached cosine coefficients ",
   mathI[rhoHatMN],
   "."}],
      "Text"],

 Cell["2.3 Advect points through the velocity field", "Subsection"],

 Cell[TextData[{"The velocity on the grid is computed by central \
differences on the current diffused grid and interpolated \
bilinearly at any query point. Each geographic point is integrated \
from ",
   mathI[RowBox[{"t", "=", "0"}]], " to ",
   mathI[tStar],
   " using an explicit Runge\[Dash]Kutta scheme (",
   StyleBox["NDSolveValue[\[CenterEllipsis], Method -> \
\"ExplicitRungeKutta\"]", "Input"],
   ") with a 5th-order stepper, where ",
   mathI[tStar],
   " is chosen from the decay constant of the slowest non-constant \
Fourier mode."}],
      "Text"],

 Cell["2.4 Pre-smoothing", "Subsection"],

 Cell[TextData[{"Because polygon-rasterisation produces step \
functions, the initial velocity field is very stiff. We pre-smooth \
the density with a small Gaussian filter (",
   mathI[RowBox[{sigmaSym, "\[TildeTilde]", "1"}]],
   " cell), which is mathematically equivalent to starting the \
diffusion at ",
   mathI[RowBox[{tZero, "=", sigmaSqDxDy}]],
   ". A small background floor (",
   StyleBox["MeanFloor", "Input"],
   ") is added to avoid division by zero inside zero-density \
cells."}],
      "Text"],

 Cell["2.5 Snapshot-cached advection (the default Speed mode)",
      "Subsection"],

 Cell[TextData[{"A naive implementation recomputes the full \
inverse DCT of ",
   mathI[rhoOfRt],
   " at every RK sub-step of the advection ODE. With ~100+ \
sub-steps on a 192\[Times]384 grid this dominates wall-clock at \
~320 s for a single world map. The package avoids it: at \
construction time we pre-compute a small (default 60) set of \
geometrically-spaced snapshots of ",
   mathI[RowBox[{"(", "\[Rho]", ",", " ",
                  RowBox[{SubscriptBox["\[PartialD]", "x"],
                          "\[Rho]"}], ",", " ",
                  RowBox[{SubscriptBox["\[PartialD]", "y"],
                          "\[Rho]"}], ")"}]],
   ", spanning ",
   mathI[RowBox[{"t", "=", "0"}]],
   " through ",
   mathI[tStar],
   " with more samples packed near ",
   mathI[RowBox[{"t", "=", "0"}]],
   " where the fast Fourier modes collapse. During integration the \
RHS linearly interpolates between the two bracketing snapshots \
instead of doing a fresh inverse DCT. Median per-point agreement \
with the Quality solver is sub-pixel at rendering resolution and \
the full world pipeline drops from ~5\[ThinSpace]min to ~45 s \
(~7\[Times] faster). Pass ",
   StyleBox["PerformanceGoal -> \"Quality\"", "Input"],
   " to any function that integrates points (",
   StyleBox["AdvectPoints", "Input"], ", ",
   StyleBox["CartogramRun", "Input"], ", ",
   StyleBox["WorldCartogram", "Input"],
   ") for the slower, bit-exact variant."}],
      "Text"],

 (* ---- Section 3: The package ---- *)

 Cell["3. The CartogramWL package", "Section"],

 Cell["3.1 Loading", "Subsection"],

 Cell[TextData[{"Clone the repository from GitHub and load the \
package directly from the ", StyleBox["CartogramWL/", "Input"],
   " directory. The package has zero external dependencies."}],
      "Text"],

 Cell[BoxData["Get[\"/path/to/wolfram/CartogramWL/CartogramWL.wl\"]"],
      "Input"],

 Cell[TextData[{"All public symbols are then available in the ",
   StyleBox["CartogramWL`", "Input"], " context:"}],
      "Text"],

 Cell[BoxData["Names[\"CartogramWL`*\"]"], "Input"],

 Cell["3.2 Low-level primitives", "Subsection"],

 Cell[TextData[{"The algorithm is exposed as a handful of \
composable functions that all operate on a density array ",
   StyleBox["rho[[i, j]]", "Input"],
   " with row index ",
   mathI[RowBox[{"i", " ", "\[LeftRightArrow]", " ", "y"}]],
   ", column index ",
   mathI[RowBox[{"j", " ", "\[LeftRightArrow]", " ", "x"}]],
   " (NumPy convention) and a bounding box ",
   StyleBox["{xmin, ymin, xmax, ymax}", "Input"], ":"}],
      "Text"],

 Cell[BoxData[
"solver = DiffusionSolver[rho, {xmin, ymin, xmax, ymax}];\n\
rhoAtT = DensityAt[solver, t];\n\
newPts = AdvectPoints[solver, {{x1, y1}, {x2, y2}, ...}];"],
      "Input"],

 Cell[TextData[{"One level up, ",
   StyleBox["Cartogram", "Input"],
   " bundles the pre-smoothing, mean floor, and \
ocean-preservation settings; ",
   StyleBox["CartogramRun", "Input"],
   " precomputes the deformed cell-centre grid (Speed mode by \
default, ~7\[Times] faster than Quality); ",
   StyleBox["CartogramTransform", "Input"],
   " maps arbitrary points through the deformation; ",
   StyleBox["CartogramTransformPolygon", "Input"],
   " additionally densifies long edges before the transform."}],
      "Text"],

 Cell[BoxData[
"cart = Cartogram[rho, bbox, \n\
   \"MeanFloor\" -> 0.02, \"BlurSigma\" -> 1.5, \n\
   \"SeaDensity\" -> \"auto\"];\n\
cart = CartogramRun[cart, PerformanceGoal -> \"Speed\"];\n\
newPoly = CartogramTransformPolygon[cart, polyRing];"], "Input"],

 Cell["3.3 The one-liner: WorldCartogram", "Subsection"],

 Cell[TextData[{"Most end users want \"given a metric, draw the \
world\". ",
   StyleBox["WorldCartogram", "Input"],
   " wraps the whole pipeline, reads polygons and values from ",
   StyleBox["CountryData", "Input"],
   ", and returns a two-panel ",
   StyleBox["GraphicsRow", "Input"],
   " (geographic + cartogram)."}], "Text"],

 Cell[BoxData[
"WorldCartogram[\"Population\"]          (* or \"GDP\", \"GDPPerCapita\", ... *)\n\
WorldCartogram[\"GDP\", \"ColorFunction\" -> ColorData[\"TemperatureMap\"]]\n\
WorldCartogram[\"GDP\", \"Background\" -> \"Satellite\"]\n\
WorldCartogram[{\"my metric\", <|\"Germany\" -> 83., ...|>}]"],
      "Input"],

 Cell[TextData[{"A complete option list is available via ",
   StyleBox["?WorldCartogram", "Input"],
   ". The most useful options are ",
   StyleBox["\"ColorFunction\"", "Input"],
   " (any ",
   StyleBox["ColorData[\[Ellipsis]]", "Input"],
   " gradient or a raw function on ",
   mathI[RowBox[{"[", "0", ",", "1", "]"}]],
   "), ",
   StyleBox["\"Background\"", "Input"],
   " (",
   StyleBox["\"Flat\"", "Input"],
   " or ",
   StyleBox["\"Satellite\"", "Input"],
   "), ",
   StyleBox["\"MissingCountries\"", "Input"],
   " (",
   StyleBox["\"Hide\"", "Input"],
   " or ",
   StyleBox["\"ShowGrey\"", "Input"],
   "), ",
   StyleBox["\"BoundingBox\"", "Input"],
   " (useful for continent-scale maps), and ",
   StyleBox["PerformanceGoal", "Input"],
   " (",
   StyleBox["\"Speed\"", "Input"],
   " default, or ",
   StyleBox["\"Quality\"", "Input"],
   " for bit-exact output)."}], "Text"],

 (* ---- Section 4: Demos ---- *)

 Cell["4. Built-in metrics", "Section"],

 Cell[TextData[{"The four strings ",
   StyleBox["\"Population\"", "Input"], ", ",
   StyleBox["\"GDP\"", "Input"], ", ",
   StyleBox["\"GDPPerCapita\"", "Input"], " and ",
   StyleBox["\"PopulationDensity\"", "Input"],
   " are recognised directly and come with sensible defaults for \
the floor/ceiling that tame heavy-tailed outliers (Monaco on \
population density, Luxembourg on GDP per capita). The sections \
below show each of them in turn; every image on this page is \
produced by the literal one-line call shown next to it."}],
      "Text"],

 Cell["4.1 Population", "Subsection"],

 Cell[TextData[{"The canonical cartogram: each country's area is \
made proportional to how many people live there. India and China \
together house ~36 % of humanity, so they swell enormously. \
Nigeria, Indonesia, Pakistan and Bangladesh together host nearly a \
billion more people across land areas that look small on a \
Mercator projection \[LongDash] all four bulge substantially. \
Russia (11 % of global land but 2 % of humanity), Canada, \
Australia, Mongolia and the Nordic countries collapse."}],
      "Text"],

 Cell[BoxData["WorldCartogram[\"Population\"]"], "Input"],

 imgCell["world.png"],

 Cell["4.2 GDP", "Subsection"],

 Cell[TextData[{"Switching to gross domestic product (total value \
of goods and services produced in a year) rewards economic output \
rather than headcount. The United States, China, Japan, and the \
major western European economies dominate; sub-Saharan Africa and \
much of South Asia \[LongDash] densely populated but far poorer \
per head \[LongDash] shrink dramatically. This is the mirror image \
of the population cartogram wherever the two quantities diverge \
most: Japan (5\[ThinSpace]th largest economy, 11\[ThinSpace]th \
most populous) expands here but stays modestly sized in 4.1, while \
the Democratic Republic of the Congo does the opposite."}],
      "Text"],

 Cell[BoxData["WorldCartogram[\"GDP\"]"], "Input"],

 imgCell["world_gdp.png"],

 Cell["4.3 GDP per capita", "Subsection"],

 Cell[TextData[{"GDP per capita (",
   mathI[Mc],
   " = country GDP divided by country population) rewards wealth \
density rather than raw wealth. Luxembourg, Switzerland, \
Singapore, the Gulf petro-states, and the Scandinavian countries \
visibly bulge; sub-Saharan Africa, Central America, and parts of \
South Asia contract. Comparing this to the plain GDP cartogram in \
4.2: China and India shrink here \[LongDash] they are large \
economies but not rich per head \[LongDash] while small \
high-income jurisdictions grow."}], "Text"],

 Cell[BoxData["WorldCartogram[\"GDPPerCapita\"]"], "Input"],

 imgCell["world_gdp_per_capita.png"],

 Cell["4.4 Population density", "Subsection"],

 Cell[TextData[{"Population density (people per ",
   mathI[SuperscriptBox["km", "2"]],
   ") asks a third question: not \"where do people live?\" but \
\"where do people crowd?\". Bangladesh, the Netherlands, Belgium, \
South Korea, and island states such as Singapore, Malta and \
Mauritius swell out of all proportion to their land area; the \
vast, sparsely populated interiors of Russia, Canada, Australia, \
Kazakhstan, and much of sub-Saharan Africa collapse. Unlike the \
absolute-population cartogram in 4.1, large countries like China \
and the United States end up modestly sized because their \
populations are diluted across a great deal of territory."}],
      "Text"],

 Cell[BoxData["WorldCartogram[\"PopulationDensity\"]"], "Input"],

 imgCell["world_population_density.png"],

 (* ---- Section 5: Colour & satellite ---- *)

 Cell["5. Colour schemes & satellite background", "Section"],

 Cell[TextData[{"Pass any ",
   StyleBox["ColorData", "Input"],
   " gradient as ",
   StyleBox["\"ColorFunction\"", "Input"],
   ":"}], "Text"],

 Cell[BoxData[
"WorldCartogram[\"PopulationDensity\",\n\
   \"ColorFunction\" -> ColorData[\"TemperatureMap\"]]"],
      "Input"],

 imgCell["world_population_density_temperaturemap.png"],

 Cell[TextData[{"Set ",
   StyleBox["\"Background\" -> \"Satellite\"", "Input"],
   " and each country is textured with the corresponding patch \
of Wolfram's ",
   StyleBox["GeoBackground -> \"Satellite\"", "Input"],
   " imagery. On the cartogram side the imagery is deformed \
along with the country geometry, giving an effect similar to a \
shaded-relief globe being inflated or deflated locally."}],
      "Text"],

 Cell[BoxData[
"WorldCartogram[\"GDP\", \"Background\" -> \"Satellite\"]"], "Input"],

 imgCell["world_gdp_satellite.png"],

 (* ---- Section 6: Distortion grid ---- *)

 Cell["6. Visualising the distortion itself", "Section"],

 Cell[TextData[{"A clean way to see the geometry of the \
transformation is to overlay a regular 15\[Degree] lat/lon grid \
on both panels. On the left the grid is undistorted; on the \
right the 325 grid intersections are pushed into the cartogram \
by the same velocity field that reshaped the countries. The \
grid therefore encodes the local stretching/compression factor \
at every point of the domain."}], "Text"],

 imgCell["distortion_grid.png",
  "Figure 2. Left: geographic frame, regular 15\[Degree] grid in \
red, countries grey. Right: the GDP cartogram, the same grid now \
warped by the same diffusion flow. Grid cells over high-wealth \
regions shrink (fewer dollars per square millimetre); cells over \
oceans and sparsely populated continents stretch."],

 Cell[TextData[{StyleBox["CartogramDistortionGrid[cart]", "Input"],
   " returns the ",
   StyleBox["Line", "Input"],
   " primitives for the warped grid; wrap in ",
   StyleBox["Graphics", "Input"],
   " to draw, overlay on a ",
   StyleBox["WorldCartogram", "Input"],
   " result to annotate a cartogram with its own distortion \
field."}], "Text"],

 (* ---- Section 7: Custom data ---- *)

 Cell["7. Bringing your own data", "Section"],

 Cell[TextData[{"Any ",
   StyleBox["Association", "Input"],
   " of country-identifier \[RightArrow] value works in place \
of a metric string. Keys may be ",
   StyleBox["CountryData", "Input"],
   " identifiers (e.g. ",
   StyleBox["\"UnitedStates\"", "Input"],
   ") or display names (",
   StyleBox["\"United States\"", "Input"],
   ")."}],
      "Text"],

 Cell[BoxData[
"g20 = <|\"UnitedStates\" -> 332., \"China\" -> 1425., \n\
        \"India\" -> 1428., \"Japan\" -> 125., \"Germany\" -> 83., \n\
        \"UnitedKingdom\" -> 67., \"France\" -> 67.5, ...|>;\n\
\n\
WorldCartogram[{\"G20 population (millions)\", g20},\n\
   \"MissingCountries\" -> \"ShowGrey\"]"], "Input"],

 imgCell["world_g20_showgrey.png",
  "Figure 3. A cartogram driven by a hand-curated Association of \
19 G20-ish economies; every other country appears in grey and \
does not contribute to the density field, but still visibly \
deforms because diffusion is global."],

 (* ---- Section 8: Notes & references ---- *)

 Cell["8. Implementation notes", "Section"],

 Cell[TextData[{
   "\[Bullet] ",
   "Only the diffusion is done in Fourier space; advection uses \
NDSolveValue's RK45 stepper over the reconstructed velocity. Modes \
below a relative tolerance (",
   mathI[SuperscriptBox["10", RowBox[{"-", "3"}]]],
   " by default) are considered extinct, giving ",
   mathI[tStar],
   ".\n\[Bullet] ",
   "Speed mode replaces the per-sub-step inverse DCT with a \
linear interpolation between 60 log-spaced density/gradient \
snapshots, giving sub-pixel agreement with Quality at ~7\[Times] \
the throughput. Switch to ",
   StyleBox["PerformanceGoal -> \"Quality\"", "Input"],
   " when you need bit-exact reproduction.\n\[Bullet] ",
   "Heavy-tailed metrics (Monaco on population density, Luxembourg \
on GDP per capita) produce enormous per-cell densities that can \
destabilise the advection. ",
   StyleBox["WorldCartogram", "Input"],
   " clamps values at a configurable fraction of the mean and caps \
the rasterised grid at a multiple of the positive-cell mean.\n\
\[Bullet] ",
   "Long country edges are densified before being transformed, \
matching ",
   StyleBox["shapely.segmentize", "Input"],
   " on the Python side of this project.\n\[Bullet] ",
   "The algorithm is grid-convention-compatible with the Python \
implementation so the two are drop-in interchangeable for \
debugging / cross-validation."}], "Text"],

 Cell["9. Related work in the Wolfram ecosystem", "Section"],

 Cell[TextData[{"This is not the first Wolfram-Language treatment of \
cartograms. Two earlier resources are worth reading alongside this \
one:"}], "Text"],

 Cell[TextData[{
   "\[Bullet] The 2019 Wolfram Community post ",
   ButtonBox["Creating a contiguous cartogram for the US",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://community.wolfram.com/groups/-/m/t/\
1729609"], None}],
   " by George Varnavides presents a ",
   StyleBox["Manipulate", "Input"],
   "-driven US-state cartogram with an interactive colour legend, \
originally written for the 2019 Wolfram Summer School. It focuses \
on a single country, uses an iterative scaling update rather than \
diffusion, and is a nice demonstration of how much can be done \
inside a notebook with no package at all."}],
      "Text"],

 Cell[TextData[{
   "\[Bullet] ",
   StyleBox["ResourceFunction[\"Cartogram\"]", "Input"],
   " \[LongDash] ",
   ButtonBox["resources.wolframcloud.com/FunctionRepository/\
resources/Cartogram",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://resources.wolframcloud.com/\
FunctionRepository/resources/Cartogram"], None}],
   " \[LongDash] is a function-repository entry that produces a \
contiguous cartogram directly from an ",
   StyleBox["EntityClass", "Input"],
   " plus a property. It is a single-call convenience similar in \
spirit to ",
   StyleBox["WorldCartogram", "Input"],
   ", runs inside the Wolfram Cloud with no installation, and is \
the right tool when you need a quick cartogram without thinking \
about parameters."}],
      "Text"],

 Cell[TextData[{"Relative to those two, the contribution of the \
present package is (a) a pure-Wolfram port of the full \
Gastner\[Dash]Newman diffusion algorithm exposed as a handful of \
composable primitives (",
   StyleBox["DiffusionSolver", "Input"], ", ",
   StyleBox["AdvectPoints", "Input"], ", ",
   StyleBox["RasterizePolygons", "Input"],
   ", \[Ellipsis]) that can be used on any rectangular domain, not \
just a country set; (b) built-in support for ocean-area \
preservation, user-supplied ",
   StyleBox["Association", "Input"],
   " input, and satellite-textured rendering; and (c) a default ",
   StyleBox["PerformanceGoal -> \"Speed\"", "Input"],
   " mode that cuts a full world cartogram from ~5 minutes to under \
a minute while staying sub-pixel-accurate."}], "Text"],

 Cell["10. References", "Section"],

 Cell[TextData[{
   "\[Bullet] M. T. Gastner & M. E. J. Newman, \"Diffusion-based \
method for producing density-equalizing maps\", ",
   StyleBox["PNAS", FontSlant -> Italic],
   " ", StyleBox["101", FontWeight -> Bold],
   "(20), 7499\[Dash]7504 (2004). ",
   ButtonBox["doi:10.1073/pnas.0400280101",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://doi.org/10.1073/pnas.0400280101"],
                    None}], ".\n\[Bullet] ",
   "Source code and the Python twin:\n        ",
   ButtonBox["github.com/mthiel74/Contiguous-Cartograms",
     BaseStyle -> "Hyperlink",
     ButtonData -> {URL["https://github.com/mthiel74/\
Contiguous-Cartograms"], None}], ".\n\[Bullet] ",
   "Relies on Wolfram's ",
   StyleBox["CountryData", "Input"],
   ", ",
   StyleBox["GeoGraphics", "Input"],
   ", ",
   StyleBox["FourierDCT", "Input"],
   ", and ",
   StyleBox["NDSolveValue", "Input"],
   "."}], "Text"]
};

nb = Notebook[cells,
  WindowSize -> {1100, 800},
  WindowMargins -> {{Automatic, 0}, {Automatic, 0}},
  FrontEndVersion -> "15.0 for Mac OS X ARM (64-bit) (January 2026)",
  StyleDefinitions -> "Default.nb",
  Saveable -> True];

outNb  = FileNameJoin[{hereDir, "cartograms.nb"}];
outPdf = FileNameJoin[{hereDir, "cartograms.pdf"}];

Print["writing notebook -> ", outNb];
Export[outNb, nb];

Print["writing pdf      -> ", outPdf];
Export[outPdf, nb];

Print["done; nb=", FileByteCount[outNb],
      "  pdf=", FileByteCount[outPdf]];
