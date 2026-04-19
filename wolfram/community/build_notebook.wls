#!/usr/bin/env wolframscript
(* Build wolfram/community/cartograms.nb — a Wolfram Notebook suitable
   for posting to the Wolfram Community. Also exports a PDF rendering
   of the notebook for layout inspection.

   The notebook walks through the Gastner-Newman diffusion cartogram
   method, references the CartogramWL package, shows all demo images
   (which must exist under wolfram/docs/images/), and includes
   runnable code cells. *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
imgDir  = FileNameJoin[{hereDir, "..", "docs", "images"}];
pkgPath = FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}];

loadImg[name_] := With[
  {p = FileNameJoin[{imgDir, name}]},
  If[FileExistsQ[p],
    Image[Import[p], ImageSize -> 500],
    Style["[missing image: " <> name <> "]", Red]]];

imgCell[name_, caption_String: ""] := Cell[
  BoxData[ToBoxes[
    If[caption === "",
      loadImg[name],
      Column[{
        loadImg[name],
        Pane[Style[caption, Italic, FontSize -> 10,
                   FontFamily -> "Helvetica"],
             500, Alignment -> Center]},
        Alignment -> Center, Spacings -> 0.6]]]],
  "Output", ShowCellLabel -> False,
  CellMargins -> {{66, 20}, {10, 14}}];

(* helper: inline math via TeX-like literal text; community post will
   render the $...$ delimiters. *)
tex[s_String] := s;

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

 Cell[TextData[{"Cartograms are maps in which the area of each region \
is deformed to be proportional to a statistical quantity \
\[LongDash] population, GDP, number of COVID cases, votes cast, and \
so on \[LongDash] while the topology (who is next to whom) is \
preserved. Gastner and Newman (2004) showed how to produce such \
maps by interpreting the quantity as a density and letting it \
diffuse to uniformity; every point of the map is advected by the \
resulting velocity field ",
      tex["$\\mathbf{v} = -\\nabla \\rho / \\rho$"],
      ". This notebook presents a compact, pure-Wolfram-Language \
port of that algorithm. It ships as a paclet-style package ",
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
      " input, and visualisation of the underlying space \
deformation via a distortion grid."}],
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

 Cell[TextData[{"Gastner & Newman's elegant insight (PNAS 2004) \
is that a map in which density ",
      tex["$\\rho(\\mathbf{r})$"],
      " is uniform is exactly what we want, and diffusion is the \
natural process that takes any non-uniform density to a uniform \
one. We therefore simulate the heat equation"}],
      "Text"],

 Cell[TextData[{tex["$$\\partial_t \\rho(\\mathbf{r}, t) = \\nabla^2 \
\\rho(\\mathbf{r}, t), \\qquad \\mathbf{n}\\cdot\\nabla \\rho = 0 \
\\text{ on } \\partial \\Omega$$"]}],
      "Text", TextAlignment -> Center],

 Cell[TextData[{"on a rectangular domain with Neumann boundary \
conditions (no flux through the box edge, so mass is conserved), \
and flow each geographic point through the instantaneous \
velocity"}],
      "Text"],

 Cell[TextData[{tex["$$\\frac{d\\mathbf{r}}{dt} = \\mathbf{v}\
(\\mathbf{r}, t) = -\\frac{\\nabla \\rho(\\mathbf{r}, t)}\
{\\rho(\\mathbf{r}, t)}.$$"]}],
      "Text", TextAlignment -> Center],

 Cell[TextData[{"As ",
      tex["$t \\to \\infty$"],
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
      tex["$M_c$"],
      " (population, GDP, \[Ellipsis]), we set the density inside \
polygon ",
      tex["$c$"],
      " to"}],
      "Text"],

 Cell[TextData[{tex["$$\\rho_0(\\mathbf{r}) = \\frac{M_c}{A_c} \
\\quad \\text{for } \\mathbf{r} \\in c,$$"]}],
      "Text", TextAlignment -> Center],

 Cell[TextData[{"where ", tex["$A_c$"],
      " is the geographic area of ", tex["$c$"],
      ". Integrated over ", tex["$c$"],
      " we recover the metric ", tex["$M_c$"],
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
conditions is diagonalised by the cosine basis. Wolfram's built-in \
",
      StyleBox["FourierDCT[\[CenterEllipsis], 2]", "Input"],
      " and ",
      StyleBox["FourierDCT[\[CenterEllipsis], 3]", "Input"],
      " implement the DCT-II and its inverse in any number of \
dimensions with an orthonormal convention, so the full 2D \
spectrum is one call each. Expanding"}],
      "Text"],

 Cell[TextData[{tex["$$\\rho(\\mathbf{r}, t) = \\sum_{m, n} \\hat\
\\rho_{mn} \\, e^{-\\lambda_{mn} t} \\, \\cos\\!\\left(\\frac{\\pi \
m x}{L_x}\\right) \\cos\\!\\left(\\frac{\\pi n y}{L_y}\\right)$$"]}],
      "Text", TextAlignment -> Center],

 Cell[TextData[{"with eigenvalues ",
      tex["$\\lambda_{mn} = \\pi^2 \\!\\left(\\frac{m^2}{L_x^2} + \
\\frac{n^2}{L_y^2}\\right)$"],
      ", the density at any time ",
      tex["$t$"],
      " is evaluated in ",
      tex["$\\mathcal{O}(N \\log N)$"],
      " from the cached cosine coefficients ",
      tex["$\\hat\\rho_{mn}$"],
      "."}],
      "Text"],

 Cell["2.3 Advect points through the velocity field", "Subsection"],

 Cell[TextData[{"The velocity on the grid is"}],
      "Text"],

 Cell[TextData[{tex["$$\\mathbf{v}(\\mathbf{r}, t) = -\\frac\
{\\nabla \\rho(\\mathbf{r}, t)}{\\rho(\\mathbf{r}, t)},$$"]}],
      "Text", TextAlignment -> Center],

 Cell[TextData[{"computed by central differences on the current \
diffused grid and interpolated bilinearly at any query point. \
Each geographic point is integrated from ",
      tex["$t = 0$"],
      " to ",
      tex["$t^{\\star}$"],
      " using an explicit Runge\[Dash]Kutta scheme ",
      "(",
      StyleBox[
        "NDSolveValue[\[CenterEllipsis], Method -> \
\"ExplicitRungeKutta\"]", "Input"],
      ") with a 5th-order stepper, where ",
      tex["$t^{\\star}$"],
      " is chosen from the decay constant of the slowest \
non-constant Fourier mode."}],
      "Text"],

 Cell["2.4 Pre-smoothing", "Subsection"],

 Cell[TextData[{"Because polygon-rasterisation produces step \
functions, the initial velocity field is very stiff. We \
pre-smooth the density with a small Gaussian filter (",
      tex["$\\sigma \\sim 1$"],
      " cell), which is mathematically equivalent to starting the \
diffusion at ",
      tex["$t_0 = \\sigma^2 \\Delta x \\Delta y / 2$"],
      ". A small background floor (",
      StyleBox["MeanFloor", "Input"],
      ") is added to avoid division by zero inside zero-density \
cells."}],
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

 Cell[BoxData[
   "Names[\"CartogramWL`*\"]"], "Input"],

 Cell["3.2 Low-level primitives", "Subsection"],

 Cell[TextData[{"The algorithm is exposed as a handful of \
composable functions that all operate on a density array ",
      StyleBox["rho[[i, j]]", "Input"],
      " with row index ",
      tex["$i \\leftrightarrow y$"],
      ", column index ",
      tex["$j \\leftrightarrow x$"],
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
      " precomputes the deformed cell-centre grid; ",
      StyleBox["CartogramTransform", "Input"],
      " maps arbitrary points through the deformation; ",
      StyleBox["CartogramTransformPolygon", "Input"],
      " additionally densifies long edges before the transform."}],
      "Text"],

 Cell[BoxData[
"cart = Cartogram[rho, bbox, \n\
   \"MeanFloor\" -> 0.02, \"BlurSigma\" -> 1.5, \n\
   \"SeaDensity\" -> \"auto\"];\n\
cart = CartogramRun[cart, \"Tol\" -> 2.*^-3];\n\
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
      tex["$[0,1]$"],
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
      "), and ",
      StyleBox["\"BoundingBox\"", "Input"],
      " (useful for continent-scale maps)."}], "Text"],

 (* ---- Section 4: Demos ---- *)

 Cell["4. Built-in metrics", "Section"],

 Cell["4.1 Population", "Subsection"],

 Cell[TextData[{"The canonical cartogram: India, China, and \
sub-Saharan Africa bloom outward; Russia, Canada, and Australia \
collapse."}], "Text"],

 imgCell["world.png"],

 Cell["4.2 GDP", "Subsection"],

 imgCell["world_gdp.png"],

 Cell["4.3 GDP per capita", "Subsection"],

 Cell[TextData[{"GDP per capita is a per-area quantity, so the \
cartogram rewards wealth density rather than raw wealth. \
Luxembourg, Switzerland, the Gulf states and Scandinavia \
visibly bulge; Africa and parts of South Asia contract."}],
      "Text"],

 imgCell["world_gdp_per_capita.png"],

 Cell["4.4 Population density", "Subsection"],

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
right the ",
      tex["$325$"],
      " grid intersections are pushed into the cartogram by the \
same velocity field that reshaped the countries. The grid \
therefore encodes the local stretching/compression factor at \
every point of the domain."}], "Text"],

 imgCell["distortion_grid.png",
  "Figure 2. Left: geographic frame, regular 15\[Degree] grid in \
red, countries grey. Right: the GDP cartogram, the same grid \
now warped by the same diffusion flow. Grid cells over \
high-wealth regions shrink (fewer dollars per square millimetre); \
cells over oceans and sparsely populated continents stretch."],

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
   "\[Bullet] ", "Only the diffusion is done in Fourier space; \
advection uses an RK45 stepper over the reconstructed velocity. \
Modes below a relative tolerance (",
   tex["$10^{-3}$"],
   " by default) are considered extinct, giving ",
   tex["$t^{\\star}$"],
   ".\n\[Bullet] ", "Heavy-tailed metrics (Monaco on population \
density, Luxembourg on GDP per capita) produce enormous \
per-cell densities that can destabilise the advection. \
",
   StyleBox["WorldCartogram", "Input"],
   " clamps values at a configurable fraction of the mean and \
caps the rasterised grid at a multiple of the positive-cell \
mean.\n\[Bullet] ", "Long country edges are densified before \
being transformed, matching ",
   StyleBox["shapely.segmentize", "Input"],
   " on the Python side of this project.\n\[Bullet] ", "The \
algorithm is grid-convention-compatible with the Python \
implementation so the two are drop-in interchangeable for \
debugging/cross-validation."}], "Text"],

 Cell["9. References", "Section"],

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
