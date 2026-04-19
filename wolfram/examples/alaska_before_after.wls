#!/usr/bin/env wolframscript
(* Build a side-by-side GDP cartogram using the two ring sources in
   turn. Each side shows both the geographic panel and the cartogram
   panel of the *full-world* render, then crops both panels to the
   North American quadrant so Alaska is clearly visible. That makes
   it obvious that Alaska stays geographic-sized when we rasterise
   from CountryData alone but grows alongside the mainland when we
   splice in the AdministrativeDivision polygons via
   CartogramWL`CountryPolygonRings.

   Output: wolfram/docs/images/alaska_before_after.png *)

hereDir = DirectoryName[If[$ScriptCommandLine =!= {}, $InputFileName,
                           NotebookFileName[]]];
Get[FileNameJoin[{hereDir, "..", "CartogramWL", "CartogramWL.wl"}]];

metric = "GDP";
bbox   = {-180., -60., 180., 85.};
{xmin, ymin, xmax, ymax} = bbox;
skip   = {"Antarctica", "FrenchSouthernTerritories", "Bouvet Island"};

toRings[poly_Polygon] := Module[{r},
  r = First[poly[[1]]];
  Map[{#[[2]], #[[1]]} &, r, {2}]];

ringsCountryDataOnly[c_] := Module[{p = CountryData[c, "Polygon"]},
  If[Head[p] === Polygon, toRings[p], {}]];

buildPanels[ringFn_] := Module[
  {countries, entries, values, floor, pairs, polys, rvals, rho,
   posMean, cart, trans, newEntries, minP, maxP, logMin, logMax,
   colorFn, draw},
  countries = Select[CountryData["Countries"],
    ! MemberQ[skip, CountryData[#, "Name"]] &];
  entries = Reap[Do[
    With[{val = CountryData[c, metric], rs = ringFn[c]},
      If[Head[val] === Quantity && rs =!= {},
        Sow[<|"name"  -> CountryData[c, "Name"],
              "value" -> N@QuantityMagnitude[val],
              "rings" -> rs|>]]],
    {c, countries}]][[2, 1]];
  values  = entries[[All, "value"]];
  floor   = 0.01 Mean[Select[values, # > 0 &]];
  entries = Map[Append[#, "value" -> Max[#["value"], floor]] &,
                entries];
  pairs = Flatten[Map[
    Function[e, Module[{rs, areas, tot},
      rs    = e["rings"];
      areas = Map[0.5 Abs[Total[MapThread[
                #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
                {#, RotateLeft[#]}]]] &, rs];
      tot   = Total[areas];
      If[tot > 0,
        MapThread[{#1, e["value"] * #2/tot} &, {rs, areas}],
        {}]]],
    entries], 1];
  {polys, rvals} = Transpose[pairs];
  rho     = CartogramWL`RasterizePolygons[polys, rvals, bbox,
              {192, 384}, 2];
  posMean = Mean[Select[Flatten[rho], # > 0 &]];
  rho     = Map[Min[#, 200. posMean] &, rho, {2}];
  cart    = CartogramWL`Cartogram[rho, bbox,
    "MeanFloor" -> 0.02, "BlurSigma" -> 1.5,
    "SeaDensity" -> "auto"];
  cart    = CartogramWL`CartogramRun[cart, "Tol" -> 3.*^-3,
              PerformanceGoal -> "Speed", "Snapshots" -> 60];
  trans[ring_] := CartogramWL`CartogramTransformPolygon[cart,
                    ring, "MaxEdge" -> 2.0];
  newEntries = Map[
    Append[#, "rings" -> Map[trans, #["rings"]]] &, entries];
  {entries, newEntries}];

Print["[1/2] pipeline without Alaska/Hawaii..."];
{entriesNo,  newEntriesNo}  = buildPanels[ringsCountryDataOnly];
Print["[2/2] pipeline with Alaska + Hawaii..."];
{entriesYes, newEntriesYes} = buildPanels[CartogramWL`CountryPolygonRings];

(* -- Use two distinct solid fills: the rest of the world in pale
      grey, the United States in a saturated orange, so Alaska and
      Hawaii are unmistakably visible when they're part of the US
      ring set. *)
usColor     = RGBColor[1, 0.5, 0.12];
otherColor  = GrayLevel[0.88];

drawWithHighlight[entryList_, plotRange_, title_] := Graphics[
  Table[Module[{e = entryList[[k]], isUS},
    isUS = e["name"] === "United States";
    {EdgeForm[Directive[Black, Thickness[0.0004]]],
     FaceForm[If[isUS, usColor, otherColor]],
     Polygon /@ e["rings"]}],
    {k, 1, Length[entryList]}],
  Frame -> True, FrameTicks -> None,
  PlotRange -> plotRange,
  PlotLabel -> title, AspectRatio -> Automatic,
  ImageSize -> 700];

geoRange  = {{-175., -50.}, {15., 75.}};

geoPanelNo  = drawWithHighlight[entriesNo,
   geoRange, "Geographic \[LongDash] CountryData polygon only"];
geoPanelYes = drawWithHighlight[entriesYes,
   geoRange,
   "Geographic \[LongDash] CountryPolygonRings (+Alaska, +Hawaii)"];

(* Cartogram panel: widen the crop so the inflated US fits. *)
cartRange = {{-200., 30.}, {-10., 90.}};

cartPanelNo  = drawWithHighlight[newEntriesNo,  cartRange,
   "Cartogram \[LongDash] Alaska + Hawaii absent from density field"];
cartPanelYes = drawWithHighlight[newEntriesYes, cartRange,
   "Cartogram \[LongDash] Alaska + Hawaii included"];

fig = GraphicsGrid[{
  {geoPanelNo,  geoPanelYes},
  {cartPanelNo, cartPanelYes}},
  ImageSize -> 1600, Spacings -> 20,
  PlotLabel -> Style[
    "North America \[LongDash] GDP cartogram \[LongDash] \
effect of splicing in Alaska + Hawaii\n\
(United States highlighted in orange)",
    FontSize -> 14, Bold]];

outPath = FileNameJoin[{hereDir, "..", "docs", "images",
                        "alaska_before_after.png"}];
Export[outPath, fig, ImageResolution -> 150];
Print["wrote ", outPath];
