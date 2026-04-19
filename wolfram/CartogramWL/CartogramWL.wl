(* ::Package:: *)

(* :Title: CartogramWL *)
(* :Summary:
   Wolfram-Language implementation of the Gastner-Newman diffusion
   cartogram algorithm. Mirrors the Python package in the parent
   directory; same math, same grid convention.
*)

BeginPackage["CartogramWL`"];

DiffusionSolver::usage = "DiffusionSolver[rho0, bbox] returns an \
association that evaluates to the diffused density at time t.";

DensityAt::usage = "DensityAt[solver, t] returns rho(.,t) on the \
solver's grid.";

DensityAndGradientAt::usage = "DensityAndGradientAt[solver, t] \
returns {rho, gx, gy} at time t.";

ConvergenceTime::usage = "ConvergenceTime[solver, tol] returns a \
heuristic upper bound t* such that the slowest non-constant Fourier \
mode has decayed below tol at t*.";

AdvectPoints::usage = "AdvectPoints[solver, pts] flows a list of \
{x,y} points through v = -grad rho / rho.";

Cartogram::usage = "Cartogram[rho, bbox, opts] builds a cartogram; \
call CartogramRun on it and then CartogramTransform to map points \
through the deformation.";

CartogramRun::usage = "CartogramRun[cart] precomputes the deformed \
grid.";

CartogramTransform::usage = "CartogramTransform[cart, pts] returns \
the cartogram coordinates of pts (an Nx2 list).";

CartogramTransformPolygon::usage = "CartogramTransformPolygon[cart, \
poly] transforms a shapely-like polygon (list of {x,y} pairs) \
through the cartogram deformation, densifying long edges first.";

RasterizePolygons::usage = "RasterizePolygons[polys, values, bbox, \
{ny, nx}, subpixel] rasters an association-less list of polygons + \
per-polygon values to a density grid.";

Begin["`Private`"];

(* =====================================================================
   Grid convention: density arrays are rho[[i, j]] with i -> y, j -> x.
   bbox = {xmin, ymin, xmax, ymax}.
===================================================================== *)

(* ---- 2-D orthonormal DCT-II / III helpers. FourierDCT in Wolfram
   is already N-D by default and uses an orthonormal convention, so
   no manual axis looping and no scaling fix-up is needed. *)
dct2[m_] := FourierDCT[m, 2];
idct2[m_] := FourierDCT[m, 3];

(* ---- Eigenvalues lambda_{mn} of -Delta on the Neumann box. *)
neumannEigenvalues[ny_, nx_, Lx_, Ly_] := With[
  {mx = Range[0, nx - 1], my = Range[0, ny - 1]},
  Outer[Plus, (Pi my/Ly)^2, (Pi mx/Lx)^2]
];

(* ---- Solver builder. Returns an Association carrying enough state
   for O(N log N) evaluation of rho(t). *)
DiffusionSolver[rho0_?(MatrixQ[#, NumericQ] &), bbox_List] := Module[
  {xmin, ymin, xmax, ymax, ny, nx, Lx, Ly, dx, dy, coeffs, lam},
  {xmin, ymin, xmax, ymax} = bbox;
  If[! (xmax > xmin && ymax > ymin),
    Message[DiffusionSolver::badbbox, bbox]; Return[$Failed]];
  If[Min[rho0] <= 0,
    Message[DiffusionSolver::nonpos]; Return[$Failed]];
  {ny, nx} = Dimensions[rho0];
  Lx = xmax - xmin; Ly = ymax - ymin;
  dx = Lx/nx; dy = Ly/ny;
  coeffs = dct2[rho0];
  lam = neumannEigenvalues[ny, nx, Lx, Ly];
  <|
    "rho0" -> rho0, "bbox" -> N@bbox,
    "ny" -> ny, "nx" -> nx,
    "Lx" -> N@Lx, "Ly" -> N@Ly,
    "dx" -> N@dx, "dy" -> N@dy,
    "coeffs" -> coeffs, "lam" -> lam,
    "mean" -> Mean[Flatten[rho0]]
  |>
];
DiffusionSolver::badbbox = "Invalid bbox `1`.";
DiffusionSolver::nonpos  = "Input density must be strictly positive; \
add a background floor first.";

DensityAt[solver_Association, t_] := idct2[
  solver["coeffs"] * Exp[-solver["lam"] N[t]]
];

DensityAndGradientAt[solver_Association, t_] := Module[
  {rho, gy, gx, dx, dy},
  rho = DensityAt[solver, t];
  dx = solver["dx"]; dy = solver["dy"];
  {gy, gx} = centralGradient[rho, dy, dx];
  {rho, gx, gy}
];

(* Second-order central differences with one-sided stencils at the
   boundary. Returns {d/dy, d/dx} to match NumPy's np.gradient. *)
centralGradient[m_, dy_, dx_] := Module[{ny, nx, gy, gx},
  {ny, nx} = Dimensions[m];
  gy = ConstantArray[0., {ny, nx}];
  gx = ConstantArray[0., {ny, nx}];
  Do[gy[[i, All]] = (m[[i + 1, All]] - m[[i - 1, All]])/(2 dy),
     {i, 2, ny - 1}];
  gy[[1, All]]  = (m[[2, All]] - m[[1, All]])/dy;
  gy[[ny, All]] = (m[[ny, All]] - m[[ny - 1, All]])/dy;
  Do[gx[[All, j]] = (m[[All, j + 1]] - m[[All, j - 1]])/(2 dx),
     {j, 2, nx - 1}];
  gx[[All, 1]]  = (m[[All, 2]] - m[[All, 1]])/dx;
  gx[[All, nx]] = (m[[All, nx]] - m[[All, nx - 1]])/dx;
  {gy, gx}
];

ConvergenceTime[solver_Association, tol_:10.^-3] := With[
  {lamMin = Min[{(Pi/solver["Lx"])^2, (Pi/solver["Ly"])^2}]},
  -Log[tol]/lamMin
];

(* =====================================================================
   Bilinear sampling on the cell-centred grid with coordinates
        x_j = xmin + (j - 1/2) dx,   j in 1..nx
        y_i = ymin + (i - 1/2) dy,   i in 1..ny
   Input: field (ny x nx), xq / yq real-valued query arrays.
===================================================================== *)
bilinearSample[field_, bbox_, dx_, dy_, xq_, yq_] := Module[
  {ny, nx, xmin, ymin, fx, fy, j0, i0, j1, i1, tx, ty,
   f00, f01, f10, f11},
  {ny, nx} = Dimensions[field];
  {xmin, ymin} = bbox[[{1, 2}]];
  fx = (xq - xmin)/dx - 0.5;
  fy = (yq - ymin)/dy - 0.5;
  fx = Clip[fx, {0., nx - 1.0000000001}];
  fy = Clip[fy, {0., ny - 1.0000000001}];
  j0 = Floor[fx]; i0 = Floor[fy];
  j1 = j0 + 1; i1 = i0 + 1;
  j1 = Clip[j1, {0, nx - 1}];
  i1 = Clip[i1, {0, ny - 1}];
  tx = fx - j0; ty = fy - i0;
  (* WL arrays are 1-indexed; add 1 when looking up. *)
  f00 = MapThread[field[[#1 + 1, #2 + 1]] &, {i0, j0}];
  f01 = MapThread[field[[#1 + 1, #2 + 1]] &, {i0, j1}];
  f10 = MapThread[field[[#1 + 1, #2 + 1]] &, {i1, j0}];
  f11 = MapThread[field[[#1 + 1, #2 + 1]] &, {i1, j1}];
  (1 - ty) ((1 - tx) f00 + tx f01) + ty ((1 - tx) f10 + tx f11)
];

(* =====================================================================
   Point advection. Integrates dy/dt = rhs(t, y) with NDSolveValue
   using an explicit RK method. y is the flat state vector of
   {x1,...,xN, y1,...,yN}.
===================================================================== *)
AdvectPoints[solver_Association, pts_?(MatrixQ[#, NumericQ] &),
             OptionsPattern[{"TMax" -> Automatic, "Tol" -> 10.^-3,
                             "AccuracyGoal" -> 5,
                             "PrecisionGoal" -> 4,
                             "MaxStepFraction" -> 1/64}]] := Module[
  {n, tMax, dx, dy, bbox, rhsFn, sol, t, Y, y0, state, final},
  n = Length[pts];
  tMax = OptionValue["TMax"];
  If[tMax === Automatic, tMax = ConvergenceTime[solver,
     OptionValue["Tol"]]];
  tMax = N@tMax;
  dx = solver["dx"]; dy = solver["dy"]; bbox = solver["bbox"];

  rhsFn[tNow_?NumericQ, stateVec_?VectorQ] := Module[
    {rho, gx, gy, xq, yq, rhoQ, gxQ, gyQ, vx, vy},
    {rho, gx, gy} = DensityAndGradientAt[solver, tNow];
    xq = stateVec[[1 ;; n]];
    yq = stateVec[[n + 1 ;; 2 n]];
    rhoQ = bilinearSample[rho, bbox, dx, dy, xq, yq];
    gxQ  = bilinearSample[gx,  bbox, dx, dy, xq, yq];
    gyQ  = bilinearSample[gy,  bbox, dx, dy, xq, yq];
    rhoQ = Map[Max[#, 1.*^-300] &, rhoQ];
    vx = -gxQ/rhoQ;
    vy = -gyQ/rhoQ;
    Join[vx, vy]
  ];

  y0 = N@Join[pts[[All, 1]], pts[[All, 2]]];

  state = NDSolveValue[
    {Y'[t] == rhsFn[t, Y[t]], Y[0] == y0},
    Y, {t, 0, tMax},
    Method -> {"ExplicitRungeKutta",
               "DifferenceOrder" -> 5,
               "StiffnessTest" -> False},
    AccuracyGoal -> OptionValue["AccuracyGoal"],
    PrecisionGoal -> OptionValue["PrecisionGoal"],
    MaxStepFraction -> OptionValue["MaxStepFraction"]
  ];

  final = state[tMax];
  Transpose[{final[[1 ;; n]], final[[n + 1 ;; 2 n]]}]
];

(* =====================================================================
   Cartogram wrapper. Applies the mean_floor + blur_sigma +
   sea_density trick, builds a solver, then on Run precomputes the
   deformed grid for fast Transform queries.
===================================================================== *)
Cartogram[rho_?(MatrixQ[#, NumericQ] &), bbox_List,
          OptionsPattern[{"MeanFloor" -> 0.005, "BlurSigma" -> 1.0,
                          "SeaDensity" -> None}]] := Module[
  {prepared, mf, bs, sd, landMask, fill},
  If[Min[rho] < 0, Message[Cartogram::neg]; Return[$Failed]];
  mf = OptionValue["MeanFloor"];
  bs = OptionValue["BlurSigma"];
  sd = OptionValue["SeaDensity"];
  If[sd =!= None,
    landMask = Map[# > 0 &, rho, {2}];
    If[! MemberQ[Flatten[landMask], True],
      Message[Cartogram::allzero]; Return[$Failed]];
    fill = If[sd === "auto" || sd === Automatic,
      Mean[Pick[Flatten[rho], Flatten[landMask]]],
      N@sd];
    prepared = MapThread[If[#2, #1, fill] &, {rho, landMask}, 2];
    prepared = prepared + mf * Mean[Flatten[prepared]]
  ,
    prepared = rho + mf * Mean[Flatten[rho]]
  ];
  If[bs > 0,
    prepared = GaussianFilter[prepared, bs]
  ];
  <|
    "bbox" -> N@bbox,
    "meanFloor" -> mf,
    "blurSigma" -> bs,
    "seaDensity" -> sd,
    "solver" -> DiffusionSolver[prepared, bbox],
    "deformedGrid" -> None
  |>
];
Cartogram::neg = "rho must be non-negative.";
Cartogram::allzero = "rho is identically zero.";

CartogramRun[cart_Association,
             OptionsPattern[{"Tol" -> 10.^-3}]] := Module[
  {solver, xs, ys, pts, tMax, moved, ny, nx, xmin, ymin, xmax, ymax,
   dx, dy},
  solver = cart["solver"];
  ny = solver["ny"]; nx = solver["nx"];
  {xmin, ymin, xmax, ymax} = solver["bbox"];
  dx = solver["dx"]; dy = solver["dy"];
  xs = Table[xmin + (j - 0.5) dx, {j, 1, nx}];
  ys = Table[ymin + (i - 0.5) dy, {i, 1, ny}];
  pts = Flatten[Table[{xs[[j]], ys[[i]]}, {i, 1, ny}, {j, 1, nx}], 1];
  tMax = ConvergenceTime[solver, OptionValue["Tol"]];
  moved = AdvectPoints[solver, pts, "TMax" -> tMax];
  moved = MapThread[{Clip[#1, {xmin, xmax}], Clip[#2, {ymin, ymax}]} &,
                    {moved[[All, 1]], moved[[All, 2]]}];
  Append[cart, "deformedGrid" -> moved]
];

CartogramTransform[cart_Association,
                   pts_?(MatrixQ[#, NumericQ] &)] := Module[
  {solver, ny, nx, dispX, dispY, dx, dy, bbox, xs, ys, moved},
  If[cart["deformedGrid"] === None,
    Message[CartogramTransform::notrun]; Return[$Failed]];
  solver = cart["solver"];
  ny = solver["ny"]; nx = solver["nx"];
  dx = solver["dx"]; dy = solver["dy"];
  bbox = solver["bbox"];
  moved = cart["deformedGrid"];
  dispX = ArrayReshape[moved[[All, 1]], {ny, nx}];
  dispY = ArrayReshape[moved[[All, 2]], {ny, nx}];
  Transpose[{
    bilinearSample[dispX, bbox, dx, dy, pts[[All, 1]], pts[[All, 2]]],
    bilinearSample[dispY, bbox, dx, dy, pts[[All, 1]], pts[[All, 2]]]
  }]
];
CartogramTransform::notrun = "Call CartogramRun first.";

(* Density polygon rasteriser. Uses a subpixel sampling pattern and
   Region[poly]'s RegionMember to count subpixel-inside fractions.
   Input polys are lists of {x,y} pairs (exterior ring only, CCW).
*)
RasterizePolygons[polys_List, values_List, bbox_List,
                  shape : {_Integer, _Integer}, subpixel_:3] := Module[
  {ny, nx, xmin, ymin, xmax, ymax, dx, dy, density, offs,
   p, v, poly, area, i, j, pMinX, pMinY, pMaxX, pMaxY,
   j0, j1, i0, i1, cx, cy, fn, hits, cellDensity, invSamples},
  {ny, nx} = shape;
  {xmin, ymin, xmax, ymax} = bbox;
  dx = (xmax - xmin)/nx; dy = (ymax - ymin)/ny;
  density = ConstantArray[0., {ny, nx}];
  offs = (Range[subpixel] - 0.5)/subpixel;
  invSamples = 1./(subpixel^2);
  Do[
    poly = polys[[p]]; v = values[[p]];
    If[v == 0 || Length[poly] < 3, Continue[]];
    (* exterior polygon area by shoelace *)
    area = 0.5 Abs[Total[MapThread[
      #1[[1]] #2[[2]] - #2[[1]] #1[[2]] &,
      {poly, RotateLeft[poly]}]]];
    If[area <= 0, Continue[]];
    cellDensity = v/area;
    {pMinX, pMaxX} = MinMax[poly[[All, 1]]];
    {pMinY, pMaxY} = MinMax[poly[[All, 2]]];
    j0 = Max[1, Ceiling[(pMinX - xmin)/dx]];
    j1 = Min[nx, Ceiling[(pMaxX - xmin)/dx]];
    i0 = Max[1, Ceiling[(pMinY - ymin)/dy]];
    i1 = Min[ny, Ceiling[(pMaxY - ymin)/dy]];
    If[j0 > j1 || i0 > i1, Continue[]];
    fn = RegionMember[Polygon[poly]];
    Do[
      cy = ymin + (i - 1) dy + offs dy;
      cx = xmin + (j - 1) dx + offs dx;
      hits = Count[Flatten[Outer[fn[{#1, #2}] &, cx, cy]], True];
      If[hits > 0,
        density[[i, j]] += cellDensity * hits * invSamples
      ],
      {i, i0, i1}, {j, j0, j1}
    ],
    {p, 1, Length[polys]}
  ];
  density
];

(* Polygon transform with edge densification (WL analogue of
   shapely.segmentize). Input poly = list of {x,y}; output new poly. *)
CartogramTransformPolygon[cart_Association, poly_List,
                          OptionsPattern[{"MaxEdge" -> Automatic}]] := Module[
  {solver, maxEdge, densified, i, a, b, d, k},
  solver = cart["solver"];
  maxEdge = OptionValue["MaxEdge"];
  If[maxEdge === Automatic,
    maxEdge = Min[solver["dx"], solver["dy"]]];
  densified = {};
  Do[
    a = poly[[i]];
    b = poly[[Mod[i, Length[poly]] + 1]];
    d = Norm[b - a];
    k = Max[1, Ceiling[d/maxEdge]];
    AppendTo[densified, a];
    If[k > 1,
      Do[AppendTo[densified, a + (b - a) (s/k)], {s, 1, k - 1}]
    ],
    {i, 1, Length[poly]}
  ];
  CartogramTransform[cart, densified]
];

End[];
EndPackage[];
