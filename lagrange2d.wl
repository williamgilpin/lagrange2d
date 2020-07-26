(* ::Package:: *)

ClearAll["findFTLE`*"];
ClearAll["findFTLE`Private`*"];

advectPoints::usage = 
"Given a collection of initial conditions, advect them forward in time. 
Returns a full set of interpolated solutions to the ODE

Arguments:\[IndentingNewLine]seedPoints : N-ist of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories\[IndentingNewLine]vfield : Pair of functions in x,y\[IndentingNewLine]	Cartesian expression of a 2D vector field
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration

Returns:
s : List of InterpolatingFunctions
	The trajectory objects solved by the integrator
"
advectPoints[ seedPoints_,vfield_,{t_,tlo_,thi_},x_,y_]:=Module[
	{seedRules,memo1,memo2,s},
	seedRules={x0->#[[1]],y0->#[[2]]}&/@seedPoints;
	memo1[x,y,t]:=vfield[[1]];
	memo2[x,y,t]:=vfield[[2]];
	s=NDSolve[{x'[t]==memo1[x[t],y[t],t],y'[t]==memo2[x[t],y[t],t],x[0]==x0,y[0]==y0},{x[t],y[t]},{t,tlo,thi}]/.seedRules;
	s
];


advectPointsFinal::usage = 
"Given a collection of initial conditions, advect them forward a fixed duration in time. 
Returns only the final locations of the points.

Arguments:\[IndentingNewLine]seedPoints : N-ist of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories\[IndentingNewLine]vfield : Pair of functions in x,y\[IndentingNewLine]	Cartesian expression of a 2D vector field
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration

Returns:
s : List of Locations
"
advectPointsFinal[ seedPoints_,vfield_,{t_,tlo_,thi_},x_,y_]:=Module[
	{seedRules,memo1,memo2,s},
	seedRules={x0->#[[1]],y0->#[[2]]}&/@seedPoints;
	memo1[x,y,t]:=vfield[[1]];
	memo2[x,y,t]:=vfield[[2]];
	s=NDSolveValue[{x'[t]==memo1[x[t],y[t],t],y'[t]==memo2[x[t],y[t],t],x[0]==x0,y[0]==y0},{x[thi],y[thi]},{t,0.999999thi,thi}]/.seedRules;
	s
];

advectPointsHalted::usage = 
"Given a collection of initial conditions, advect them forward in time until 
they exit the domain.

Arguments:\[IndentingNewLine]seedPoints : N-ist of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories\[IndentingNewLine]vfield : Pair of functions in x,y\[IndentingNewLine]	Cartesian expression of a 2D vector field

tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration

Returns:
s : List of InterpolatingFunctions
	The trajectory objects solved by the integrator
"
advectPointsHalted[seedPoints_,vfield_,
	{t_,tlo_,thi_},
	{x_,xlo_,xhi_},
	{y_,ylo_,yhi_}]:=Module[
	{seedRules,memo1,memo2,s,stopCondition},
	stopCondition=WhenEvent[(y[t]<ylo)||(y[t]>yhi)||(x[t]<xlo)||(x[t]>xhi),"StopIntegration"];
	seedRules={x0->#[[1]],y0->#[[2]]}&/@seedPoints;
	memo1[x,y,t]:=vfield[[1]];
	memo2[x,y,t]:=vfield[[2]];
	s=NDSolve[{x'[t]==memo1[x[t],y[t],t],y'[t]==memo2[x[t],y[t],t],x[0]==x0,y[0]==y0,stopCondition},{x[t],y[t]},{t,tlo,thi}]/.seedRules;
	s
];

LinSpace::usage = 
"Given the start and end points of an interval, as well as an integer number of points, 
compute a list of evenly-spaced values along the interval."
LinSpace[start_, stop_,n_]:=Array[#&,n,{start,stop}];

makeMesh::usage = 
"Given the start and end points of two intervals, as well as two integer numbers of points, 
computes an evenly spaced mesh along the interval."
makeMesh[xRange_,yRange_,numX_,numY_]:=Flatten[
	Outer[
	{#1,#2}&,
	LinSpace[xRange[[1]],xRange[[-1]],numX],
	LinSpace[yRange[[1]],yRange[[-1]],numY]
	],1]



pathPlot::usage = 
"Given velocity field data, generate a series of equal-time trajectories from uniform 
initial conditions

Arguments:
xlo: Real
	The lower bound on x
xhi: Real
	The upper bound on x
ylo: Real
	The lower bound on y
yhi: Real
	The upper bound on y
tlo: Real
	The starting value for t
thi: Real
	The ending value for t

Arguments (Optional):
n : Integer
	The number of points to use to discretize the domain 
seeds1 : List of {Real, Real}
	An explicit list of starting points to use for the integration
Options[ParametricPlot] : Any plotting options that would normally be passed to ParametricPlot

"
Options[pathPlot] =Options[ParametricPlot];
pathPlot[vfield_,{x_,xlo_,xhi_},{y_,ylo_,yhi_},{t_,tlo_,thi_},
	EvalOptions_List:{a->5},
	opts:OptionsPattern[pathPlot]]:=Module[
	{EvalOptionsDefault,EvalOptionsFull,nPts,seeds0,seeds,num,numX,numY,allTraj,um},
	
	EvalOptionsDefault={n->100, seeds1->0};
	EvalOptionsFull=Join[EvalOptions,EvalOptionsDefault];
	{nPts,seeds0}={n, seeds1}/.EvalOptionsFull;
	(*Print[EvalOptionsFull];*)
	num=nPts;
	{numX,numY}=Round@Sqrt[{num,num}];
	seeds=If[
	seeds0===0,
	makeMesh[{xlo,xhi},{ylo,yhi},numX,numY],
	seeds0];

	allTraj=advectPoints[ seeds,vfield,{t,tlo,thi},x,y];
	ParametricPlot[Evaluate[{x[t],y[t]}/.allTraj],{t,tlo,thi},opts]
];



findFTLEField::usage = 
"Given a velocity field and a set of coordinate points, compute the FTLE field 
at a given timepoint

Arguments:\[IndentingNewLine]vfield : Pair of functions in x,y\[IndentingNewLine]	Cartesian expression of a 2D vector field
seeds : N-ist of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration

Returns:
ftleField : List of 3-Lists 
	A list of {x,y,{lamba1,lambda2} points denoting the two finite time
	Lyapunov exponents at each location
"
Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];
findFTLEField[vfield_,seeds_,{t_,tlo_,thi_},x_,y_]:=Module[
	{allTraj,startPts,stopPts,xMapPts,yMapPts,
	xMap,yMap,fMat,FtF,ftf,ftle,
	xloVal,xhiVal,yloVal,yhiVal,ftleVals,ftleField},
	
	
	(*allTraj=advectPoints[ seeds,vfield,{t,tlo,thi},x,y];
	{x[t],y[t]}/.allTraj;
	startPts=Evaluate[{x[t],y[t]}/.allTraj]/.t->tlo;
	stopPts = Evaluate[{x[t],y[t]}/.allTraj]/.t->thi;
	startPts=startPts[[1;;All,1,1;;All]];
	stopPts=stopPts[[1;;All,1,1;;All]];*)
	startPts=seeds;
	stopPts=advectPointsFinal[ seeds,vfield,{t,tlo,thi},x,y];
	{x[t],y[t]}/.allTraj;
	
	xMapPts=MapThread[{#1,#2}&,{startPts,Transpose[stopPts][[1]]}];
	yMapPts=MapThread[{#1,#2}&,{startPts,Transpose[stopPts][[2]]}];
	xMap = Interpolation[xMapPts];
	yMap = Interpolation[yMapPts];

	fMat=Evaluate[{\!\(
\*SubscriptBox[\(\[PartialD]\), \(p\)]\({xMap[p, q], yMap[p, q]}\)\),\!\(
\*SubscriptBox[\(\[PartialD]\), \(q\)]\({xMap[p, q], yMap[p, q]}\)\)}];
	ftf=Transpose[fMat].fMat;
	
	FtF=(ftf/.{p->#[[1]],q->#[[2]]})&/@seeds;
	ftleVals=(1/Abs[thi-tlo] Log[Sqrt@(*Max@*)Eigenvalues[#]])&/@FtF;
	
	ftleField=Transpose[{Transpose[startPts][[1]],Transpose[startPts][[2]],ftleVals}];
	
	ftleField
]

findMaxFTLEField::usage = 
"Given a velocity field and a set of coordinate points, compute the max FTLE field 
at a given timepoint

Arguments:\[IndentingNewLine]vfield : Pair of functions in x,y\[IndentingNewLine]	Cartesian expression of a 2D vector field
seeds : N-ist of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration

Returns:
ftleField : List of 3-Lists 
	A list of {x,y,lambda} points denoting the max finite time
	Lyapunov exponents at each location
";
findMaxFTLEField[vfield_,seeds_,{t_,tlo_,thi_},x_,y_]:=Module[
	{ftleAll},
	ftleAll=findFTLEField[vfield,seeds,{t,tlo,thi},x,y];
	{#[[1]],#[[2]],Max[#[[3]]]}&/@ftleAll
	]


findStretchlines::usage = 
"Given a velocity field and a set of coordinate points, compute the vector field 
corresponding to the maximum or minimum stretching direction

Arguments:\[IndentingNewLine]vfield : Pair of functions in x,y\[IndentingNewLine]	Cartesian expression of a 2D vector field
seeds : N-ist of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration
ordering : String
	Whether to compute the field for maximum or minimum streching

Returns:
stretchField : List of 3-Lists 
	A list of {x,y,{v1, v2}} points denoting the stretching field
	at each location
";
findStretchlines[vfield_,seeds_,{t_,tlo_,thi_},x_,y_,ordering_:"Positive"]:=Module[
	{allTraj,startPts,stopPts,xMapPts,yMapPts,
	xMap,yMap,fMat,FtF,ftf,ftle,sortEigVecs,
	xloVal,xhiVal,yloVal,yhiVal,stretchVecs,stretchField},
	allTraj=advectPoints[ seeds,vfield,{t,tlo,thi},x,y];
	{x[t],y[t]}/.allTraj;
	(* add a cleaned outputs step *)
	startPts=Evaluate[{x[t],y[t]}/.allTraj]/.t->tlo;
	stopPts = Evaluate[{x[t],y[t]}/.allTraj]/.t->thi;
	startPts=startPts[[1;;All,1,1;;All]];
	stopPts=stopPts[[1;;All,1,1;;All]];
	xMapPts=MapThread[{#1,#2}&,{startPts,Transpose[stopPts][[1]]}];
	yMapPts=MapThread[{#1,#2}&,{startPts,Transpose[stopPts][[2]]}];
	xMap = Interpolation[xMapPts];
	yMap = Interpolation[yMapPts];

	fMat=Evaluate[{\!\(
\*SubscriptBox[\(\[PartialD]\), \(p\)]\({xMap[p, q], yMap[p, q]}\)\),\!\(
\*SubscriptBox[\(\[PartialD]\), \(q\)]\({xMap[p, q], yMap[p, q]}\)\)}];
	ftf=Transpose[fMat].fMat;
	
	FtF=(ftf/.{p->#[[1]],q->#[[2]]})&/@seeds;
	
	sortEigVecs=If[ordering==="Positive",
	(#[[2]][[Reverse@Ordering[#[[1]]]]])[[1]]&,
	(#[[2]][[Ordering[#[[1]]]]])[[1]]&];
	stretchVecs=sortEigVecs/@(Eigensystem/@FtF);
	stretchField=Transpose[{Transpose[startPts][[1]],Transpose[startPts][[2]],stretchVecs}];
	
	stretchField
]


calculateKYDim::usage=
"
Given a List of lyapunov exponents, calculate the associated Kaplan-Yorke
fractal dimension.
"
calculateKYDim[ptList_]:=Module[
	{cumSumList,sortedList,flipInd,kyd},
	sortedList=Sort[ptList,Greater];
	kyd =If[Length[ptList]>2,
	cumSumList=Accumulate@Sort[ptList];
	flipInd=(Flatten@Position[cumSumList,_?Positive])[[-1]];
	flipInd + Total[sortedList[[;;flipInd]]]/Abs[sortedList[[flipInd+1]]],
	1+sortedList[[1]]/Abs[sortedList[[2]]]
	];
	kyd
]

findKYDim::usage=
"
Given a list of Lyapunov exponents, compute the Kaplan-Yorke fractal dimension

Arguments:
ptList : A List of {x,y,{lamba1,lambda2}} points denoting the two Lyapunov 
	exponents at each location


Returns:
kyField : List of 3-Lists 
	A list of {x,y,K} points denoting the Kaplan-Yorke dimension at
	each location
";
findKYDim[ptList_]:=Module[
	{kyvals,kyField},
	kyvals=calculateKYDim/@(Transpose[ptList][[3]]);
	kyField=Transpose[{Transpose[ptList][[1]],Transpose[ptList][[2]],kyvals}];
	kyField
]



flushingTimes::usage = 
"Given a velocity field and a domain, calculate the flushing time field

Arguments:\[IndentingNewLine]vfield : Pair of functions in x,y,t\[IndentingNewLine]	Cartesian expression of a 2D time-dependent vector field
xlo: Real
	The lower bound on x
xhi: Real
	The upper bound on x
ylo: Real
	The lower bound on y
yhi: Real
	The upper bound on y
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration (the maximum flushing time)

Arguments (Optional):
n : Integer
	The number of points to use to discretize the domain 
seeds1 : List of {Real, Real}
	An explicit list of starting points to use for the integration
Options[ParametricPlot] : Any plotting options that would normally be
					      passed to ParametricPlot

Returns:
flushField : List of 3-Lists
	A list of {x,y,t} points denoting the flushing time at 
	each spatial location
"
Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];
flushingTimes[vfield_,{x_,xlo_,xhi_},{y_,ylo_,yhi_},{t_,tlo_,thi_},
	EvalOptions_List:{a->5}]:=Module[
	{EvalOptionsDefault,EvalOptionsFull,nPts,seeds0,seeds,num,numX,numY,allTraj,um,flushField,allMaxVals},
	
	EvalOptionsDefault={n->100, seeds1->0};
	EvalOptionsFull=Join[EvalOptions,EvalOptionsDefault];
	{nPts,seeds0}={n, seeds1}/.EvalOptionsFull;
	(*Print[EvalOptionsFull];*)
	num=nPts;
	{numX,numY}=Round@Sqrt[{num,num}];
	seeds=If[
	seeds0===0,
	makeMesh[{xlo,xhi},{ylo,yhi},numX,numY],
	seeds0];

	allTraj=advectPointsHalted[ seeds,vfield,{t,tlo,thi},{x,xlo,xhi},{y,ylo,yhi}];
	allMaxVals=Transpose[Flatten[InterpolatingFunctionDomain[Head@(x[t]/.#[[1]])]&/@allTraj,1]][[2]];
	flushField=Transpose[{Transpose[seeds][[1]],Transpose[seeds][[2]],allMaxVals}];
	flushField
];


animateFlow::usage = 
"
Given a velocity field and a domain, create a video of the flow as it 
evolves in time

Arguments:\[IndentingNewLine]vfield : Pair of functions in x,y,t\[IndentingNewLine]	Cartesian expression of a 2D time-dependent vector field
seeds : Length N List of 2-Lists\[IndentingNewLine]	A list of starting points at which to initialize trajectories
xlo: Real
	The lower bound on x
xhi: Real
	The upper bound on x
ylo: Real
	The lower bound on y
yhi: Real
	The upper bound on y
tlo : Real
	The time to start integration\[IndentingNewLine]thi : Real\[IndentingNewLine]	The time to stop integration (the maximum flushing time)
frameRate : Integer
	The number of frames per second of the video
playbackSpeed : Integer
	The playback rate, 1x, 2x, etc.
filename : String
	The name of the video. The extension .avi may also be used

Returns: None
"

animateFlow[vfield_,seeds_,{t_,tlo_,thi_},{x_,xlo_,xhi_},{y_,ylo_,yhi_},frameRate_:10,
			playbackSpeed_:1,filename_:"animation.mov"]:=Module[
	{advectedPoints,qq,m,allTraj,startPts},
	advectedPoints=advectPoints[ seeds,vfield,{t,tlo,thi},x,y];
	(* Pass Plot options *)
	m=Manipulate[
	Show[{
	ListPlot[Flatten[({x[t],y[t]}/.advectedPoints)/.t->tval,1],PlotStyle->Black]
	},PlotRange->{{xlo,xhi},{ylo,yhi}},AspectRatio->Automatic,Axes->False],
	{tval,tlo,thi},
	ControlType->None,
	AutorunSequencing->{{1,1/playbackSpeed}}
	];
	Export[filename,m,"FrameRate"->frameRate]
]


fitVField::usage = 
"
Given an experimental velocity field dataset, fit a smooth velocity field function
 using interpolation

Arguments:\[IndentingNewLine]dataset : 5xN List of velocity field values {t, x, y, vx, vy} corresponding to
		a 2D time-dependent vector field

Returns:
{vx, vy} : List of 2 functions
	Two functions corresponding to the x and y components of the velocity field,
	with call structures vx[t,x,y], vy[t,x,y]
"

fitVField[dataset_]:=Module[{formattedDataUX,formattedDataUY},
	formattedDataUX={{#[[1]],#[[2]],#[[3]]},#[[4]]}&/@Transpose[dataset[[{1,2,3,4},1,All]]];
	formattedDataUY={{#[[1]],#[[2]],#[[3]]},#[[4]]}&/@Transpose[dataset[[{1,2,3,5},1,All]]];
	{vx,vy}=Interpolation/@{formattedDataUX,formattedDataUY}
]
