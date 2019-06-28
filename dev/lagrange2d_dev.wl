(* ::Package:: *)

(*BeginPackage["findFTLE`"]*)

(* 
These lines allow automatic re-loading. If there are issues with 
saving persistent kernel configurations, etc, these lines may be 
responsible:
*)
(*Unprotect["findFTLE`*"];*)
ClearAll["findFTLE`*"];
ClearAll["findFTLE`Private`*"];

findFTLE::usage = 
"Given velocity field data, generate an estimate of 
the FTLE field over a  timescale tau"

pathPlot::usage = 
"Given velocity field data, generate a series of equal-time trajectories from uniform 
initial conditions"

streamPlot3::usage = 
"Given velocity field data, generate a series of equal-time trajectories from uniform 
initial conditions"

LinSpace::usage = 
"Given the start and end points of an interval, as well as an integer number of points, 
compute a list of evenly-spaced values along the interval."

makeMesh::usage = 
"Given the start and end points of two intervals, as well as two integer number sof points, 
computes an evenly spaced mesh along the interval."


findFTLEField::usage = 
"Given a velocity field and a set of coordinate points, compute the FTLE field 
at a given timepoint"

findKYDim::usage=
"Given a list of Lyapunov exponents, compute the Kaplan-Yorke fractal dimension";



advectPoints::usage = 
"Given a collection of initial conditions, advect them forward in time

Arguments:\[IndentingNewLine]seedPoints : N-ist of 2-Lists;\[IndentingNewLine]	A list of starting points at which to initialize trajectories\[IndentingNewLine]\[IndentingNewLine]vfield : Pair of functions in x,y;\[IndentingNewLine]	Cartesian expression of a 2D vector field\[IndentingNewLine]\[IndentingNewLine]timeLimit : Real;\[IndentingNewLine]	The amount of time to integrate

"
advectPoints[ seedPoints_,vfield_,{t_,tlo_,thi_},x_,y_]:=Module[
	{seedRules,memo1,memo2,s},
	seedRules={x0->#[[1]],y0->#[[2]]}&/@seedPoints;
	memo1[x,y,t]:=vfield[[1]];
	memo2[x,y,t]:=vfield[[2]];
	s=NDSolve[{x'[t]==memo1[x[t],y[t],t],y'[t]==memo2[x[t],y[t],t],x[0]==x0,y[0]==y0},{x[t],y[t]},{t,tlo,thi}]/.seedRules;
	s
];


advectPointsHalted::usage = 
"Given a collection of initial conditions, advect them forward in time until 
they exit the domain.

Arguments:\[IndentingNewLine]seedPoints : N-ist of 2-Lists;\[IndentingNewLine]	A list of starting points at which to initialize trajectories\[IndentingNewLine]\[IndentingNewLine]vfield : Pair of functions in x,y;\[IndentingNewLine]	Cartesian expression of a 2D vector field\[IndentingNewLine]\[IndentingNewLine]timeLimit : Real;\[IndentingNewLine]	The amount of time to integrate

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

(*Begin["`Private`"]*)

LinSpace[start_, stop_,n_]:=Array[#&,n,{start,stop}];

makeMesh[xRange_,yRange_,numX_,numY_]:=Flatten[
	Outer[
	{#1,#2}&,
	LinSpace[xRange[[1]],xRange[[-1]],numX],
	LinSpace[yRange[[1]],yRange[[-1]],numY]
	],1]




(*(*Unprotect[streamPoints];*)
Options[pathPlot] =Append[Options[ParametricPlot],streamPoints\[Rule]20];
(*Options[pathPlot] ={streamPoints\[Rule]20};*)
pathPlot[vfield_,{x_,xlo_,xhi_},{y_,ylo_,yhi_},{t_,tlo_,thi_},
	seeds0_:0,
	OptionsPattern[pathPlot]]:=Module[
	{seeds,num,numX,numY,allTraj,um},
	(*Print[Options[pathPlot]];*)
	(* Default starting locations *)
	Print[OptionValue[streamPoints]];
	num=OptionValue[streamPoints];
	{numX,numY}=Round@Sqrt[{num,num}];
	seeds=If[
	seeds0===0,
	makeMesh[{xlo,xhi},{ylo,yhi},numX,numY],
	seeds0];
	Print["there"];
	
	allTraj=advectPoints[ seeds,vfield,{t,tlo,thi},x,y];
	ParametricPlot[Evaluate[{x[t],y[t]}/.allTraj],{t,tlo,thi}]
];*)

(*Options[streamPlot3] =Options[ParametricPlot];
streamPlot3[ seedPoints_,vfield_,timeLimit_,opts1_:{nPts->5},opts:OptionsPattern[streamPlot3]]:=Module[{opts2},
(*
seedPoints : N-ist of 2-Lists;
	A list of starting points at which to initialize trajectories;

vfield : Pair of functions in x,y;
	Cartesian expression of a 2D vector field;

timeLimit : Real;
	The amount of time to integrate;
*)
Print["trii"];
Print[nPts/.opts1];
seedRules={x0->#[[1]],y0->#[[2]]}&/@seedPoints;
memo1[x,y]:=vfield[[1]];
memo2[x,y]:=vfield[[2]];
s=NDSolve[{x'[t]==memo1[x[t],y[t]],y'[t]==memo2[x[t],y[t]],x[0]==x0,y[0]==y0},{x[t],y[t]},{t,0,timeLimit}]/.seedRules;
ParametricPlot[Evaluate[{x[t],y[t]}/.s],{t,0,timeLimit},opts]
];*)

(*Unprotect[streamPoints];*)
(*Options[pathPlot] =Options[ParametricPlot],streamPoints\[Rule]20];
Options[pathPlot] ={streamPoints\[Rule]20};*)
(*Options[pathPlot] ={streamPoints\[Rule]20};*)
Options[pathPlot] =Options[ParametricPlot];
pathPlot[vfield_,{x_,xlo_,xhi_},{y_,ylo_,yhi_},{t_,tlo_,thi_},
	EvalOptions_List:{a->5},
	opts:OptionsPattern[pathPlot]]:=Module[
	{EvalOptionsDefault,EvalOptionsFull,nPts,seeds0,seeds,num,numX,numY,allTraj,um},
	
	EvalOptionsDefault={n->100, seeds1->0};
	EvalOptionsFull=Join[EvalOptions,EvalOptionsDefault];
	{nPts,seeds0}={n, seeds1}/.EvalOptionsFull;
	Print[EvalOptionsFull];
	num=nPts;
	{numX,numY}=Round@Sqrt[{num,num}];
	seeds=If[
	seeds0===0,
	makeMesh[{xlo,xhi},{ylo,yhi},numX,numY],
	seeds0];

	allTraj=advectPoints[ seeds,vfield,{t,tlo,thi},x,y];
	ParametricPlot[Evaluate[{x[t],y[t]}/.allTraj],{t,tlo,thi},opts]
];

Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];
(*Needs["DivergentColorMaps`"];*)
findFTLEField[vfield_,seeds_,{t_,tlo_,thi_},x_,y_]:=Module[
	{allTraj,startPts,stopPts,xMapPts,yMapPts,
	xMap,yMap,fMat,FtF,ftf,ftle,
	xloVal,xhiVal,yloVal,yhiVal,ftleVals,ftleField},
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

	(*fMat=Evaluate[{\!\(
\*SubscriptBox[\(\[PartialD]\), \(p\)]\({xMap[p, q], yMap[p, q]}\)\),\!\(
\*SubscriptBox[\(\[PartialD]\), \(q\)]\({xMap[p, q], yMap[p, q]}\)\)}];
	(*ftf=MatrixPower[Transpose[fMat].fMat,1];*)
	ftf=Transpose[fMat].fMat;*)
	
	fMat=Evaluate[{\!\(
\*SubscriptBox[\(\[PartialD]\), \(p\)]\({xMap[p, q], yMap[p, q]}\)\),\!\(
\*SubscriptBox[\(\[PartialD]\), \(q\)]\({xMap[p, q], yMap[p, q]}\)\)}];
	ftf=Transpose[fMat].fMat;
	
	FtF=(ftf/.{p->#[[1]],q->#[[2]]})&/@seeds;
	ftleVals=(1/Abs[thi-tlo] Log[Sqrt@(*Max@*)Eigenvalues[#]])&/@FtF;
	
	ftleField=Transpose[{Transpose[startPts][[1]],Transpose[startPts][[2]],ftleVals}];
	
	ftleField
]

findMaxFTLEField[vfield_,seeds_,{t_,tlo_,thi_},x_,y_]:=Module[
	{ftleAll},
	ftleAll=findFTLEField[vfield,seeds,{t,tlo,thi},x,y];
	{#[[1]],#[[2]],Max[#[[3]]]}&/@ftleAll
	]


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
	Print["yes"];
	
	stretchField=Transpose[{Transpose[startPts][[1]],Transpose[startPts][[2]],stretchVecs}];
	
	stretchField
]

(*End[]
Protect["findFTLE`*"];
EndPackage[]*)


(*findStopTime[functionList_]:=Module[{allMaxVals,maxVal},
allMaxVals=Transpose[Flatten[InterpolatingFunctionDomain[Head@(x[t]/.#[[1]])]&/@functionList,1]][[2]];
maxVal=Median[allMaxVals];
maxVal
]

findStopTime::usage="Given a list of trajectories produced by a simulation, find the max integration time";
cleanTrajectories[trajList_,tCutoff_]:=Module[
{jj,mm,qq,cleanedOutputs,allXRange,allYRange,bothRange,fullRangeTraj,tMax},
jj=Head/@trajList;
mm=#\[Equal]NDSolve&/@jj;
qq =If[#===True,#,False]&/@mm;
cleanedOutputs=MapThread[If[!#2,#1]&,{trajList,qq}];
cleanedOutputs = DeleteCases[cleanedOutputs,Null];
tMax=findStopTime[trajList];
(* Get rid of trajectories that terminate early*)
allYRange=InterpolatingFunctionDomain[Head[y[t]/.#[[1]]]]&/@cleanedOutputs;
allXRange=InterpolatingFunctionDomain[Head[x[t]/.#[[1]]]]&/@cleanedOutputs;
{allXRange,allYRange} = {allXRange,allYRange}\[LeftDoubleBracket]1;;All,1;;All,1,1;;All\[RightDoubleBracket];
bothRange=MapThread[{#1,#2}&,{allXRange,allYRange}];
fullRangeTraj=Select[Transpose[{cleanedOutputs,bothRange}],#1\[LeftDoubleBracket]2\[RightDoubleBracket]\[LeftDoubleBracket]1\[RightDoubleBracket]\[LeftDoubleBracket]2\[RightDoubleBracket]\[Equal]tCutoff&&#1\[LeftDoubleBracket]2\[RightDoubleBracket]\[LeftDoubleBracket]2\[RightDoubleBracket]\[LeftDoubleBracket]2\[RightDoubleBracket]\[Equal]tCutoff&];
cleanedOutputs=#\[LeftDoubleBracket]1\[RightDoubleBracket]&/@fullRangeTraj;
cleanedOutputs
]
cleanSquirmerTraj::usage="Given a list of trajectories produced by a simulation, and an integration time, remove all trajectories that did not successfully run for the full integration time.";
*)

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

findKYDim[ptList_]:=Module[
	{kyvals},
	kyvals=calculateKYDim/@(Transpose[ptList][[3]]);
	Transpose[{Transpose[ptList][[1]],Transpose[ptList][[2]],kyvals}]
]




Needs["DifferentialEquations`InterpolatingFunctionAnatomy`"];
flushingTimes[vfield_,{x_,xlo_,xhi_},{y_,ylo_,yhi_},{t_,tlo_,thi_},
	EvalOptions_List:{a->5},
	opts:OptionsPattern[pathPlot]]:=Module[
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
