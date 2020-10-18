(* ::Package:: *)

BeginPackage["ScalarWaveEquation`"];


ScalarWaveEquation::usage = " ScalarWaveEquation[Options] gives the (massive) scalar wave equation, on a spacetime with metric tensor g
(\!\(\*FractionBox[\(1\), SqrtBox[\(\(|\)\(Det[g]\)\(|\)\)]]\) \!\(\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\)(\!\(\*SqrtBox[\(\(|\)\(Det[g]\)\(|\)\)]\)\!\(\*SuperscriptBox[\(g\), \(\[Mu]\[Nu]\)]\)\!\(\*SubscriptBox[\(\[PartialD]\), \(\[Nu]\)]\))- \!\(\*SuperscriptBox[\(m\), \(2\)]\) ) \[CapitalPsi] = 0. One can provide the data of the metric through the options 'spacetime', 'diagonal_metric' or 'metric'. Please check the documentation for more information.
	";

metricDictionary::usage = " Dictionary with spacetime names and their metrics. Options of spacetimes are: 'Minkowski Spherical n=4', 'Schwarzschild n=4','Kerr n=4', 'static spherically','AdS Global n=4','Rindler-AdS n=4','Schwarzschild-AdS n=4','Static BTZ n=3','Rindler-AdS n=3.";


Options[ScalarWaveEquation]={
"spacetime"-> {},
"diagonal_metric"-> {},
"metric"-> {},
 "coordinates"-> {t,r,\[Theta],\[CurlyPhi]},
"Lorentzian"-> -1,
"ansatz"-> {}
};


(*Begin["`Private`"];*)


ScalarWaveEquation[OptionsPattern[]]:=( 

metricDictionary={
{"Minkowski Spherical n=4", DiagonalMatrix[{-1,1,r^2,r^2 Sin[\[Theta]]^2}]},  
{"Schwarzschild n=4", DiagonalMatrix[{-(1-(2 M)/r),(1-(2 M)/r)^-1,r^2,r^2 Sin[\[Theta]]^2}] },
{"Kerr n=4",{{1-(2 M r)/(r^2+a^2 Cos[\[Theta]]^2),0,0,-((2 a M r Sin[\[Theta]]^2)/(r^2+a^2 Cos[\[Theta]]^2))},{0,(r^2+a^2 Cos[\[Theta]]^2)/(a^2-2 M r+r^2),0,0},{0,0,r^2+a^2 Cos[\[Theta]]^2,0},{-((2 a M r Sin[\[Theta]]^2)/(r^2+a^2 Cos[\[Theta]]^2)),0,0,Sin[\[Theta]]^2 (a^2+r^2+(2 a^2 M r Sin[\[Theta]]^2)/(r^2+a^2 Cos[\[Theta]]^2))}}}, 
{"static spherically symmetric n=4", DiagonalMatrix[{-f[r],f[r]^-1,r^2,r^2 Sin[\[Theta]]^2}]},
{"AdS Global n=4", DiagonalMatrix[{-(1+r^2),(1+r^2)^-1,r^2,r^2 Sin[\[Theta]]^2}]} ,
{"Rindler-AdS n=4", DiagonalMatrix[{-(r^2/L^2),(1 +r^2/L^2)^-1,(1 +r^2/L^2),(1 +r^2/L^2)L^2 Sinh[\[Theta]/L]^2}]},
{"Schwarzschild-AdS n=4", DiagonalMatrix[{-(1-(2M)/r+k r^2),(1-(2M)/r+k r^2)^-1,r^2,r^2 Sin[\[Theta]]^2} ]},
{"Static BTZ n=3", DiagonalMatrix[{-(r^2-M),(r^2-M)^-1,r^2}]} ,
{"Rindler-AdS n=3", DiagonalMatrix[{-(r^2-rh^2),(r^2-rh^2)^-1,r^2}]}};


g=If[OptionValue["metric"]=!={},OptionValue["metric"],If[OptionValue["diagonal_metric"]=!={},DiagonalMatrix[ OptionValue["diagonal_metric"]], If[OptionValue["spacetime"]=!={}||Position[metricDictionary,OptionValue["spacetime"]]=!={},metricDictionary[[Position[metricDictionary,OptionValue["spacetime"]][[1,1]],2]], Print["Your waves have to travel somewhere!"];Abort[]] ]] ;


n=Length[g];
ginv=Inverse[g];
moduleDet = OptionValue["Lorentzian"]*Det[g];
coordinates=Table[OptionValue[ "coordinates"][[i]],{i,1,n}];
If[OptionValue["ansatz"]==={},ansatz=1;Sol=\[CapitalPsi][coordinates],ansatz=Product[OptionValue["ansatz"][[i]],{i,1,Length[OptionValue["ansatz"]]} ];Sol=ansatz];
If[g-DiagonalMatrix[Diagonal[g]]=!=Table[0,{i,1,n},{j,1,n}], 
waveEq= 1/ansatz*(1/moduleDet Sum[Sum[D[moduleDet*ginv[[\[Mu],\[Nu]]]D[Sol,coordinates[[\[Nu]]]], coordinates[[\[Mu]]]],{\[Nu], 1,n}],{\[Mu],1,n}] -m^2 Sol),
waveEq= FullSimplify[ 1/ansatz*(1/moduleDet Sum[Sum[D[moduleDet*ginv[[\[Mu],\[Nu]]]D[Sol,coordinates[[\[Nu]]]], coordinates[[\[Mu]]]],{\[Nu], 1,n}],{\[Mu],1,n}] -m^2 Sol)]];

Collect[waveEq,OptionValue["ansatz"]]

);


Null


(*End[];*)
EndPackage[];
