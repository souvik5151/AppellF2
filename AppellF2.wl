(* ::Package:: *)

BeginPackage["AppellF2`"]


Print["AppellF2.wl v1.3\n","Authors : Souvik Bera & Tanay Pathak"];

(*
What is new?

1. Accelerated technique is used. The double sum is rewritten with 
the Mathematica 2F1 command, which results only one sum to be evaluated.
This makes the sum to run faster compared to the double sum.

2. AppellF2 command is now replaced by F2. This is due to the fact that:
there is inbuild AppellF2 command in Mathematica version 13.3

3. Only the AppellF2 command in the previous version is changed. 
Others are kept intact. Check others before using!!!! <<< IMPORTANT >>>


*)


F2::usage="The command gives the numerical value of the Appell Function F2.
 F2[a, b1, b2, c1, c2,x,y,precision,terms, F2show-> True]";
F2ROC::usage="The command gives the region of convergence (ROC) of the series along with the given point
 F2ROC[{x,y},series_number,{plot range}]";
F2findall::usage="The command finds all the analytic continuations that are valid for the given point
 F2findall[{x,y}]";
F2evaluate::usage="The command gives the value of the series at a given point and Pochhammer parameters
 F2evaluate[series_number,{a,b1,b2,c1,c2,x,y}, precision, terms]";
F2expose::usage="The command exposes the region of convergence (ROC) and the expression of the analytic continuation of F2[a,b1,b2,c1,c2,x,y,m,n]
 F2expose[series_number]";


Begin["`Private`"]


F2expose[list_]:=Module[{i},
If[list>24,Print["Total number of ACs is 24"];Abort[];];
ClearAll[Global`x,Global`y,Global`a1,Global`a2,Global`b1,Global`b2,Global`c,Global`m,Global`m];
Return[{F2serieswsum[Global`a1,Global`a2,Global`b1,Global`b2,Global`c,Global`x,Global`y,Global`m,Global`n][[list[[1]],1]],
F2serieswsum[Global`a1,Global`a2,Global`b1,Global`b2,Global`c,Global`x,Global`y,Global`m,Global`n][[list[[1]],2]](*/.Gamma[x_. +1]-> x Gamma[x]*)/.Re[x_]-> x (*//Simplify//Expand*)}]];



F2findall[{x_?NumberQ, y_?NumberQ}]:=Module[{test,pos,a, b1, b2, c1, c2,m},

test=F2series[a, b1, b2, c1, c2,Rationalize[x],Rationalize[y],m][[All,1]];
pos=Position[test,True];
Return[Flatten[pos]];


]


F2evaluate[list_/;(IntegerQ[list]&&list>0), para_List, p_,terms_]:=Module[{a, b1, b2, c1, c2,x,y,result,p0,seriesac,m},
p0=p+10;

If[Length[para]=!=7,Abort[];];
{a, b1, b2, c1, c2,x,y}=SetPrecision[para,p0];



If[F2series[a, b1, b2, c1, c2,Rationalize[x],Rationalize[y],m][[list,1]],{},Print["The point does not lie in ", list];Abort[];];
seriesac[m1_]:= F2series[a, b1, b2, c1, c2,x,y,m1][[list,2]];

Block[{sum,m1,n1},{sum=0;Monitor[For[m1=0,m1<=terms,m1++,
{sum=sum+SetPrecision[seriesac[m1],p0];};(*Pause[10^-6];*)],{SetPrecision[sum,p0],m1}];result=sum;}];


Return[SetPrecision[Chop[result],p]]]


F2ROC[point___List,list_/;(IntegerQ[list]&&list>0),range_List]:=Module[{i,a,b,c,d,e,x,y,m,roc},
roc[x_,y_]=F2series[a,b,c,d,e,x,y,m][[list,1]];
Show[ListPlot[{point},PlotStyle->Directive[PointSize[Medium],Red],PlotRange->{{range[[1]],range[[2]]},{range[[1]],range[[2]]}}],RegionPlot[roc[x,y],{x,range[[1]]-1,range[[2]]+1},{y,range[[1]]-1,range[[2]]+1},PlotPoints->100],PlotRange->{{range[[1]],range[[2]]},{range[[1]],range[[2]]}},AspectRatio-> 1]
]


Options[F2]={F2show->False};
F2::poch="Pochhammer parameters do not obey the constraint";
F2::val="The point is not covered";
F2[args___] := f2core[F2,args];

f2core[F2,a0_,b0_,c0_,d0_,e0_,x0_,y0_,p_,terms_,OptionsPattern[F2]]:=Module[{aa,bb1,bb2,cc1,cc2,xx,yy,a, b1, b2, c1, c2,x1,y1,
x,y,m,n,m0,m1,n0,n1,v,i,test1,test2,pos1,pos2,
result,p0,integrand2,integrand3,seriesac,seriesselect,allroc,acs,parameters,
eachserlist,convrate1,eachser,selectedseries},
p0=p+10;m1=10;
(*Print[{a0,b0,c0,d0,e0,x0,y0}];*)
(*Label[start];
Print["start"];*)

(*a=SetPrecision[a0,p0]; b1=SetPrecision[b0,p0]; b2=SetPrecision[c0,p0];
c1=SetPrecision[d0,p0]; c2=SetPrecision[e0,p0];x=SetPrecision[x0,p0]; 
y=SetPrecision[y0,p0];*)
result=0;
parameters = SetPrecision[{a0,b0,c0,d0,e0,x0,y0},p0];


(* If x=0 or y=0*)

If[x0===0.0||x0===0,result=Hypergeometric2F1[a,b2,c2,y];Goto[end];];
If[y0===0.0||y0===0,result=Hypergeometric2F1[a,b1,c1,x];Goto[end];];






Off[General::infy,General::indet];
(* check the region *)
allroc[x_,y_]=F2series[a, b1, b2, c1, c2,x,y,m][[All,1]];
acs[{a_, b1_, b2_, c1_, c2_,x_,y_},m_]= F2series[a, b1, b2, c1, c2,x,y,m][[All,2]];
test1=allroc[Rationalize[x0],Rationalize[y0]];
(*Print[test1];*)

PrintTemporary["Finding a proper AC..."];

pos1=Flatten[Position[test1,True]];

If[pos1==={},Message[F2::val];Abort[];];

Print["valid series : ",pos1];



(*test2=Table[{Max[Norm[#]&/@(Abs[{Limit[Simplify[(#1/.{m->m+ 1})/(#1)],{m-> m0},PerformanceGoal-> "Speed"]}/.{m0-> m1}/.{aa-> a,bb1-> b1,bb2-> b2,cc1-> c1, cc2-> c2,xx-> x,yy-> y}]&/@If[Head[#]===Plus,List@@#,{#}]&@Expand[F2series[aa,bb1,bb2,cc1,cc2,xx,yy,m][[pos1[[1]],2]]])],pos1[[i]][[1]]},{i,Length[pos1]}];
*)

test2=Table[eachser = acs[{a, b1, b2, c1, c2, x, y},m][[pos1[[i]]]];
eachserlist = If[Head[#]===Plus,List@@#,{#}]&@eachser;
convrate1 = Abs[N[((eachserlist/.m->m+ 1)/(eachserlist))/.{m-> m1}
/.MapThread[Rule,{{a, b1, b2, c1, c2, x, y},parameters}],10]];
{Max[convrate1],pos1[[i]]}
,{i,1,Length[pos1]}];

(*Print[test2];*)

seriesselect=SortBy[test2,First];
Print["convergence rates :",seriesselect];

For[i=1,i<= Length[test2],i++,
pos2=seriesselect[[i]][[2]];
seriesac[m11_]= F2series[a, b1, b2, c1, c2,x,y,m11][[pos2,2]];
If[MemberQ[{"0","Indeterminate"},ToString[Sum[seriesac[m11],{m11,0,1}]]],Continue[],Break[]]];


Print["selected series : ",pos2];

PrintTemporary["Evaluating sum..."];

selectedseries[{a_, b1_, b2_, c1_,c2_,x_,y_},m_]=F2series[a, b1, b2, c1,c2, x, y,m][[pos2,2]]; 
DistributeDefinitions[terms,p,parameters,p0];
result = N[ParallelSum[selectedseries[SetPrecision[parameters,p0],m],{m,0,terms}],p0];

(*Block[{sum,m11},{ result=ParallelSum[seriesac[m11],{m11,0,terms}];}];
*)
(*NO need of F2show option*)
(*If[OptionValue[F2show]===True, Block[{sum,m11,n11},{sum=0;Monitor[For[m11=0,m11<=terms,m11++,
{sum=sum+SetPrecision[F2series[a, b1, b2, c1, c2,x,y,m11][[pos2,2]],p0];}];,{SetPrecision[sum,p0],m11}];result=sum;}];,
;


];
*)


Label[end];
Return[SetPrecision[Chop[result],p]]
]



(* ::Subchapter:: *)
(*Set of 44 ACs*)


F2series[a_,b1_,b2_,c1_,c2_,x_,y_,m_]:= {{Abs[x]+Abs[y]<1&&Abs[x]<1,(x^m HypergeometricPFQ[{b2,a+m},{c2},y] Pochhammer[a,m] Pochhammer[b1,m])/(m! Pochhammer[c1,m])},{Abs[x/(1-y)]<1&&Abs[1-y]<1,((x/(1-y))^m (1-y)^(-a-b2+c2) Gamma[a+b2-c2] Gamma[c2] HypergeometricPFQ[{-b2+c2,-a+c2-m},{1-a-b2+c2-m},1-y] Pochhammer[b1,m] Pochhammer[a+b2-c2,m])/(m! Gamma[a] Gamma[b2] Pochhammer[c1,m])+(x^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{b2,a+m},{1+a+b2-c2+m},1-y] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])},{Abs[y/(1-x)]<1&&Abs[1-x]<1,((1-x)^(-a-b1+c1) (y/(1-x))^m Gamma[a+b1-c1] Gamma[c1] HypergeometricPFQ[{-b1+c1,-a+c1-m},{1-a-b1+c1-m},1-x] Pochhammer[b2,m] Pochhammer[a+b1-c1,m])/(m! Gamma[a] Gamma[b1] Pochhammer[c2,m])+(y^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{b1,a+m},{1+a+b1-c1+m},1-x] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])},{Abs[(1-y)/x]<1&&Abs[x]<1,(x^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{b2,a+m},{1+a+b2-c2+m},1-y] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])+((1-y)^(-a-b2+c2) ((1-y)/x)^m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a+b1-b2+c2] HypergeometricPFQ[{-b2+c2,1-b2-m,-a+b1-b2+c2-m},{1-a-b2+c2-m,-a-b2+c1+c2-m},x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a+b2-c2,m] Pochhammer[1+a+b2-c1-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b2+c1+c2] Pochhammer[1+a-b1+b2-c2,m])+((1-y)^(-a-b2+c2) ((1-y)/x)^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] HypergeometricPFQ[{-b2+c2,-a+b1+c2+m},{1-a+b1-b2+c2+m},1-y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b1,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m])},{Abs[(1-x)/y]<1&&Abs[y]<1,((1-x)^(-a-b1+c1) ((1-x)/y)^m Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] HypergeometricPFQ[{-b1+c1,-a+b2+c1+m},{1-a-b1+b2+c1+m},1-x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m])+((1-x)^(-a-b1+c1) ((1-x)/y)^m Gamma[a+b1-c1] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] HypergeometricPFQ[{-b1+c1,1-b1-m,-a-b1+b2+c1-m},{1-a-b1+c1-m,-a-b1+c1+c2-m},y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a+b1-c1,m] Pochhammer[1+a+b1-c1-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b1+c1+c2] Pochhammer[1+a+b1-b2-c1,m])+(y^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{b1,a+m},{1+a+b1-c1+m},1-x] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])},{1/Abs[x]<1&&Abs[1-y]<1,((1/x)^m (-x)^-a Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{b2,a+m,1+a-c1+m},{1+a-b1+m,c2+m},(1-y)/x] Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])+((1/x)^m (-x)^-b1 Gamma[a-b1] Gamma[c1] Gamma[c2] Gamma[-a+b1-b2+c2] HypergeometricPFQ[{b2,a-b1-m},{1+a-b1+b2-c2-m},1-y] Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1-b2+c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-a+b1+c2] Gamma[-b2+c2] Pochhammer[1-a+b1,m] Pochhammer[-a+b1+c2,m])+((-x)^-b1 (1-y)^(-a+b1-b2+c2) ((1-y)/x)^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] HypergeometricPFQ[{-b2+c2,-a+b1+c2+m},{1-a+b1-b2+c2+m},1-y] Pochhammer[b1,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m])},{1/Abs[y]<1&&Abs[1-x]<1,((1/y)^m (-y)^-a Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{b1,a+m,1+a-c2+m},{1+a-b2+m,c1+m},(1-x)/y] Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])+((1/y)^m (-y)^-b2 Gamma[a-b2] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] HypergeometricPFQ[{b1,a-b2-m},{1+a+b1-b2-c1-m},1-x] Pochhammer[b2,m] Pochhammer[-a-b1+b2+c1,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1-a+b2,m] Pochhammer[-a+b2+c1,m])+((1-x)^(-a-b1+b2+c1) ((1-x)/y)^m (-y)^-b2 Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] HypergeometricPFQ[{-b1+c1,-a+b2+c1+m},{1-a-b1+b2+c1+m},1-x] Pochhammer[b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m])},{Abs[x]<1&&1/Abs[y]<1&&(1+Abs[x])/Abs[y]<1,(x^m (-y)^-b2 Gamma[a-b2] Gamma[c2] HypergeometricPFQ[{b2,1+b2-c2},{1-a+b2-m},1/y] Pochhammer[b1,m] Pochhammer[a-b2,m])/(m! Gamma[a] Gamma[-b2+c2] Pochhammer[c1,m])+((-(x/y))^m (-y)^-a Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a-b2+m},1/y] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{Abs[y]<1&&1/Abs[x]<1&&(1+Abs[y])/Abs[x]<1,((-x)^-b1 y^m Gamma[a-b1] Gamma[c1] HypergeometricPFQ[{b1,1+b1-c1},{1-a+b1-m},1/x] Pochhammer[a-b1,m] Pochhammer[b2,m])/(m! Gamma[a] Gamma[-b1+c1] Pochhammer[c2,m])+((-x)^-a (-(y/x))^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a-b1+m},1/x] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])},{1/Abs[x]<1&&Abs[x/y]<1&&Abs[x/y]<Abs[x]/(1+Abs[x])&&1/Abs[y]<1&&Abs[x/y]+1/Abs[y]<1,(x^m (-y)^(-a-m) Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a-b2+m},1/y] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])+((-x)^-b1 (-y)^-b2 y^-m Gamma[a-b1-b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b1,1+b1-c1},{1-a+b1+b2+m},1/x] Pochhammer[b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m])+((-1)^m (-x)^(-a+b2) x^m (-y)^-b2 y^-m Gamma[a-b2] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{a-b2-m,1+a-b2-c1-m},{1+a-b1-b2-m},1/x] Pochhammer[b2,m] Pochhammer[-a+b1+b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1-a+b2,m] Pochhammer[-a+b2+c1,m])},{1/Abs[y]<1&&Abs[y/x]<1&&Abs[y/x]<Abs[y]/(1+Abs[y])&&1/Abs[x]<1&&1/Abs[x]+Abs[y/x]<1,((-x)^-b1 (-y)^-b2 y^-m Gamma[a-b1-b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b1,1+b1-c1},{1-a+b1+b2+m},1/x] Pochhammer[b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m])+((-x)^(-a-m) y^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a-b1+m},1/x] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])+((-1)^m (-x)^-b1 x^-m (-y)^(-a+b1) y^m Gamma[a-b1] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{a-b1-m,1+a-b1-c2-m},{1+a-b1-b2-m},1/y] Pochhammer[b1,m] Pochhammer[-a+b1+b2,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a+b1+c2] Pochhammer[1-a+b1,m] Pochhammer[-a+b1+c2,m])},{Abs[x/(-1+x)]+Abs[y/(1-x)]<1&&Abs[x/(-1+x)]<1,((1-x)^-a (x/(-1+x))^m HypergeometricPFQ[{b2,a+m},{c2},y/(1-x)] Pochhammer[a,m] Pochhammer[-b1+c1,m])/(m! Pochhammer[c1,m])},{Abs[x/(-1+x+y)]<1&&Abs[(-1+x+y)/(-1+x)]<1,((1-x)^-a (x/(-1+x+y))^m Gamma[a+b2-c2] Gamma[c2] HypergeometricPFQ[{-b2+c2,-a+c2-m},{1-a-b2+c2-m},(-1+x+y)/(-1+x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[-b1+c1,m] Pochhammer[a+b2-c2,m])/(m! Gamma[a] Gamma[b2] Pochhammer[c1,m])+((1-x)^-a (x/(-1+x))^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{b2,a+m},{1+a+b2-c2+m},(-1+x+y)/(-1+x)] Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])},{Abs[1-x]>1&&Abs[y]<1,((1-x)^-b1 y^m Gamma[a-b1] Gamma[c1] HypergeometricPFQ[{b1,-a+c1-m},{1-a+b1-m},1/(1-x)] Pochhammer[a-b1,m] Pochhammer[b2,m])/(m! Gamma[a] Gamma[-b1+c1] Pochhammer[c2,m])+((1-x)^-a (y/(1-x))^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{-b1+c1,a+m},{1+a-b1+m},1/(1-x)] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])},{Abs[(-1+x+y)/x]<1&&Abs[x/(-1+x)]<1,((1-x)^-a (x/(-1+x))^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{b2,a+m},{1+a+b2-c2+m},(-1+x+y)/(-1+x)] Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])+((1-x)^-a ((-1+x+y)/x)^m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a-b1-b2+c1+c2] HypergeometricPFQ[{-b2+c2,1-b2-m,-a-b1-b2+c1+c2-m},{1-a-b2+c2-m,-a-b2+c1+c2-m},x/(-1+x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a+b2-c2,m] Pochhammer[1+a+b2-c1-c2,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a-b2+c1+c2] Pochhammer[1+a+b1+b2-c1-c2,m])+((1-x)^-a ((-1+x+y)/x)^m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] HypergeometricPFQ[{-b2+c2,-a-b1+c1+c2+m},{1-a-b1-b2+c1+c2+m},(-1+x+y)/(-1+x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m])},{Abs[y]>1&&Abs[y/(1-x)]<1,((1-x)^-b1 (1/y)^m (-y)^(-a+b1) Gamma[a-b1] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b1,-a+b1+b2-m,1+b1-c1-m},{1-a+b1-m,-a+b1+c2-m},y/(1-x)] Pochhammer[a-b1,m] Pochhammer[1+a-b1-c2,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a+b1+c2] Pochhammer[1+a-b1-b2,m])+((1-x)^-b1 (1/y)^m (-y)^-b2 Gamma[a-b1-b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b1,-a+b2+c1+m},{1-a+b1+b2+m},1/(1-x)] Pochhammer[b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m])+((1-x)^-a (y/(1-x))^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{-b1+c1,a+m},{1+a-b1+m},1/(1-x)] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])},{Abs[(-1+x)/x]<1&&Abs[(-1+x+y)/(-1+x)]<1,((-1+1/x)^a (1-x)^-a ((-1+x)/x)^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{b2,a+m,1+a-c1+m},{1+a+b1-c1+m,c2+m},(-1+x+y)/x] Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])+((-1+1/x)^(-b1+c1) (1-x)^-a ((-1+x)/x)^m Gamma[a+b1-c1] Gamma[c1] Gamma[c2] Gamma[-a-b1-b2+c1+c2] HypergeometricPFQ[{b2,a+b1-c1-m},{1+a+b1+b2-c1-c2-m},(-1+x+y)/(-1+x)] Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-a-b1-b2+c1+c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Gamma[-a-b1+c1+c2] Pochhammer[1-a-b1+c1,m] Pochhammer[-a-b1+c1+c2,m])+((-1+1/x)^(-b1+c1) (1-x)^-a ((-1+x+y)/x)^m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] HypergeometricPFQ[{-b2+c2,-a-b1+c1+c2+m},{1-a-b1-b2+c1+c2+m},(-1+x+y)/(-1+x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b1", "+", "b2", "-", "c1", "-", "c2"}]], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "-", "b2", "+", "c1", "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m])},{Abs[1-x]>1&&Abs[(1-x)/y]<1,((1-x)^-a ((1-x)/y)^m Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{-b1+c1,a+m,1+a-c2+m},{1+a-b2+m,c1+m},1/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "a"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "a"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])+((1-x)^(-b1-b2) (1/y)^m Gamma[a-b1-b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b1,-a+b2+c1+m},{1-a+b1+b2+m},1/(1-x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m])+((1-x)^-a ((1-x)/y)^m Gamma[a-b2] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{-b1+c1,a-b2-m},{1+a-b1-b2-m},1/(1-x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b2,m] Pochhammer[-a+b1+b2,m] Pochhammer[1+b2-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1-a+b2,m] Pochhammer[-a+b2+c1,m])},{Abs[x/(-1+x)]<1&&Abs[(1-x)/y]<1&&(1+Abs[x/(-1+x)]) Abs[(1-x)/y]<1,((1-x)^-a (x/(-1+x))^m Gamma[a-b2] Gamma[c2] HypergeometricPFQ[{b2,1+b2-c2},{1-a+b2-m},-((-1+x)/y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b2,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[-b2+c2] Pochhammer[c1,m])+((1-x)^-a (x/y)^m Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a-b2+m},(1-x)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "a"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "a"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{Abs[y/(1-x)]<1&&Abs[(-1+x)/x]<1&&Abs[(-1+x)/x] (1+Abs[y/(1-x)])<1,((-1+1/x)^(-b1+c1) (1-x)^-a (y/(1-x))^m Gamma[a+b1-c1] Gamma[c1] HypergeometricPFQ[{1-b1,-b1+c1},{1-a-b1+c1-m},1-1/x] Pochhammer[b2,m] Pochhammer[a+b1-c1,m])/(m! Gamma[a] Gamma[b1] Pochhammer[c2,m])+((-1+1/x)^a (1-x)^-a (y/x)^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a+b1-c1+m},(-1+x)/x] Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])},{1/Abs[x/(-1+x)]<1&&Abs[x/y]<1&&Abs[x/(-1+x)] (-1+Abs[x/y])+Abs[x/y]<0&&1/Abs[y/(1-x)]<1&&Abs[x/y]+1/Abs[y/(1-x)]<1,((-1+1/x)^(a-b2) (1-x)^-a (x/(-1+x))^-m Gamma[a-b2] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] HypergeometricPFQ[{b2,1+b2-c2,-a-b1+b2+c1-m},{1-a+b2-m,-a+b2+c1-m},-(((1-x) x)/((-1+x) y))] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b2,m] Pochhammer[1+a-b2-c1,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a+b1-b2-c1,m])+((-1+1/x)^(-b1+c1) (1-x)^-a (x/(-1+x))^-m Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,1+b2-c2},{1-a-b1+b2+c1+m},(1-x)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m])+((1-x)^-a (x/(-1+x))^m (y/(-1+x))^-m Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a-b2+m},(1-x)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "a"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "a"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{1/Abs[y/(1-x)]<1&&Abs[y/x]<1&&Abs[y/(1-x)] (-1+Abs[y/x])+Abs[y/x]<0&&1/Abs[x/(-1+x)]<1&&1/Abs[x/(-1+x)]+Abs[y/x]<1,((-1+1/x)^a (1-x)^-a (x/(-1+x))^-m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{b2,a+m,1+a-c1+m},{c2,1+a+b1-c1+m},-(((-1+x) y)/((1-x) x))] Pochhammer[a,m] Pochhammer[1+a-c1,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m])+((-1+1/x)^(-b1+c1) (1-x)^-a (x/(-1+x))^-m Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,1+b2-c2},{1-a-b1+b2+c1+m},(1-x)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], "b2"], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m])+((-1)^m (-1+1/x)^(-b1+c1) (1-x)^-a (x/(-1+x))^-m (y/(1-x))^m Gamma[a+b1-c1] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] HypergeometricPFQ[{a+b1-c1-m,1+a+b1-c1-c2-m},{1+a+b1-b2-c1-m},(1-x)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], "y"], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
RowBox[{"x", ">=", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-a-b1+b2+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b1+c1+c2] Pochhammer[1-a-b1+c1,m] Pochhammer[-a-b1+c1+c2,m])},{Abs[x/(1-y)]+Abs[y/(-1+y)]<1&&Abs[x/(1-y)]<1,((x/(1-y))^m (1-y)^-a HypergeometricPFQ[{-b2+c2,a+m},{c2},y/(-1+y)] Pochhammer[a,m] Pochhammer[b1,m])/(m! Pochhammer[c1,m])},{Abs[x]<1&&Abs[1-y]>1,(x^m (1-y)^-b2 Gamma[a-b2] Gamma[c2] HypergeometricPFQ[{b2,-a+c2-m},{1-a+b2-m},1/(1-y)] Pochhammer[b1,m] Pochhammer[a-b2,m])/(m! Gamma[a] Gamma[-b2+c2] Pochhammer[c1,m])+((x/(1-y))^m (1-y)^-a Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{-b2+c2,a+m},{1+a-b2+m},1/(1-y)] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{Abs[y/(-1+x+y)]<1&&Abs[(-1+x+y)/(-1+y)]<1,((1-y)^-a (y/(-1+x+y))^m Gamma[a+b1-c1] Gamma[c1] HypergeometricPFQ[{-b1+c1,-a+c1-m},{1-a-b1+c1-m},(-1+x+y)/(-1+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<=", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a+b1-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[b1] Pochhammer[c2,m])+((1-y)^-a (y/(-1+y))^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{b1,a+m},{1+a+b1-c1+m},(-1+x+y)/(-1+y)] Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])},{Abs[x]>1&&Abs[x/(1-y)]<1,((1/x)^m (-x)^-b1 (1-y)^-b2 Gamma[a-b1-b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,-a+b1+c2+m},{1-a+b1+b2+m},1/(1-y)] Pochhammer[b1,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m])+((1/x)^m (-x)^(-a+b2) (1-y)^-b2 Gamma[a-b2] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,-a+b1+b2-m,1+b2-c2-m},{1-a+b2-m,-a+b2+c1-m},x/(1-y)] Pochhammer[a-b2,m] Pochhammer[1+a-b2-c1,m])/(m! Gamma[a] Gamma[b1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a-b1-b2,m])+((x/(1-y))^m (1-y)^-a Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{-b2+c2,a+m},{1+a-b2+m},1/(1-y)] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{Abs[(-1+x+y)/y]<1&&Abs[y/(-1+y)]<1,((1-y)^-a ((-1+x+y)/y)^m Gamma[a+b1-c1] Gamma[c1] Gamma[c2] Gamma[-a-b1-b2+c1+c2] HypergeometricPFQ[{-b1+c1,1-b1-m,-a-b1-b2+c1+c2-m},{1-a-b1+c1-m,-a-b1+c1+c2-m},y/(-1+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<=", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a+b1-c1,m] Pochhammer[1+a+b1-c1-c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Gamma[-a-b1+c1+c2] Pochhammer[1+a+b1+b2-c1-c2,m])+((1-y)^-a (y/(-1+y))^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{b1,a+m},{1+a+b1-c1+m},(-1+x+y)/(-1+y)] Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])+((1-y)^-a ((-1+x+y)/y)^m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] HypergeometricPFQ[{-b1+c1,-a-b2+c1+c2+m},{1-a-b1-b2+c1+c2+m},(-1+x+y)/(-1+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<=", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m])},{Abs[1-y]>1&&Abs[(1-y)/x]<1,((1/x)^m (1-y)^(-b1-b2) Gamma[a-b1-b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,-a+b1+c2+m},{1-a+b1+b2+m},1/(1-y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b1,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m])+((1-y)^-a ((1-y)/x)^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{-b2+c2,a+m,1+a-c1+m},{1+a-b1+m,c2+m},1/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])+((1-y)^-a ((1-y)/x)^m Gamma[a-b1] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] HypergeometricPFQ[{-b2+c2,a-b1-m},{1+a-b1-b2-m},1/(1-y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b1,m] Pochhammer[-a+b1+b2,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a+b1+c2] Pochhammer[1-a+b1,m] Pochhammer[-a+b1+c2,m])},{Abs[(-1+y)/y]<1&&Abs[(-1+x+y)/(-1+y)]<1,((-1+1/y)^a (1-y)^-a ((-1+y)/y)^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{b1,a+m,1+a-c2+m},{c1+m,1+a+b2-c2+m},(-1+x+y)/y] Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])+((-1+1/y)^(-b2+c2) (1-y)^-a ((-1+y)/y)^m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a-b1-b2+c1+c2] HypergeometricPFQ[{b1,a+b2-c2-m},{1+a+b1+b2-c1-c2-m},(-1+x+y)/(-1+y)] Pochhammer[1-b2,m] Pochhammer[-b2+c2,m] Pochhammer[-a-b1-b2+c1+c2,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a-b2+c1+c2] Pochhammer[1-a-b2+c2,m] Pochhammer[-a-b2+c1+c2,m])+((-1+1/y)^(-b2+c2) (1-y)^-a ((-1+x+y)/y)^m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] HypergeometricPFQ[{-b1+c1,-a-b2+c1+c2+m},{1-a-b1-b2+c1+c2+m},(-1+x+y)/(-1+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{"a", "+", "b1", "+", "b2", "-", "c1", "-", "c2"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<=", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "-", "b2", "+", "c1", "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m])},{Abs[x/(1-y)]<1&&Abs[(-1+y)/y]<1&&(1+Abs[x/(1-y)]) Abs[(-1+y)/y]<1,((-1+1/y)^(-b2+c2) (x/(1-y))^m (1-y)^-a Gamma[a+b2-c2] Gamma[c2] HypergeometricPFQ[{1-b2,-b2+c2},{1-a-b2+c2-m},1-1/y] Pochhammer[b1,m] Pochhammer[a+b2-c2,m])/(m! Gamma[a] Gamma[b2] Pochhammer[c1,m])+((-1+1/y)^a (1-y)^-a (x/y)^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a+b2-c2+m},(-1+y)/y] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])},{Abs[y/(-1+y)]<1&&Abs[(1-y)/x]<1&&Abs[(1-y)/x] (1+Abs[y/(-1+y)])<1,((1-y)^-a (y/(-1+y))^m Gamma[a-b1] Gamma[c1] HypergeometricPFQ[{b1,1+b1-c1},{1-a+b1-m},-((-1+y)/x)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b1,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[-b1+c1] Pochhammer[c2,m])+((1-y)^-a (y/x)^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a-b1+m},(1-y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])},{1/Abs[x/(1-y)]<1&&Abs[x/y]<1&&Abs[x/(1-y)] (-1+Abs[x/y])+Abs[x/y]<0&&1/Abs[y/(-1+y)]<1&&Abs[x/y]+1/Abs[y/(-1+y)]<1,{((-1+1/y)^a (x/(1-y))^m (1-y)^-a (-(y/(-1+y)))^-m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a+b2-c2+m},(-1+y)/y] Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m]),((-1)^m (-1+1/y)^(-b2+c2) (x/(1-y))^m (1-y)^-a (y/(-1+y))^-m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a+b1-b2+c2] HypergeometricPFQ[{a+b2-c2-m,1+a+b2-c1-c2-m},{1+a-b1+b2-c2-m},(1-y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m] Pochhammer[-a+b1-b2+c2,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b2+c1+c2] Pochhammer[1-a-b2+c2,m] Pochhammer[-a-b2+c1+c2,m]),((-1+1/y)^(-b2+c2) (x/(1-y))^-m (1-y)^-a Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] HypergeometricPFQ[{1-b2,-b2+c2},{1-a+b1-b2+c2+m},(-1+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b1,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m])}},{1/Abs[y/(-1+y)]<1&&Abs[y/x]<1&&Abs[y/x] (1+Abs[y/(-1+y)])<Abs[y/(-1+y)]&&1/Abs[x/(1-y)]<1&&1/Abs[x/(1-y)]+Abs[y/x]<1,((1-y)^-a (x/(-1+y))^-m (y/(-1+y))^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a-b1+m},(1-y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])+((-1)^m (-1+1/y)^(a-b1) (x/(1-y))^-m (1-y)^-a (y/(-1+y))^m Gamma[a-b1] Gamma[c1] Gamma[c2] Gamma[-a+b1-b2+c2] HypergeometricPFQ[{a-b1-m,1+a-b1-c2-m},{1+a-b1+b2-c2-m},(-1+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1-b2+c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-a+b1+c2] Gamma[-b2+c2] Pochhammer[1-a+b1,m] Pochhammer[-a+b1+c2,m])+((-1+1/y)^(-b2+c2) (x/(1-y))^-m (1-y)^-a Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] HypergeometricPFQ[{1-b2,-b2+c2},{1-a+b1-b2+c2+m},(-1+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"-", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], ">", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], "x"], ")"}], "b1"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[b1,m] Pochhammer[1+b1-c1,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m])},{Abs[x/(-1+x+y)]+Abs[y/(-1+x+y)]<1&&Abs[x/(-1+x+y)]<1,((1-x-y)^-a (x/(-1+x+y))^m HypergeometricPFQ[{-b2+c2,a+m},{c2},y/(-1+x+y)] Pochhammer[a,m] Pochhammer[-b1+c1,m])/(m! Pochhammer[c1,m])},{Abs[x/(-1+x)]<1&&Abs[(-1+x)/(-1+x+y)]<1,((x/(-1+x))^m (1-x-y)^-a Gamma[a-b2] Gamma[c2] HypergeometricPFQ[{b2,-a+c2-m},{1-a+b2-m},(-1+x)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"a", "-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b2,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[-b2+c2] Pochhammer[c1,m])+((1-x-y)^-a (x/(-1+x+y))^m Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{-b2+c2,a+m},{1+a-b2+m},(-1+x)/(-1+x+y)] Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{Abs[y/(-1+y)]<1&&Abs[(-1+y)/(-1+x+y)]<1,((1-x-y)^-a (y/(-1+y))^m Gamma[a-b1] Gamma[c1] HypergeometricPFQ[{b1,-a+c1-m},{1-a+b1-m},(-1+y)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"a", "-", "b1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b1,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[-b1+c1] Pochhammer[c2,m])+((1-x-y)^-a (y/(-1+x+y))^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{-b1+c1,a+m},{1+a-b1+m},(-1+y)/(-1+x+y)] Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])},{Abs[(-1+x)/x]<1&&Abs[x/(-1+x+y)]<1,((-1+1/x)^(a-b2) ((-1+x)/x)^m (1-x-y)^-a Gamma[a-b2] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] HypergeometricPFQ[{b2,-a-b1+b2+c1-m,1+b2-c2-m},{1-a+b2-m,-a+b2+c1-m},x/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"a", "-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b2,m] Pochhammer[1+a-b2-c1,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a+b1-b2-c1,m])+((-1+1/x)^(-b1+c1) ((-1+x)/x)^m (1-x-y)^-a Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,-a-b1+c1+c2+m},{1-a-b1+b2+c1+m},(-1+x)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"a", "-", "b2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m])+((1-x-y)^-a (x/(-1+x+y))^m Gamma[-a+b2] Gamma[c2] HypergeometricPFQ[{-b2+c2,a+m},{1+a-b2+m},(-1+x)/(-1+x+y)] Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m] Pochhammer[c1,m])},{Abs[(-1+y)/y]<1&&Abs[y/(-1+x+y)]<1,((-1+1/y)^(a-b1) (1-x-y)^-a ((-1+y)/y)^m Gamma[a-b1] Gamma[c1] Gamma[c2] Gamma[-a+b1-b2+c2] HypergeometricPFQ[{b1,1+b1-c1-m,-a+b1-b2+c2-m},{1-a+b1-m,-a+b1+c2-m},y/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"a", "-", "b1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a-b1,m] Pochhammer[1+a-b1-c2,m])/(m! Gamma[a] Gamma[-b1+c1] Gamma[-a+b1+c2] Gamma[-b2+c2] Pochhammer[1+a-b1+b2-c2,m])+((1-x-y)^-a (y/(-1+x+y))^m Gamma[-a+b1] Gamma[c1] HypergeometricPFQ[{-b1+c1,a+m},{1+a-b1+m},(-1+y)/(-1+x+y)] Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m] Pochhammer[c2,m])+((-1+1/y)^(-b2+c2) (1-x-y)^-a ((-1+y)/y)^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] HypergeometricPFQ[{b1,-a-b2+c1+c2+m},{1-a+b1-b2+c2+m},(-1+y)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"a", "-", "b1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m])},{Abs[(-1+x+y)/x]<1&&Abs[(-1+x)/(-1+x+y)]<1,(((-1+x)/x)^m (1-x-y)^-a Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] HypergeometricPFQ[{b2,-a-b1+c1+c2+m},{1-a-b1+b2+c1+m},(-1+x)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "b2", "+", "c1"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x"}]], ")"}], 
RowBox[{"a", "+", "b1", "-", "b2", "-", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m])+((1-x-y)^-a ((-1+x+y)/x)^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{-b2+c2,a+m,1+a-c1+m},{1+a+b1-c1+m,c2+m},(-1+x)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[b2,m] Pochhammer[1+a-c1,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])+((1-x-y)^-a ((-1+x+y)/x)^m Gamma[a+b1-c1] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] HypergeometricPFQ[{-b2+c2,a+b1-c1-m},{1+a+b1-b2-c1-m},(-1+x)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-a-b1+b2+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b1+c1+c2] Pochhammer[1-a-b1+c1,m] Pochhammer[-a-b1+c1+c2,m])},{Abs[(-1+x+y)/y]<1&&Abs[(-1+y)/(-1+x+y)]<1,((1-x-y)^-a ((-1+x+y)/y)^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{-b1+c1,a+m,1+a-c2+m},{c1+m,1+a+b2-c2+m},(-1+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[b1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])+((1-x-y)^-a ((-1+y)/y)^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] HypergeometricPFQ[{b1,-a-b2+c1+c2+m},{1-a+b1-b2+c2+m},(-1+y)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "+", "b1", "-", "b2", "+", "c2"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], 
RowBox[{
RowBox[{"-", "1"}], "+", "y"}]], ")"}], 
RowBox[{"a", "-", "b1", "+", "b2", "-", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m])+((1-x-y)^-a ((-1+x+y)/y)^m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a+b1-b2+c2] HypergeometricPFQ[{-b1+c1,a+b2-c2-m},{1+a-b1+b2-c2-m},(-1+y)/(-1+x+y)] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m] Pochhammer[-a+b1-b2+c2,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b2+c1+c2] Pochhammer[1-a-b2+c2,m] Pochhammer[-a-b2+c1+c2,m])},{Abs[x/(-1+x+y)]<1&&Abs[(-1+x+y)/y]<1&&(1+Abs[x/(-1+x+y)]) Abs[(-1+x+y)/y]<1,((1-x-y)^-a (x/(-1+x+y))^m Gamma[a+b2-c2] Gamma[c2] HypergeometricPFQ[{1-b2,-b2+c2},{1-a-b2+c2-m},(-1+x+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[-b1+c1,m] Pochhammer[a+b2-c2,m])/(m! Gamma[a] Gamma[b2] Pochhammer[c1,m])+((1-x-y)^-a (-(x/y))^m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a+b2-c2+m},(-1+x+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])},{Abs[y/(-1+x+y)]<1&&Abs[(-1+x+y)/x]<1&&(1+Abs[y/(-1+x+y)]) Abs[(-1+x+y)/x]<1,((1-x-y)^-a (y/(-1+x+y))^m Gamma[a+b1-c1] Gamma[c1] HypergeometricPFQ[{1-b1,-b1+c1},{1-a-b1+c1-m},(-1+x+y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a+b1-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[a] Gamma[b1] Pochhammer[c2,m])+((1-x-y)^-a (-(y/x))^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a+b1-c1+m},(-1+x+y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])},{1/Abs[x/(-1+x+y)]<1&&Abs[x/y]<1&&Abs[x/y] (1+Abs[x/(-1+x+y)])<Abs[x/(-1+x+y)]&&1/Abs[y/(-1+x+y)]<1&&Abs[x/y]+1/Abs[y/(-1+x+y)]<1,((1-x-y)^-a (x/(-1+x+y))^m (-(y/(-1+x+y)))^-m Gamma[c2] Gamma[-a-b2+c2] HypergeometricPFQ[{a+m,1+a-c2+m},{1+a+b2-c2+m},(-1+x+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m])+((-1)^m (1-x-y)^-a (x/(-1+x+y))^m (y/(-1+x+y))^-m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a-b1-b2+c1+c2] HypergeometricPFQ[{a+b2-c2-m,1+a+b2-c1-c2-m},{1+a+b1+b2-c1-c2-m},(-1+x+y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b2", "+", "c2"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{"a", "+", "b2", "-", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b2,m] Pochhammer[-b2+c2,m] Pochhammer[-a-b1-b2+c1+c2,m])/(m! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a-b2+c1+c2] Pochhammer[1-a-b2+c2,m] Pochhammer[-a-b2+c1+c2,m])+((1-x-y)^-a (x/(-1+x+y))^-m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] HypergeometricPFQ[{1-b2,-b2+c2},{1-a-b1-b2+c1+c2+m},(-1+x+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m])},{1/Abs[y/(-1+x+y)]<1&&Abs[y/x]<1&&Abs[y/x] (1+Abs[y/(-1+x+y)])<Abs[y/(-1+x+y)]&&1/Abs[x/(-1+x+y)]<1&&Abs[y/x]+1/Abs[x/(-1+x+y)]<1,((1-x-y)^-a (-(x/(-1+x+y)))^-m (y/(-1+x+y))^m Gamma[c1] Gamma[-a-b1+c1] HypergeometricPFQ[{a+m,1+a-c1+m},{1+a+b1-c1+m},(-1+x+y)/x] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"-", "a"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], "a"], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[a,m] Pochhammer[1+a-c1,m] Pochhammer[-b2+c2,m])/(m! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m] Pochhammer[c2,m])+((-1)^m (1-x-y)^-a (x/(-1+x+y))^-m (y/(-1+x+y))^m Gamma[a+b1-c1] Gamma[c1] Gamma[c2] Gamma[-a-b1-b2+c1+c2] HypergeometricPFQ[{a+b1-c1-m,1+a+b1-c1-c2-m},{1+a+b1+b2-c1-c2-m},(-1+x+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{
RowBox[{"-", "a"}], "-", "b1", "+", "c1"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{"a", "+", "b1", "-", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-a-b1-b2+c1+c2,m])/(m! Gamma[a] Gamma[b1] Gamma[-b2+c2] Gamma[-a-b1+c1+c2] Pochhammer[1-a-b1+c1,m] Pochhammer[-a-b1+c1+c2,m])+((1-x-y)^-a (x/(-1+x+y))^-m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] HypergeometricPFQ[{1-b2,-b2+c2},{1-a-b1-b2+c1+c2+m},(-1+x+y)/y] (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["x", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b1", "-", "c1"}]], 
RowBox[{
RowBox[{"1", "+", "x"}], "<", "y"}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "x"]}], ")"}], 
RowBox[{
RowBox[{"-", "b1"}], "+", "c1"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) (\!\(\*
TagBox[GridBox[{
{"\[Piecewise]", GridBox[{
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox["y", 
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}]]}], ")"}], 
RowBox[{"b2", "-", "c2"}]], 
RowBox[{"x", ">", 
RowBox[{"1", "+", "y"}]}]},
{
SuperscriptBox[
RowBox[{"(", 
RowBox[{"-", 
FractionBox[
RowBox[{
RowBox[{"-", "1"}], "+", "x", "+", "y"}], "y"]}], ")"}], 
RowBox[{
RowBox[{"-", "b2"}], "+", "c2"}]], 
TagBox["True",
"PiecewiseDefault",
AutoDelete->True]}
},
AllowedDimensions->{2, Automatic},
Editable->True,
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.84]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}},
Selectable->True]}
},
GridBoxAlignment->{"Columns" -> {{Left}}, "Rows" -> {{Baseline}}},
GridBoxItemSize->{"Columns" -> {{Automatic}}, "Rows" -> {{1.}}},
GridBoxSpacings->{"Columns" -> {Offset[0.27999999999999997`], {Offset[0.35]}, Offset[0.27999999999999997`]}, "Rows" -> {Offset[0.2], {Offset[0.4]}, Offset[0.2]}}],
"Piecewise",
DeleteWithContents->True,
Editable->False,
SelectWithContents->True,
Selectable->False,
StripWrapperBoxes->True]\)) Pochhammer[1-b1,m] Pochhammer[-b1+c1,m])/(m! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m])}};


F2serieswsum[a1_,a2_,b1_,b2_,c_,x_,y_,m_,n_]:= {{Abs[x]<1&&Abs[y]<1,(x^m y^n Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n])/(m! n! Pochhammer[c,m+n])},{Abs[1-x]<1&&Abs[y]<1&&Abs[y]<-1+1/Abs[1-x]&&Abs[1-x]<1&&Abs[(-1+x) y]<1&&Abs[(-1+x) y]<1-Abs[1-x],((-1)^n (1-x)^m y^n Gamma[c] Gamma[-a1-b1+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n])/(m! n! Gamma[-a1+c] Gamma[-b1+c] Pochhammer[1+a1+b1-c,m-n] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n])+((-1)^n (1-x)^(-a1-b1+c+m+n) y^n Gamma[a1+b1-c] Gamma[c] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[-a1+c,m+n] Pochhammer[-b1+c,m+n])/(m! n! Gamma[a1] Gamma[b1] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n] Pochhammer[1-a1-b1+c,m+n])},{Abs[1-y]<1&&Abs[x]<1&&Abs[x]<-1+1/Abs[1-y]&&Abs[1-y]<1&&Abs[x (-1+y)]<1&&Abs[x (-1+y)]<1-Abs[1-y],((-1)^m x^m (1-y)^n Gamma[c] Gamma[-a2-b2+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n])/(m! n! Gamma[-a2+c] Gamma[-b2+c] Pochhammer[1+a2+b2-c,-m+n] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m])+((-1)^m x^m (1-y)^(-a2-b2+c+m+n) Gamma[a2+b2-c] Gamma[c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[-a2+c,m+n] Pochhammer[-b2+c,m+n])/(m! n! Gamma[a2] Gamma[b2] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m] Pochhammer[1-a2-b2+c,m+n])},{1/Abs[x]<1&&Abs[y]<1&&Abs[y]<1/(1+1/Abs[x])&&1/Abs[x]<1&&Abs[y]<1&&Abs[y]<1/(1+1/Abs[x]),((-1)^n (-x)^-a1 x^-m y^n Gamma[-a1+b1] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[1+a1-c,m-n])/(m! n! Gamma[b1] Gamma[-a1+c] Pochhammer[1+a1-b1,m])+((-1)^n (-x)^-b1 x^-m y^n Gamma[a1-b1] Gamma[c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c,m-n])/(m! n! Gamma[a1] Gamma[-b1+c] Pochhammer[1-a1+b1,m])},{1/Abs[y]<1&&Abs[x]<1&&Abs[x]<1/(1+1/Abs[y])&&1/Abs[y]<1&&Abs[x]<1&&Abs[x]<1/(1+1/Abs[y]),((-1)^m x^m (-y)^-a2 y^-n Gamma[-a2+b2] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a2-c,-m+n])/(m! n! Gamma[b2] Gamma[-a2+c] Pochhammer[1+a2-b2,n])+((-1)^m x^m (-y)^-b2 y^-n Gamma[a2-b2] Gamma[c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b2-c,-m+n])/(m! n! Gamma[a2] Gamma[-b2+c] Pochhammer[1-a2+b2,n])},{1/Abs[x]<1&&1/Abs[y]<1&&1/Abs[y]<1-1/Abs[x],((-x)^-a1 x^-m (-y)^-a2 y^-n Gamma[-a1+b1] Gamma[-a2+b2] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,m+n])/(m! n! Gamma[b1] Gamma[b2] Gamma[-a1-a2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a2-b2,n])+((-x)^-b1 x^-m (-y)^-a2 y^-n Gamma[a1-b1] Gamma[-a2+b2] Gamma[c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a2+b1-c,m+n])/(m! n! Gamma[a1] Gamma[b2] Gamma[-a2-b1+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2-b2,n])+((-x)^-a1 x^-m (-y)^-b2 y^-n Gamma[-a1+b1] Gamma[a2-b2] Gamma[c] Pochhammer[a1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,m+n])/(m! n! Gamma[a2] Gamma[b1] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1-a2+b2,n])+((-x)^-b1 x^-m (-y)^-b2 y^-n Gamma[a1-b1] Gamma[a2-b2] Gamma[c] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1+b2-c,m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1-a2+b2,n])},{Abs[1-x]<1&&1/Abs[y]<1&&1+1/Abs[y]<1/Abs[y-x y]&&Abs[(-1+x) y]<1&&Abs[1-x]+Abs[(-1+x) y]<1&&Abs[y-x y]<1,((1-x)^m (-y)^-a2 y^-n Gamma[-a2+b2] Gamma[c] Gamma[-a1-a2-b1+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a1+a2-c,n] Pochhammer[1+a2+b1-c,n])/(m! n! Gamma[b2] Gamma[-a1-a2+c] Gamma[-a2-b1+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2+b1-c,m+n])+((1-x)^m (-y)^-b2 y^-n Gamma[a2-b2] Gamma[c] Gamma[-a1-b1-b2+c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,n] Pochhammer[1+b1+b2-c,n])/(m! n! Gamma[a2] Gamma[-a1-b2+c] Gamma[-b1-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b1+b2-c,m+n])+((-1)^n (1-x)^(-a1-b1+c+m+n) y^n Gamma[a1+b1-c] Gamma[c] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[-a1+c,m+n] Pochhammer[-b1+c,m+n])/(m! n! Gamma[a1] Gamma[b1] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n] Pochhammer[1-a1-b1+c,m+n])+((-1)^m (1-x)^m (-y)^(a1+b1-c+m) y^-n Gamma[a1+a2+b1-c] Gamma[a1+b1+b2-c] Gamma[c] Gamma[-a1-b1+c] Pochhammer[1-a1,-m+n] Pochhammer[a1,m] Pochhammer[1-b1,-m+n] Pochhammer[b1,m] Pochhammer[-a1-b1+c,-m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1-a1-a2-b1+c,-m+n] Pochhammer[1-a1-b1-b2+c,-m+n])},{Abs[1-y]<1&&1/Abs[x]<1&&1+1/Abs[x]<1/Abs[x-x y]&&Abs[x (-1+y)]<1&&Abs[1-y]+Abs[x (-1+y)]<1&&Abs[x-x y]<1,((-x)^-a1 x^-m (1-y)^n Gamma[-a1+b1] Gamma[c] Gamma[-a1-a2-b2+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[1+a1+a2-c,m] Pochhammer[1+a1+b2-c,m])/(m! n! Gamma[b1] Gamma[-a1-a2+c] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2+b2-c,m+n])+((-x)^-b1 x^-m (1-y)^n Gamma[a1-b1] Gamma[c] Gamma[-a2-b1-b2+c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a2+b1-c,m] Pochhammer[1+b1+b2-c,m])/(m! n! Gamma[a1] Gamma[-a2-b1+c] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1+b2-c,m+n])+((-1)^m x^m (1-y)^(-a2-b2+c+m+n) Gamma[a2+b2-c] Gamma[c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[-a2+c,m+n] Pochhammer[-b2+c,m+n])/(m! n! Gamma[a2] Gamma[b2] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m] Pochhammer[1-a2-b2+c,m+n])+((-1)^n (-x)^(a2+b2-c+n) x^-m (1-y)^n Gamma[a1+a2+b2-c] Gamma[a2+b1+b2-c] Gamma[c] Gamma[-a2-b2+c] Pochhammer[1-a2,m-n] Pochhammer[a2,n] Pochhammer[1-b2,m-n] Pochhammer[b2,n] Pochhammer[-a2-b2+c,m-n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1-a1-a2-b2+c,m-n] Pochhammer[1-a2-b1-b2+c,m-n])},{Abs[1-x]>1&&Abs[y]<1&&1/Abs[1-x]+Abs[y]<1,((1-x)^(-b1-m) y^n Gamma[a1-b1] Gamma[c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a1+c,m+n])/(m! n! Gamma[a1] Gamma[-b1+c] Pochhammer[1-a1+b1,m] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n])+((1-x)^(-a1-m) y^n Gamma[-a1+b1] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[-b1+c,m+n])/(m! n! Gamma[b1] Gamma[-a1+c] Pochhammer[1+a1-b1,m] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n])},{Abs[1-y]>1&&Abs[x]<1&&Abs[x]+1/Abs[1-y]<1,(x^m (1-y)^(-b2-n) Gamma[a2-b2] Gamma[c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a2+c,m+n])/(m! n! Gamma[a2] Gamma[-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m])+(x^m (1-y)^(-a2-n) Gamma[-a2+b2] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[-b2+c,m+n])/(m! n! Gamma[b2] Gamma[-a2+c] Pochhammer[1+a2-b2,n] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m])},{1/Abs[-1+x]<1&&1/Abs[y]<1&&1/Abs[y]<1/(1+1/Abs[-1+x]),((-1)^m (1-x)^(-b1-m) (-y)^-a2 y^-n Gamma[a1-b1] Gamma[-a2+b2] Gamma[c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a1+a2-c,n] Pochhammer[1+a2+b1-c,n])/(m! n! Gamma[a1] Gamma[b2] Gamma[-a2-b1+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2-c,-m+n])+((-1)^m (1-x)^(-a1-m) (-y)^-a2 y^-n Gamma[-a1+b1] Gamma[-a2+b2] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,n] Pochhammer[1+a2+b1-c,n])/(m! n! Gamma[b1] Gamma[b2] Gamma[-a1-a2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a2-b2,n] Pochhammer[1+a2+b1-c,-m+n])+((-1)^m (1-x)^(-b1-m) (-y)^-b2 y^-n Gamma[a1-b1] Gamma[a2-b2] Gamma[c] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,n] Pochhammer[1+b1+b2-c,n])/(m! n! Gamma[a1] Gamma[a2] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b2-c,-m+n])+((-1)^m (1-x)^(-a1-m) (-y)^-b2 y^-n Gamma[-a1+b1] Gamma[a2-b2] Gamma[c] Pochhammer[a1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,n] Pochhammer[1+b1+b2-c,n])/(m! n! Gamma[a2] Gamma[b1] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1-a2+b2,n] Pochhammer[1+b1+b2-c,-m+n])},{1/Abs[-1+y]<1&&1/Abs[x]<1&&1/Abs[x]<1/(1+1/Abs[-1+y]),((-1)^n (-x)^-a1 x^-m (1-y)^(-b2-n) Gamma[-a1+b1] Gamma[a2-b2] Gamma[c] Pochhammer[a1,m] Pochhammer[b2,n] Pochhammer[1+a1+a2-c,m] Pochhammer[1+a1+b2-c,m])/(m! n! Gamma[a2] Gamma[b1] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+a2-c,m-n])+((-1)^n (-x)^-a1 x^-m (1-y)^(-a2-n) Gamma[-a1+b1] Gamma[-a2+b2] Gamma[c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,m] Pochhammer[1+a1+b2-c,m])/(m! n! Gamma[b1] Gamma[b2] Gamma[-a1-a2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+b2-c,m-n])+((-1)^n (-x)^-b1 x^-m (1-y)^(-b2-n) Gamma[a1-b1] Gamma[a2-b2] Gamma[c] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a2+b1-c,m] Pochhammer[1+b1+b2-c,m])/(m! n! Gamma[a1] Gamma[a2] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1-a2+b2,n] Pochhammer[1+a2+b1-c,m-n])+((-1)^n (-x)^-b1 x^-m (1-y)^(-a2-n) Gamma[a1-b1] Gamma[-a2+b2] Gamma[c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a2+b1-c,m] Pochhammer[1+b1+b2-c,m])/(m! n! Gamma[a1] Gamma[b2] Gamma[-a2-b1+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2-b2,n] Pochhammer[1+b1+b2-c,m-n])},{Abs[1-x]<1&&Abs[-1+x]<1&&1+Abs[-1+x]<Abs[y-x y]&&1/Abs[y-x y]<1&&1/Abs[y-x y]<1/(1+Abs[-1+x])&&1/Abs[y]<1,((-1)^n (1/((-1+x) y))^m (-y)^(a1+b1-c) y^-n (y-x y)^n Gamma[-a2+b2] Gamma[1+a1+b1-c] Gamma[a1+a2+b1-c] Gamma[c] Gamma[-a1-b1+c] If[-1+Re[x]+Re[y]>0,(y-x y)^(-a1-a2-b1+c),(1/(y-x y))^(a1+a2+b1-c)] Pochhammer[a2,m] Pochhammer[1+a1+a2-c,m] Pochhammer[1+a2+b1-c,m] Pochhammer[a1+a2+b1-c,m-n])/(m! n! Gamma[a1] Gamma[1-a2] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1+a2-b2,m] Pochhammer[1+a1+a2-c,m-n] Pochhammer[1+a2+b1-c,m-n])+((-1)^m (1-x)^(-a1-b1+c+m) (1/((-1+x) y))^n Gamma[-a2+b2] Gamma[a1+b1-c] Gamma[c] Gamma[1-a1-b1+c] If[-1+Re[x]+Re[y]>0,(y-x y)^-a2,(1/(y-x y))^a2] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,n] Pochhammer[1+a2+b1-c,n] Pochhammer[a1+a2+b1-c,-m+n])/(m! n! Gamma[a1] Gamma[b1] Gamma[b2] Gamma[1-a1-a2-b1+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2-c,-m+n] Pochhammer[1+a2+b1-c,-m+n])+((1-x)^m (-y)^-a2 y^-n Gamma[-a2+b2] Gamma[c] Gamma[-a1-a2-b1+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a1+a2-c,n] Pochhammer[1+a2+b1-c,n])/(m! n! Gamma[b2] Gamma[-a1-a2+c] Gamma[-a2-b1+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2+b1-c,m+n])+((-1)^n (1/((-1+x) y))^m (-y)^(a1+b1-c) y^-n (y-x y)^n Gamma[a2-b2] Gamma[1+a1+b1-c] Gamma[a1+b1+b2-c] Gamma[c] Gamma[-a1-b1+c] If[-1+Re[x]+Re[y]>0,(y-x y)^(-a1-b1-b2+c),(1/(y-x y))^(a1+b1+b2-c)] Pochhammer[b2,m] Pochhammer[1+a1+b2-c,m] Pochhammer[1+b1+b2-c,m] Pochhammer[a1+b1+b2-c,m-n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[1-b2] Gamma[b2] Pochhammer[1-a2+b2,m] Pochhammer[1+a1+b2-c,m-n] Pochhammer[1+b1+b2-c,m-n])+((-1)^m (1-x)^(-a1-b1+c+m) (1/((-1+x) y))^n Gamma[a2-b2] Gamma[a1+b1-c] Gamma[c] Gamma[1-a1-b1+c] If[-1+Re[x]+Re[y]>0,(y-x y)^-b2,(1/(y-x y))^b2] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,n] Pochhammer[1+b1+b2-c,n] Pochhammer[a1+b1+b2-c,-m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[1-a1-b1-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b2-c,-m+n] Pochhammer[1+b1+b2-c,-m+n])+((1-x)^m (-y)^-b2 y^-n Gamma[a2-b2] Gamma[c] Gamma[-a1-b1-b2+c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,n] Pochhammer[1+b1+b2-c,n])/(m! n! Gamma[a2] Gamma[-a1-b2+c] Gamma[-b1-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b1+b2-c,m+n])},{Abs[1-y]<1&&Abs[-1+y]<1&&1+Abs[-1+y]<Abs[x-x y]&&1/Abs[x-x y]<1&&1/Abs[x-x y]<1/(1+Abs[-1+y])&&1/Abs[x]<1,((-1)^n (1-y)^(-a2-b2+c+n) (1/(x (-1+y)))^m Gamma[-a1+b1] Gamma[a2+b2-c] Gamma[c] Gamma[1-a2-b2+c] If[-1+Re[x]+Re[y]>0,(x-x y)^-a1,(1/(x-x y))^a1] Pochhammer[a1,m] Pochhammer[1+a1+a2-c,m] Pochhammer[1+a1+b2-c,m] Pochhammer[a1+a2+b2-c,m-n])/(m! n! Gamma[a2] Gamma[b1] Gamma[b2] Gamma[1-a1-a2-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2-c,m-n] Pochhammer[1+a1+b2-c,m-n])+((-1)^m (-x)^(a2+b2-c) x^-m (1/(x (-1+y)))^n (x-x y)^m Gamma[-a1+b1] Gamma[1+a2+b2-c] Gamma[a1+a2+b2-c] Gamma[c] Gamma[-a2-b2+c] If[-1+Re[x]+Re[y]>0,(x-x y)^(-a1-a2-b2+c),(1/(x-x y))^(a1+a2+b2-c)] Pochhammer[a1,n] Pochhammer[1+a1+a2-c,n] Pochhammer[1+a1+b2-c,n] Pochhammer[a1+a2+b2-c,-m+n])/(m! n! Gamma[1-a1] Gamma[a1] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1+a1-b1,n] Pochhammer[1+a1+a2-c,-m+n] Pochhammer[1+a1+b2-c,-m+n])+((-x)^-a1 x^-m (1-y)^n Gamma[-a1+b1] Gamma[c] Gamma[-a1-a2-b2+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[1+a1+a2-c,m] Pochhammer[1+a1+b2-c,m])/(m! n! Gamma[b1] Gamma[-a1-a2+c] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2+b2-c,m+n])+((-1)^n (1-y)^(-a2-b2+c+n) (1/(x (-1+y)))^m Gamma[a1-b1] Gamma[a2+b2-c] Gamma[c] Gamma[1-a2-b2+c] If[-1+Re[x]+Re[y]>0,(x-x y)^-b1,(1/(x-x y))^b1] Pochhammer[b1,m] Pochhammer[1+a2+b1-c,m] Pochhammer[1+b1+b2-c,m] Pochhammer[a2+b1+b2-c,m-n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b2] Gamma[1-a2-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1-c,m-n] Pochhammer[1+b1+b2-c,m-n])+((-1)^m (-x)^(a2+b2-c) x^-m (1/(x (-1+y)))^n (x-x y)^m Gamma[a1-b1] Gamma[1+a2+b2-c] Gamma[a2+b1+b2-c] Gamma[c] Gamma[-a2-b2+c] If[-1+Re[x]+Re[y]>0,(x-x y)^(-a2-b1-b2+c),(1/(x-x y))^(a2+b1+b2-c)] Pochhammer[b1,n] Pochhammer[1+a2+b1-c,n] Pochhammer[1+b1+b2-c,n] Pochhammer[a2+b1+b2-c,-m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[1-b1] Gamma[b1] Gamma[b2] Pochhammer[1-a1+b1,n] Pochhammer[1+a2+b1-c,-m+n] Pochhammer[1+b1+b2-c,-m+n])+((-x)^-b1 x^-m (1-y)^n Gamma[a1-b1] Gamma[c] Gamma[-a2-b1-b2+c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a2+b1-c,m] Pochhammer[1+b1+b2-c,m])/(m! n! Gamma[a1] Gamma[-a2-b1+c] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1+b2-c,m+n])},{Abs[x/(-1+x)]<1&&Abs[y]<1,((1-x)^-b1 (x/(-1+x))^m y^n Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a1+c,m+n])/(m! n! Pochhammer[c,m+n] Pochhammer[-a1+c,n])},{Abs[y/(-1+y)]<1&&Abs[x]<1,(x^m (1-y)^-b2 (y/(-1+y))^n Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a2+c,m+n])/(m! n! Pochhammer[c,m+n] Pochhammer[-a2+c,m])},{Abs[x/(-1+x)]<1&&1/Abs[y]<1&&Abs[x/(-1+x)]<1&&1/Abs[y]<1,((1-x)^-b1 (x/(-1+x))^m (-y)^-a2 y^-n Gamma[-a2+b2] Gamma[c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a2-c,-m+n] Pochhammer[1+a1+a2-c,n])/(m! n! Gamma[b2] Gamma[-a2+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2-c,-m+n])+((1-x)^-b1 (x/(-1+x))^m (-y)^-b2 y^-n Gamma[a2-b2] Gamma[c] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b2-c,-m+n] Pochhammer[1+a1+b2-c,n])/(m! n! Gamma[a2] Gamma[-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b2-c,-m+n])},{Abs[y/(-1+y)]<1&&1/Abs[x]<1&&Abs[y/(-1+y)]<1&&1/Abs[x]<1,((-x)^-a1 x^-m (1-y)^-b2 (y/(-1+y))^n Gamma[-a1+b1] Gamma[c] Pochhammer[a1,m] Pochhammer[b2,n] Pochhammer[1+a1-c,m-n] Pochhammer[1+a1+a2-c,m])/(m! n! Gamma[b1] Gamma[-a1+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2-c,m-n])+((-x)^-b1 x^-m (1-y)^-b2 (y/(-1+y))^n Gamma[a1-b1] Gamma[c] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c,m-n] Pochhammer[1+a2+b1-c,m])/(m! n! Gamma[a1] Gamma[-b1+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1-c,m-n])},{Abs[(-1+x)/x]<1&&Abs[y]>1&&Abs[((-1+x) y)/x]<1,((1/x)^a1 ((-1+x)/x)^m (-y)^-a2 y^-n Gamma[-a2+b2] Gamma[c] Gamma[-a1-a2-b1+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,m+n] Pochhammer[1+a2+b1-c,n])/(m! n! Gamma[b2] Gamma[-a1-a2+c] Gamma[-a2-b1+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2+b1-c,m+n])+((1/x)^a1 ((-1+x)/x)^m (-y)^-b2 y^-n Gamma[a2-b2] Gamma[c] Gamma[-a1-b1-b2+c] Pochhammer[a1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,m+n] Pochhammer[1+b1+b2-c,n])/(m! n! Gamma[a2] Gamma[-a1-b2+c] Gamma[-b1-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b1+b2-c,m+n])+((-1)^n (1-x)^(-a1-b1+c+n) (1/x)^(-a1+c) ((-1+x)/x)^m x^-n y^n Gamma[a1+b1-c] Gamma[c] Pochhammer[1-a1,m] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[-a1+c,m+n])/(m! n! Gamma[a1] Gamma[b1] Pochhammer[-a1+c,n] Pochhammer[1-a1-b1+c,m+n])+((-1)^m (1/x)^a1 ((-1+x)/x)^m (-y)^(a1+b1-c+m) y^-n Gamma[a1+a2+b1-c] Gamma[a1+b1+b2-c] Gamma[c] Gamma[-a1-b1+c] Pochhammer[1-a1,-m+n] Pochhammer[a1,m] Pochhammer[1-b1,n] Pochhammer[-a1-b1+c,-m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1-a1-a2-b1+c,-m+n] Pochhammer[1-a1-b1-b2+c,-m+n])},{Abs[(-1+y)/y]<1&&Abs[x]>1&&Abs[(x (-1+y))/y]<1,((-x)^-a1 x^-m (1/y)^a2 ((-1+y)/y)^n Gamma[-a1+b1] Gamma[c] Gamma[-a1-a2-b2+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,m+n] Pochhammer[1+a1+b2-c,m])/(m! n! Gamma[b1] Gamma[-a1-a2+c] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2+b2-c,m+n])+((-x)^-b1 x^-m (1/y)^a2 ((-1+y)/y)^n Gamma[a1-b1] Gamma[c] Gamma[-a2-b1-b2+c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a2+b1-c,m+n] Pochhammer[1+b1+b2-c,m])/(m! n! Gamma[a1] Gamma[-a2-b1+c] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1+b2-c,m+n])+((-1)^m x^m (1-y)^(-a2-b2+c+m) (1/y)^(-a2+c) ((-1+y)/y)^n y^-m Gamma[a2+b2-c] Gamma[c] Pochhammer[a1,m] Pochhammer[1-a2,n] Pochhammer[b1,m] Pochhammer[-a2+c,m+n])/(m! n! Gamma[a2] Gamma[b2] Pochhammer[-a2+c,m] Pochhammer[1-a2-b2+c,m+n])+((-1)^n (-x)^(a2+b2-c+n) x^-m (1/y)^a2 ((-1+y)/y)^n Gamma[a1+a2+b2-c] Gamma[a2+b1+b2-c] Gamma[c] Gamma[-a2-b2+c] Pochhammer[1-a2,m-n] Pochhammer[a2,n] Pochhammer[1-b2,m] Pochhammer[-a2-b2+c,m-n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1-a1-a2-b2+c,m-n] Pochhammer[1-a2-b1-b2+c,m-n])},{Abs[x]>1&&Abs[((-1+x) y)/x]<1&&1/Abs[x]+Abs[((-1+x) y)/x]<1,(((-1+x)/x)^n (-x)^-b1 x^-m (x/(-1+x))^(a1+b1-c) y^n Gamma[a1-b1] Gamma[c] Pochhammer[1-a1,m] Pochhammer[a2,n] Pochhammer[b2,n] Pochhammer[-a1+c,m+n])/(m! n! Gamma[a1] Gamma[-b1+c] Pochhammer[1-a1+b1,m] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n])+(((-1+x)/x)^n (-x)^-a1 x^-m (x/(-1+x))^(a1+b1-c) y^n Gamma[-a1+b1] Gamma[c] Pochhammer[a2,n] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[-b1+c,m+n])/(m! n! Gamma[b1] Gamma[-a1+c] Pochhammer[1+a1-b1,m] Pochhammer[-a1+c,n] Pochhammer[-b1+c,n])},{Abs[y]>1&&Abs[(x (-1+y))/y]<1&&Abs[(x (-1+y))/y]+1/Abs[y]<1,(x^m ((-1+y)/y)^m (-y)^-b2 y^-n (y/(-1+y))^(a2+b2-c) Gamma[a2-b2] Gamma[c] Pochhammer[a1,m] Pochhammer[1-a2,n] Pochhammer[b1,m] Pochhammer[-a2+c,m+n])/(m! n! Gamma[a2] Gamma[-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m])+(x^m ((-1+y)/y)^m (-y)^-a2 y^-n (y/(-1+y))^(a2+b2-c) Gamma[-a2+b2] Gamma[c] Pochhammer[a1,m] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[-b2+c,m+n])/(m! n! Gamma[b2] Gamma[-a2+c] Pochhammer[1+a2-b2,n] Pochhammer[-a2+c,m] Pochhammer[-b2+c,m])},{Abs[(-1+x)/x]<1&&Abs[y]>1&&Abs[x/(y-x y)]<1,((-1)^n (1/x)^a1 (x/((-1+x) y))^m (-y)^(a1+b1-c) y^-n ((-1+1/x) y)^n Gamma[-a2+b2] Gamma[1+a1+b1-c] Gamma[a1+a2+b1-c] Gamma[c] Gamma[-a1-b1+c] If[-Re[x]+Re[x]^2+Re[y]>0,((-1+1/x) y)^(-a1-a2-b1+c),(1/((-1+1/x) y))^(a1+a2+b1-c)] Pochhammer[a2,m] Pochhammer[1-b1,n] Pochhammer[1+a2+b1-c,m] Pochhammer[a1+a2+b1-c,m-n])/(m! n! Gamma[a1] Gamma[1-a2] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1+a2-b2,m] Pochhammer[1+a2+b1-c,m-n])+((1-x)^(-a1-b1+c) (1/x)^(-a1+c) ((-1+x)/x)^m (x/((-1+x) y))^n Gamma[-a2+b2] Gamma[a1+b1-c] Gamma[c] Gamma[1-a1-b1+c] If[-Re[x]+Re[x]^2+Re[y]>0,((-1+1/x) y)^-a2,(1/((-1+1/x) y))^a2] Pochhammer[1-a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,n] Pochhammer[a1+a2+b1-c,-m+n])/(m! n! Gamma[a1] Gamma[b1] Gamma[b2] Gamma[1-a1-a2-b1+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2-c,-m+n])+((1/x)^a1 ((-1+x)/x)^m (-y)^-a2 y^-n Gamma[-a2+b2] Gamma[c] Gamma[-a1-a2-b1+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,m+n] Pochhammer[1+a2+b1-c,n])/(m! n! Gamma[b2] Gamma[-a1-a2+c] Gamma[-a2-b1+c] Pochhammer[1+a2-b2,n] Pochhammer[1+a1+a2+b1-c,m+n])+((-1)^n (1/x)^a1 (x/((-1+x) y))^m (-y)^(a1+b1-c) y^-n ((-1+1/x) y)^n Gamma[a2-b2] Gamma[1+a1+b1-c] Gamma[a1+b1+b2-c] Gamma[c] Gamma[-a1-b1+c] If[-Re[x]+Re[x]^2+Re[y]>0,((-1+1/x) y)^(-a1-b1-b2+c),(1/((-1+1/x) y))^(a1+b1+b2-c)] Pochhammer[1-b1,n] Pochhammer[b2,m] Pochhammer[1+b1+b2-c,m] Pochhammer[a1+b1+b2-c,m-n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[1-b2] Gamma[b2] Pochhammer[1-a2+b2,m] Pochhammer[1+b1+b2-c,m-n])+((1-x)^(-a1-b1+c) (1/x)^(-a1+c) ((-1+x)/x)^m (x/((-1+x) y))^n Gamma[a2-b2] Gamma[a1+b1-c] Gamma[c] Gamma[1-a1-b1+c] If[-Re[x]+Re[x]^2+Re[y]>0,((-1+1/x) y)^-b2,(1/((-1+1/x) y))^b2] Pochhammer[1-a1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,n] Pochhammer[a1+b1+b2-c,-m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b1] Gamma[1-a1-b1-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b2-c,-m+n])+((1/x)^a1 ((-1+x)/x)^m (-y)^-b2 y^-n Gamma[a2-b2] Gamma[c] Gamma[-a1-b1-b2+c] Pochhammer[a1,m] Pochhammer[b2,n] Pochhammer[1+a1+b2-c,m+n] Pochhammer[1+b1+b2-c,n])/(m! n! Gamma[a2] Gamma[-a1-b2+c] Gamma[-b1-b2+c] Pochhammer[1-a2+b2,n] Pochhammer[1+a1+b1+b2-c,m+n])},{Abs[(-1+y)/y]<1&&Abs[x]>1&&Abs[y/(x-x y)]<1,((1-y)^(-a2-b2+c) (1/y)^(-a2+c) ((-1+y)/y)^n (y/(x (-1+y)))^m Gamma[-a1+b1] Gamma[a2+b2-c] Gamma[c] Gamma[1-a2-b2+c] If[Re[x]-Re[y]+Re[y]^2>0,(x (-1+1/y))^-a1,(1/(x (-1+1/y)))^a1] Pochhammer[a1,m] Pochhammer[1-a2,n] Pochhammer[1+a1+a2-c,m] Pochhammer[a1+a2+b2-c,m-n])/(m! n! Gamma[a2] Gamma[b1] Gamma[b2] Gamma[1-a1-a2-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2-c,m-n])+((-1)^m (-x)^(a2+b2-c) x^-m (x (-1+1/y))^m (1/y)^a2 (y/(x (-1+y)))^n Gamma[-a1+b1] Gamma[1+a2+b2-c] Gamma[a1+a2+b2-c] Gamma[c] Gamma[-a2-b2+c] If[Re[x]-Re[y]+Re[y]^2>0,(x (-1+1/y))^(-a1-a2-b2+c),(1/(x (-1+1/y)))^(a1+a2+b2-c)] Pochhammer[a1,n] Pochhammer[1-b2,m] Pochhammer[1+a1+b2-c,n] Pochhammer[a1+a2+b2-c,-m+n])/(m! n! Gamma[1-a1] Gamma[a1] Gamma[a2] Gamma[b1] Gamma[b2] Pochhammer[1+a1-b1,n] Pochhammer[1+a1+b2-c,-m+n])+((-x)^-a1 x^-m (1/y)^a2 ((-1+y)/y)^n Gamma[-a1+b1] Gamma[c] Gamma[-a1-a2-b2+c] Pochhammer[a1,m] Pochhammer[a2,n] Pochhammer[1+a1+a2-c,m+n] Pochhammer[1+a1+b2-c,m])/(m! n! Gamma[b1] Gamma[-a1-a2+c] Gamma[-a1-b2+c] Pochhammer[1+a1-b1,m] Pochhammer[1+a1+a2+b2-c,m+n])+((1-y)^(-a2-b2+c) (1/y)^(-a2+c) ((-1+y)/y)^n (y/(x (-1+y)))^m Gamma[a1-b1] Gamma[a2+b2-c] Gamma[c] Gamma[1-a2-b2+c] If[Re[x]-Re[y]+Re[y]^2>0,(x (-1+1/y))^-b1,(1/(x (-1+1/y)))^b1] Pochhammer[1-a2,n] Pochhammer[b1,m] Pochhammer[1+a2+b1-c,m] Pochhammer[a2+b1+b2-c,m-n])/(m! n! Gamma[a1] Gamma[a2] Gamma[b2] Gamma[1-a2-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1-c,m-n])+((-1)^m (-x)^(a2+b2-c) x^-m (x (-1+1/y))^m (1/y)^a2 (y/(x (-1+y)))^n Gamma[a1-b1] Gamma[1+a2+b2-c] Gamma[a2+b1+b2-c] Gamma[c] Gamma[-a2-b2+c] If[Re[x]-Re[y]+Re[y]^2>0,(x (-1+1/y))^(-a2-b1-b2+c),(1/(x (-1+1/y)))^(a2+b1+b2-c)] Pochhammer[b1,n] Pochhammer[1-b2,m] Pochhammer[1+b1+b2-c,n] Pochhammer[a2+b1+b2-c,-m+n])/(m! n! Gamma[a1] Gamma[a2] Gamma[1-b1] Gamma[b1] Gamma[b2] Pochhammer[1-a1+b1,n] Pochhammer[1+b1+b2-c,-m+n])+((-x)^-b1 x^-m (1/y)^a2 ((-1+y)/y)^n Gamma[a1-b1] Gamma[c] Gamma[-a2-b1-b2+c] Pochhammer[a2,n] Pochhammer[b1,m] Pochhammer[1+a2+b1-c,m+n] Pochhammer[1+b1+b2-c,m])/(m! n! Gamma[a1] Gamma[-a2-b1+c] Gamma[-b1-b2+c] Pochhammer[1-a1+b1,m] Pochhammer[1+a2+b1+b2-c,m+n])}};


End[]


EndPackage[]
