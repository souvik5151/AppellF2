(* ::Package:: *)

BeginPackage["AppellF2`"]


Print["AppellF2.wl v1.0\n","Authors : Souvik Bera & Tanay Pathak"];


AppellF2::usage="The command gives the numerical value of the Appell Function F2.
 F2[a,b1,b2,c1,c2,x,y,precision,terms, F2show-> True]";
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
If[list>44,Print["Total number of ACs is 44"];Abort[];];
ClearAll[Global`x,Global`y,Global`a,Global`b1,Global`b2,Global`c1,Global`c2,Global`m,Global`n];

Return[{Subscript[Roc, list][Global`x,Global`y],Subscript[AC, list][Global`a,Global`b1,Global`b2,Global`c1,Global`c2,Global`x,Global`y,Global`m,Global`n]/.Gamma[x_. +1]-> x Gamma[x]/.Re[x_]-> x //Simplify//Expand}]];



F2findall[{x_?NumberQ, y_?NumberQ}]:=Module[{test,pos},

test=Table[Subscript[Roc, i][Rationalize[x],Rationalize[y]],{i,1,44}];
pos=Position[test,True];
Return[Flatten[pos]];


]


F2evaluate[list_/;(IntegerQ[list]&&list>0), para_List, p_,terms_]:=Module[{a,b1,b2,c1,c2,x,y,result,p0},
p0=p+10;

If[Length[para]=!=7,Abort[];];
{a,b1,b2,c1,c2,x,y}=SetPrecision[para,p0];



If[Subscript[Roc, list][Rationalize[x],Rationalize[y]],{},Print["The point does not lie in ", list];Abort[];];


Block[{sum,m1,n1},{sum=0;Monitor[For[m1=0,m1<=terms,m1++,{For[n1=0,n1<=terms,n1++,
{sum=sum+SetPrecision[Subscript[AC, list][a,b1,b2,c1,c2,x,y,m1,n1],p0];}]};(*Pause[10^-6];*)],{SetPrecision[sum,p0],m1,n1}];result=sum;}];


Return[SetPrecision[Chop[result],p]]]


F2ROC[point___List,list_/;(IntegerQ[list]&&list>0),range_List]:=Module[{i},
Show[ListPlot[{point},PlotStyle->Directive[PointSize[Medium],Red],PlotRange->{{range[[1]],range[[2]]},{range[[1]],range[[2]]}}],RegionPlot[Subscript[Roc,list][x,y],{x,range[[1]]-1,range[[2]]+1},{y,range[[1]]-1,range[[2]]+1},PlotPoints->100],PlotRange->{{range[[1]],range[[2]]},{range[[1]],range[[2]]}},AspectRatio-> 1]
]


Options[AppellF2]={F2show->False};
AppellF2::poch="Pochhammer parameters do not obey the constraint";
AppellF2::val="The point is not covered";
AppellF2[args___] := f2core[AppellF2,args];

f2core[AppellF2,a0_,b0_,c0_,d0_,e0_,x0_,y0_,p_,terms_,OptionsPattern[AppellF2]]:=Module[{aa,bb1,bb2,cc1,cc2,xx,yy,a,b1,b2,c1,c2,x1,y1,x,y,m,n,m0,m1,n0,n1,v,i,test1,test2,pos1,pos2,result,p0,integrand2,integrand3,seriesselect},
p0=p+10;m1=100;n1=100;
(*Print[{a0,b0,c0,d0,e0,x0,y0}];*)
(*Label[start];
Print["start"];*)

a=SetPrecision[a0,p0]; b1=SetPrecision[b0,p0]; b2=SetPrecision[c0,p0]; c1=SetPrecision[d0,p0]; c2=SetPrecision[e0,p0];x=SetPrecision[x0,p0]; y=SetPrecision[y0,p0];
 result=0;



(* If x=0 or y=0*)

If[x0===0.0||x0===0,result=Hypergeometric2F1[a,b2,c2,y];Goto[end];];
If[y0===0.0||y0===0,result=Hypergeometric2F1[a,b1,c1,x];Goto[end];];






Off[General::infy,General::indet];
(* check the region *)
test1=Table[Subscript[Roc, i][Rationalize[x0],Rationalize[y0]],{i,1,44}];
(*Print[test1];*)

PrintTemporary["Finding a proper AC..."];

pos1=Position[test1,True];

If[pos1==={},Message[AppellF2::val];Abort[];];

Print["valid series : ",pos1];





test2=Table[{Max[Norm[#]&/@(Abs[{Limit[Simplify[(#1/.{m-> 0,n->n+ 1})/(#1/.m-> 0)],{n-> n0},PerformanceGoal-> "Speed"],Limit[Simplify[(#1/.{n-> 0,m->m+ 1})/(#1/.n-> 0)],{m-> m0},PerformanceGoal-> "Speed"]}/.{m0-> m1,n0-> n1}/.{aa-> a,bb1-> b1,bb2-> b2,cc1-> c1,cc2-> c2,xx-> x,yy-> y}]&/@If[Head[#]===Plus,List@@#,{#}]&@Expand[Subscript[AC, pos1[[i]][[1]]][aa,bb1,bb2,cc1,cc2,xx,yy,m,n]])],pos1[[i]][[1]]},{i,Length[pos1]}];



(*Print[test2];*)

seriesselect=SortBy[test2,First];
Print["convergence rates :",seriesselect];

For[i=1,i<= Length[test2],i++,
pos2=seriesselect[[i]][[2]];
If[MemberQ[{"0","Indeterminate"},ToString[Sum[Subscript[AC, pos2][a,b1,b2,c1,c2,x,y,m1,n1],{m1,0,1},{n1,0,1}]]],Continue[],Break[]];];


Print["selected series : ",pos2];

PrintTemporary["Evaluating sum..."];


If[OptionValue[F2show]===True, Block[{sum,m11,n11},{sum=0;Monitor[For[m11=0,m11<=terms,m11++,{For[n11=0,n11<=terms,n11++,
{sum=sum+SetPrecision[Subscript[AC, pos2][a,b1,b2,c1,c2,x,y,m11,n11],p0];}]};],{SetPrecision[sum,p0],m11,n11}];result=sum;}];,
Block[{sum,m11,n11},{sum=0;For[m11=0,m11<=terms,m11++,{For[n11=0,n11<=terms,n11++,
{sum=sum+SetPrecision[Subscript[AC, pos2][a,b1,b2,c1,c2,x,y,m11,n11],p0];}]}];result=sum;}];


];



Label[end];
Return[SetPrecision[Chop[result],p]]
]



(* ::Subchapter:: *)
(*Set of 44 ACs*)


(* ::Input::Initialization:: *)
Subscript[AC, 1][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= (Pochhammer[a,m+n]Pochhammer[b1,m]Pochhammer[b2,n])/(Pochhammer[c1,m]Pochhammer[c2,n]m! n!) x^m y^n;


Subscript[Roc, 1][x_,y_]:= Abs[x]+ Abs[y]<1 && Abs[x]<1;


Subscript[AC, 2][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=(x^m (1-y)^n Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a-c2,m])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((x/(1-y))^m (1-y)^(-a-b2+c2+n) Gamma[a+b2-c2] Gamma[c2] Pochhammer[b1,m] Pochhammer[1+a-c2,m] Pochhammer[-a+c2,-m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Pochhammer[c1,m] Pochhammer[1-a-b2+c2,-m+n]);


Subscript[Roc, 2][x_,y_]:= Abs[x/(1-y)]<1&&Abs[1-y]<1;


Subscript[AC, 3][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x)^m y^n Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a-c1,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n])+((1-x)^(-a-b1+c1+m) (y/(1-x))^n Gamma[a+b1-c1] Gamma[c1] Pochhammer[b2,n] Pochhammer[1+a-c1,n] Pochhammer[-a+c1,m-n] Pochhammer[-b1+c1,m])/(m! n! Gamma[a] Gamma[b1] Pochhammer[1-a-b1+c1,m-n] Pochhammer[c2,n]);


Subscript[Roc, 3][x_,y_]:=Abs[y/(1-x)]<1&&Abs[1-x]<1;


Subscript[AC, 4][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= (x^m (1-y)^n Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a-c2,m])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+(x^n (1-y)^(-a-b2+c2) ((1-y)/x)^m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a+b1-b2+c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^(-a-b2+c2),((-1+y)/x)^(a+b2-c2)] Pochhammer[b2,m] Pochhammer[a+b2-c2,m-n] Pochhammer[1+a+b2-c1-c2,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b2+c1+c2] Pochhammer[b2,m-n] Pochhammer[1+a-b1+b2-c2,m-n])+((1-y)^(-a-b2+c2+n) ((1-y)/x)^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1+c2,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[-a+b1+c2,m] Pochhammer[1-a+b1-b2+c2,m+n]);


Subscript[Roc, 4][x_,y_]:=Abs[(1-y)/x]<1&&Abs[x]<1;


Subscript[AC, 5][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^(-a-b1+c1+m) ((1-x)/y)^n Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[-a+b2+c1,m+n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[-a+b2+c1,n] Pochhammer[1-a-b1+b2+c1,m+n])+((1-x)^(-a-b1+c1) ((1-x)/y)^n y^m Gamma[a+b1-c1] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^(-a-b1+c1),((-1+x)/y)^(a+b1-c1)] Pochhammer[b1,n] Pochhammer[a+b1-c1,-m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a+b1-c1-c2,-m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b1+c1+c2] Pochhammer[b1,-m+n] Pochhammer[1+a+b1-b2-c1,-m+n])+((1-x)^m y^n Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a-c1,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n]);


Subscript[Roc, 5][x_,y_]:= Abs[(1-x)/y]<1&&Abs[y]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 6][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1/x)^m (-x)^-a ((1-y)/x)^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,m])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,m+n])+((1/x)^m (-x)^-b1 (1-y)^n Gamma[a-b1] Gamma[c1] Gamma[1+a-b1-c2] Gamma[a-b1+b2-c2] Gamma[c2] Gamma[1-a+b1-b2+c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1-b2+c2,m-n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[a-b1-c2] Gamma[1+a-b1+b2-c2] Gamma[1-a+b1+c2] Gamma[-b2+c2] Pochhammer[1-a+b1,m-n] Pochhammer[-a+b1+c2,m])+((-x)^-b1 (1-y)^(-a+b1-b2+c2+n) ((1-y)/x)^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1+c2,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[-a+b1+c2,m] Pochhammer[1-a+b1-b2+c2,m+n]);


(* ::Input::Initialization:: *)
Subscript[Roc, 6][x_,y_]:=1/Abs[x]<1&&Abs[1-y]<1;


Subscript[AC, 7][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1/y)^n ((1-x)/y)^m (-y)^-a Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[-b1+c1,n] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m+n])+((1-x)^m (1/y)^n (-y)^-b2 Gamma[a-b2] Gamma[1+a-b2-c1] Gamma[a+b1-b2-c1] Gamma[c1] Gamma[1-a-b1+b2+c1] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a-b1+b2+c1,-m+n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[a-b2-c1] Gamma[1+a+b1-b2-c1] Gamma[-b1+c1] Gamma[1-a+b2+c1] Gamma[-b2+c2] Pochhammer[1-a+b2,-m+n] Pochhammer[-a+b2+c1,n])+((1-x)^(-a-b1+b2+c1+m) ((1-x)/y)^n (-y)^-b2 Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[-a+b2+c1,m+n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[-a+b2+c1,n] Pochhammer[1-a-b1+b2+c1,m+n]);


Subscript[Roc, 7][x_,y_]:=1/Abs[y]<1&&Abs[1-x]<1;


Subscript[AC, 8][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1/y)^n (-(x/y))^m (-y)^-a Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m])+(x^m (-(1/y))^n (-y)^-b2 Gamma[a-b2] Gamma[c2] Pochhammer[b1,m] Pochhammer[a-b2,m-n] Pochhammer[b2,n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b2+c2] Pochhammer[c1,m]);


Subscript[Roc, 8][x_,y_]:=Abs[x]<1&&1/Abs[y]<1&&(1+Abs[x])/Abs[y]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 9][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1/x)^m (-x)^-a (-(y/x))^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,m+n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n])+((-(1/x))^m (-x)^-b1 y^n Gamma[a-b1] Gamma[c1] Pochhammer[a-b1,-m+n] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c1,m])/(m! n! Gamma[a] Gamma[-b1+c1] Pochhammer[c2,n]);


(* ::Input::Initialization:: *)
Subscript[Roc, 9][x_,y_]:=Abs[y]<1&&1/Abs[x]<1&&(1+Abs[y])/Abs[x]<1;


Subscript[AC, 10][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=(x^m (-y)^(-a-m) y^-n Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m])+((-x)^-b1 x^-m (-y)^-b2 y^-n Gamma[a-b1-b2] Gamma[c1] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c1,m] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m+n])+((-x)^(-a+b2) x^(-m+n) (-y)^-b2 y^-n Gamma[a-b2] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] Pochhammer[a-b2,m-n] Pochhammer[b2,n] Pochhammer[1+a-b2-c1,m-n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a-b1-b2,m-n]);


Subscript[Roc, 10][x_,y_]:=1/Abs[x]<1&&Abs[x/y]<1&&Abs[x/y]<Abs[x]/(1+Abs[x])&&1/Abs[y]<1&&Abs[x/y]+1/Abs[y]<1;


Subscript[AC, 11][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((-x)^-b1 x^-m (-y)^(-a+b1) y^(m-n) Gamma[a-b1] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] Pochhammer[a-b1,-m+n] Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[1+a-b1-c2,-m+n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a+b1+c2] Pochhammer[1+a-b1-b2,-m+n])+((-x)^-b1 x^-m (-y)^-b2 y^-n Gamma[a-b1-b2] Gamma[c1] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c1,m] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m+n])+((-x)^(-a-n) x^-m y^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,m+n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n]);


Subscript[Roc, 11][x_,y_]:=1/Abs[y]<1&&Abs[y/x]<1&&Abs[y/x]<Abs[y]/(1+Abs[y])&&1/Abs[x]<1&&1/Abs[x]+Abs[y/x]<1;


(*ACs using ET1 are below from 12-22*)


(* ::Input::Initialization:: *)
Subscript[AC, 12][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^-a (x/(-1+x))^m (y/(1-x))^n Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[-b1+c1,m])/(m! n! Pochhammer[c1,m] Pochhammer[c2,n]);


(* ::Input::Initialization:: *)
Subscript[Roc, 12][x_,y_]:=Abs[x/(-1+x)]+Abs[y/(1-x)]<1&&Abs[x/(-1+x)]<1;


Subscript[AC, 13][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x)^-a (x/(-1+x))^m ((-1+x+y)/(-1+x))^n Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((1-x)^-a (x/(-1+x+y))^m ((-1+x+y)/(-1+x))^n Gamma[a+b2-c2] Gamma[c2] If[1-Re[x]+Re[y]>0,((-1+x+y)/(-1+x))^(-a-b2+c2),((-1+x)/(-1+x+y))^(a+b2-c2)] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m] Pochhammer[-a+c2,-m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Pochhammer[c1,m] Pochhammer[1-a-b2+c2,-m+n]);


Subscript[Roc, 13][x_,y_]:=Abs[x/(-1+x+y)]<1&&Abs[(-1+x+y)/(-1+x)]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 14][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1/(1-x))^m (1-x)^-b1 y^n Gamma[a-b1] Gamma[c1] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a-c1,n] Pochhammer[-a+c1,m-n])/(m! n! Gamma[a] Gamma[-b1+c1] Pochhammer[1-a+b1,m-n] Pochhammer[c2,n])+((1/(1-x))^m (1-x)^-a (y/(1-x))^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,n] Pochhammer[-b1+c1,m])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n]);


(* ::Input::Initialization:: *)
Subscript[Roc, 14][x_,y_]:=Abs[1-x]>1&&Abs[y]<1;


Subscript[AC, 15][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x)^-a (x/(-1+x))^m ((-1+x+y)/(-1+x))^n Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((1-x)^-a (x/(-1+x))^n ((-1+x+y)/x)^m Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a-b1-b2+c1+c2] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(-a-b2+c2),(-((-1+x+y)/x))^(a+b2-c2)] If[1-Re[x]+Re[y]>0,((-1+x+y)/(-1+x))^(-a-b2+c2),((-1+x)/(-1+x+y))^(a+b2-c2)] Pochhammer[b2,m] Pochhammer[a+b2-c2,m-n] Pochhammer[1+a+b2-c1-c2,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a-b2+c1+c2] Pochhammer[b2,m-n] Pochhammer[1+a+b1+b2-c1-c2,m-n])+((1-x)^-a ((-1+x+y)/(-1+x))^n ((-1+x+y)/x)^m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] If[1-Re[x]+Re[y]>0,((-1+x+y)/(-1+x))^(-a-b2+c2),((-1+x)/(-1+x+y))^(a+b2-c2)] Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n] Pochhammer[-a-b1+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[-a-b1+c1+c2,m] Pochhammer[1-a-b1-b2+c1+c2,m+n]);


Subscript[Roc, 15][x_,y_]:=Abs[(-1+x+y)/x]<1&&Abs[x/(-1+x)]<1;


Subscript[AC, 16][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x)^-b1 (1/y)^n (-y)^(-a+b1) (y/(1-x))^m Gamma[a-b1] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] Pochhammer[a-b1,-m+n] Pochhammer[b1,m] Pochhammer[-b1+c1,n] Pochhammer[1+a-b1-c2,-m+n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a+b1+c2] Pochhammer[1+a-b1-b2,-m+n] Pochhammer[-b1+c1,-m+n])+((1/(1-x))^m (1-x)^-b1 (1/y)^n (-y)^-b2 Gamma[a-b1-b2] Gamma[c1] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a+b2+c1,m+n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m+n] Pochhammer[-a+b2+c1,n])+((1/(1-x))^m (1-x)^-a (y/(1-x))^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,n] Pochhammer[-b1+c1,m])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n]);


Subscript[Roc, 16][x_,y_]:=Abs[y]>1&&Abs[y/(1-x)]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 17][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^-a (-((-1+x)/x))^a ((-1+x)/x)^m ((-1+x+y)/x)^n Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,m])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,m+n])+((1-x)^-a (-((-1+x)/x))^(-b1+c1) ((-1+x)/x)^m ((-1+x+y)/(-1+x))^n Gamma[a+b1-c1] Gamma[c1] Gamma[1+a+b1-c1-c2] Gamma[a+b1+b2-c1-c2] Gamma[c2] Gamma[1-a-b1-b2+c1+c2] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[-a-b1-b2+c1+c2,m-n])/(m! n! Gamma[a] Gamma[b1] Gamma[a+b1-c1-c2] Gamma[1+a+b1+b2-c1-c2] Gamma[-b2+c2] Gamma[1-a-b1+c1+c2] Pochhammer[1-a-b1+c1,m-n] Pochhammer[-a-b1+c1+c2,m])+((1-x)^-a (-((-1+x)/x))^(-b1+c1) ((-1+x+y)/(-1+x))^n ((-1+x+y)/x)^m Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] If[1-Re[x]+Re[y]>0,((-1+x+y)/(-1+x))^(-a-b1-b2+c1+c2),((-1+x)/(-1+x+y))^(a+b1+b2-c1-c2)] Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n] Pochhammer[-a-b1+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[-a-b1+c1+c2,m] Pochhammer[1-a-b1-b2+c1+c2,m+n]);


(* ::Input::Initialization:: *)
Subscript[Roc, 17][x_,y_]:=Abs[(-1+x)/x]<1&&Abs[(-1+x+y)/(-1+x)]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 18][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^-a (1/y)^m ((1-x)/y)^n Gamma[-a+b2] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-a,((-1+x)/y)^a] Pochhammer[a,m+n] Pochhammer[b1,n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m+n])+((1/(1-x))^m (1-x)^-a ((1-x)/y)^n Gamma[a-b2] Gamma[a-b1-b2] Gamma[1-a+b1+b2] Gamma[1+a-b2-c1] Gamma[c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[b2,n] Pochhammer[-a+b1+b2,-m+n] Pochhammer[-b1+c1,m] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[1+a-b1-b2] Gamma[a-b2-c1] Gamma[1-a+b2+c1] Gamma[-b2+c2] Pochhammer[1-a+b2,-m+n] Pochhammer[-a+b2+c1,n])+((1/(1-x))^m (1-x)^(-b1-b2) (1/y)^n Gamma[a-b1-b2] Gamma[c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[-a+b2+c1,m+n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m+n] Pochhammer[-a+b2+c1,n]);


(* ::Input::Initialization:: *)
Subscript[Roc, 18][x_,y_]:= Abs[1-x]>1&&Abs[(1-x)/y]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 19][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^-a ((1-x)/y)^n (x/y)^m Gamma[-a+b2] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-a,((-1+x)/y)^a] Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m])+((1-x)^-a (x/(-1+x))^m ((-1+x)/y)^n Gamma[a-b2] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[a-b2,m-n] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b2+c2] Pochhammer[c1,m]);


Subscript[Roc, 19][x_,y_]:= Abs[x/(-1+x)]<1&&Abs[(1-x)/y]<1&&(1+Abs[x/(-1+x)]) Abs[(1-x)/y]<1;


(* ::Input::Initialization:: *)
Subscript[AC, 20][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^-a (-((-1+x)/x))^a ((-1+x)/x)^m (y/x)^n Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,m+n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n])+((-1+1/x)^m (1-x)^-a (-((-1+x)/x))^(-b1+c1) (y/(1-x))^n Gamma[a+b1-c1] Gamma[c1] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[a+b1-c1,-m+n] Pochhammer[-b1+c1,m])/(m! n! Gamma[a] Gamma[b1] Pochhammer[c2,n]);


Subscript[Roc, 20][x_,y_]:= Abs[y/(1-x)]<1&&Abs[(-1+x)/x]<1&&Abs[(-1+x)/x] (1+Abs[y/(1-x)])<1;


(* ::Input::Initialization:: *)
Subscript[AC, 21][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:= ((1-x)^-a (x/(-1+x))^m (y/(1-x))^-n (y/(-1+x))^-m Gamma[-a+b2] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-a,((-1+x)/y)^a] Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m])+((1-x)^-a (-((-1+x)/x))^(a-b2) (x/(-1+x))^(-m+n) (y/(1-x))^-n Gamma[a-b2] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[a-b2,m-n] Pochhammer[b2,n] Pochhammer[1+a-b2-c1,m-n] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a+b1-b2-c1,m-n])+((1-x)^-a (-((-1+x)/x))^(-b1+c1) (x/(-1+x))^-m (y/(1-x))^-n Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m+n]);


Subscript[Roc, 21][x_,y_]:= 1/Abs[x/(-1+x)]<1&&Abs[x/y]<1&&Abs[x/(-1+x)] (-1+Abs[x/y])+Abs[x/y]<0&&1/Abs[y/(1-x)]<1&&Abs[x/y]+1/Abs[y/(1-x)]<1;


Subscript[AC, 22][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x)^-a (-((-1+x)/x))^(-b1+c1) (x/(-1+x))^-m (y/(1-x))^-n Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^-b2,((-1+x)/y)^b2] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[1+b2-c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m+n])+((1-x)^-a (-((-1+x)/x))^(-b1+c1) (x/(-1+x))^-m (y/(1-x))^(m-n) Gamma[a+b1-c1] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] If[1-Re[x]+Re[y]>0,(y/(-1+x))^(-a-b1+c1),((-1+x)/y)^(a+b1-c1)] Pochhammer[1-b1,m] Pochhammer[a+b1-c1,-m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a+b1-c1-c2,-m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b1+c1+c2] Pochhammer[1+a+b1-b2-c1,-m+n])+((1-x)^-a (-((-1+x)/x))^a (-(x/(-1+x)))^-n (x/(-1+x))^-m (y/(1-x))^n Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b2,n] Pochhammer[1+a-c1,m+n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n]);


Subscript[Roc, 22][x_,y_]:= 1/Abs[y/(1-x)]<1&&Abs[y/x]<1&&Abs[y/(1-x)] (-1+Abs[y/x])+Abs[y/x]<0&&1/Abs[x/(-1+x)]<1&&1/Abs[x/(-1+x)]+Abs[y/x]<1;


(*ACs obtained using ET2 are numbered from 23-33*)


Subscript[AC, 23][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((x/(1-y))^m (1-y)^-a (y/(-1+y))^n Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[-b2+c2,n])/(m! n! Pochhammer[c1,m] Pochhammer[c2,n]);


Subscript[Roc, 23][x_,y_]:= Abs[x/(1-y)]+Abs[y/(-1+y)]<1&&Abs[x/(1-y)]<1;


Subscript[AC, 24][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=(x^m (1/(1-y))^n (1-y)^-b2 Gamma[a-b2] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+a-c2,m] Pochhammer[-a+c2,-m+n])/(m! n! Gamma[a] Gamma[-b2+c2] Pochhammer[1-a+b2,-m+n] Pochhammer[c1,m])+((1/(1-y))^n (x/(1-y))^m (1-y)^-a Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c2,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m]);


Subscript[Roc, 24][x_,y_]:= Abs[x]<1&&Abs[1-y]>1;


Subscript[AC, 25][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-y)^-a (y/(-1+y))^n ((-1+x+y)/(-1+y))^m Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c1,n] Pochhammer[-b2+c2,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n])+((1-y)^-a (y/(-1+x+y))^n ((-1+x+y)/(-1+y))^m Gamma[a+b1-c1] Gamma[c1] If[1+Re[x]-Re[y]>0,((-1+x+y)/(-1+y))^(-a-b1+c1),((-1+y)/(-1+x+y))^(a+b1-c1)] Pochhammer[1+a-c1,n] Pochhammer[-a+c1,m-n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Pochhammer[1-a-b1+c1,m-n] Pochhammer[c2,n]);


Subscript[Roc, 25][x_,y_]:= Abs[y/(-1+x+y)]<1&&Abs[(-1+x+y)/(-1+y)]<1;


Subscript[AC, 26][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1/x)^m (-x)^-b1 (1/(1-y))^n (1-y)^-b2 Gamma[a-b1-b2] Gamma[c1] Gamma[c2] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1+c2,m+n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m+n] Pochhammer[-a+b1+c2,m])+((1/x)^m (-x)^(-a+b2) (x/(1-y))^n (1-y)^-b2 Gamma[a-b2] Gamma[-a+b1+b2] Gamma[c1] Gamma[c2] Pochhammer[a-b2,m-n] Pochhammer[b2,n] Pochhammer[1+a-b2-c1,m-n] Pochhammer[-b2+c2,m])/(m! n! Gamma[a] Gamma[b1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a-b1-b2,m-n] Pochhammer[-b2+c2,m-n])+((1/(1-y))^n (x/(1-y))^m (1-y)^-a Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c2,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m]);


Subscript[Roc, 26][x_,y_]:= Abs[x]>1&&Abs[x/(1-y)]<1;


Subscript[AC, 27][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-y)^-a (y/(-1+y))^m ((-1+x+y)/y)^n Gamma[a+b1-c1] Gamma[c1] Gamma[c2] Gamma[-a-b1-b2+c1+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(-a-b1+c1),(-((-1+x+y)/y))^(a+b1-c1)] If[1+Re[x]-Re[y]>0,((-1+x+y)/(-1+y))^(-a-b1+c1),((-1+y)/(-1+x+y))^(a+b1-c1)] Pochhammer[b1,n] Pochhammer[a+b1-c1,-m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a+b1-c1-c2,-m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Gamma[-a-b1+c1+c2] Pochhammer[b1,-m+n] Pochhammer[1+a+b1+b2-c1-c2,-m+n])+((1-y)^-a (y/(-1+y))^n ((-1+x+y)/(-1+y))^m Gamma[c1] Gamma[-a-b1+c1] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c1,n] Pochhammer[-b2+c2,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n])+((1-y)^-a ((-1+x+y)/(-1+y))^m ((-1+x+y)/y)^n Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] If[1+Re[x]-Re[y]>0,((-1+x+y)/(-1+y))^(-a-b1+c1),((-1+y)/(-1+x+y))^(a+b1-c1)] Pochhammer[1-b2,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n] Pochhammer[-a-b2+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[-a-b2+c1+c2,n] Pochhammer[1-a-b1-b2+c1+c2,m+n]);


Subscript[Roc, 27][x_,y_]:= Abs[(-1+x+y)/y]<1&&Abs[y/(-1+y)]<1;


Subscript[AC, 28][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1/x)^m (1/(1-y))^n (1-y)^(-b1-b2) Gamma[a-b1-b2] Gamma[c1] Gamma[c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[b1,m] Pochhammer[b2,n] Pochhammer[1+b1-c1,m] Pochhammer[-a+b1+c2,m+n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-b2+c2] Pochhammer[1-a+b1+b2,m+n] Pochhammer[-a+b1+c2,m])+((1/x)^n (1-y)^-a ((1-y)/x)^m Gamma[-a+b1] Gamma[c1] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-a,((-1+y)/x)^a] Pochhammer[a,m+n] Pochhammer[b2,m] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,m+n])+((1/(1-y))^n (1-y)^-a ((1-y)/x)^m Gamma[a-b1] Gamma[a-b1-b2] Gamma[1-a+b1+b2] Gamma[c1] Gamma[1+a-b1-c2] Gamma[c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[b1,m] Pochhammer[-a+b1+b2,m-n] Pochhammer[1+b1-c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[1+a-b1-b2] Gamma[b2] Gamma[-b1+c1] Gamma[a-b1-c2] Gamma[1-a+b1+c2] Pochhammer[1-a+b1,m-n] Pochhammer[-a+b1+c2,m]);


Subscript[Roc, 28][x_,y_]:= Abs[1-y]>1&&Abs[(1-y)/x]<1;


Subscript[AC, 29][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-y)^-a (-((-1+y)/y))^a ((-1+y)/y)^n ((-1+x+y)/y)^m Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[-b1+c1,n] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m+n] Pochhammer[1+a+b2-c2,m+n])+((1-y)^-a (-((-1+y)/y))^(-b2+c2) ((-1+y)/y)^n ((-1+x+y)/(-1+y))^m Gamma[c1] Gamma[a+b2-c2] Gamma[1+a+b2-c1-c2] Gamma[a+b1+b2-c1-c2] Gamma[c2] Gamma[1-a-b1-b2+c1+c2] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[-b2+c2,n] Pochhammer[-a-b1-b2+c1+c2,-m+n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[a+b2-c1-c2] Gamma[1+a+b1+b2-c1-c2] Gamma[1-a-b2+c1+c2] Pochhammer[1-a-b2+c2,-m+n] Pochhammer[-a-b2+c1+c2,n])+((1-y)^-a (-((-1+y)/y))^(-b2+c2) ((-1+x+y)/(-1+y))^m ((-1+x+y)/y)^n Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] If[1+Re[x]-Re[y]>0,((-1+x+y)/(-1+y))^(-a-b1-b2+c1+c2),((-1+y)/(-1+x+y))^(a+b1+b2-c1-c2)] Pochhammer[1-b2,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n] Pochhammer[-a-b2+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[-a-b2+c1+c2,n] Pochhammer[1-a-b1-b2+c1+c2,m+n]);


Subscript[Roc, 29][x_,y_]:=Abs[(-1+y)/y]<1&&Abs[(-1+x+y)/(-1+y)]<1;


Subscript[AC, 30][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-y)^-a (x/y)^m (-((-1+y)/y))^a ((-1+y)/y)^n Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((-1+1/y)^n (x/(1-y))^m (1-y)^-a (-((-1+y)/y))^(-b2+c2) Gamma[a+b2-c2] Gamma[c2] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[a+b2-c2,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Pochhammer[c1,m]);


Subscript[Roc, 30][x_,y_]:=Abs[x/(1-y)]<1&&Abs[(-1+y)/y]<1&&(1+Abs[x/(1-y)]) Abs[(-1+y)/y]<1;


Subscript[AC, 31][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-y)^-a ((1-y)/x)^m (y/x)^n Gamma[-a+b1] Gamma[c1] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-a,((-1+y)/x)^a] Pochhammer[a,m+n] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n])+((1-y)^-a ((-1+y)/x)^m (y/(-1+y))^n Gamma[a-b1] Gamma[c1] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[a-b1,-m+n] Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Pochhammer[c2,n]);


Subscript[Roc, 31][x_,y_]:=Abs[y/(-1+y)]<1&&Abs[(1-y)/x]<1&&Abs[(1-y)/x] (1+Abs[y/(-1+y)])<1;


Subscript[AC, 32][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((x/(1-y))^m (1-y)^-a (-((-1+y)/y))^a (-(y/(-1+y)))^-m (y/(-1+y))^-n Gamma[c2] Gamma[-a-b2+c2] Pochhammer[a,m+n] Pochhammer[b1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((x/(1-y))^(-m+n) (1-y)^-a (-((-1+y)/y))^(-b2+c2) (y/(-1+y))^-n Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a+b1-b2+c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^(-a-b2+c2),((-1+y)/x)^(a+b2-c2)] Pochhammer[1-b2,n] Pochhammer[a+b2-c2,m-n] Pochhammer[1+a+b2-c1-c2,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Gamma[-a-b2+c1+c2] Pochhammer[1+a-b1+b2-c2,m-n])+((x/(1-y))^-m (1-y)^-a (-((-1+y)/y))^(-b2+c2) (y/(-1+y))^-n Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[1+b1-c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m+n]);


Subscript[Roc, 32][x_,y_]:=1/Abs[x/(1-y)]<1&&Abs[x/y]<1&&Abs[x/(1-y)] (-1+Abs[x/y])+Abs[x/y]<0&&1/Abs[y/(-1+y)]<1&&Abs[x/y]+1/Abs[y/(-1+y)]<1;


Subscript[AC, 33][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((x/(1-y))^-m (1-y)^-a (-((-1+y)/y))^(a-b1) (y/(-1+y))^(m-n) Gamma[a-b1] Gamma[c1] Gamma[c2] Gamma[-a+b1-b2+c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[a-b1,-m+n] Pochhammer[b1,m] Pochhammer[1+b1-c1,m] Pochhammer[1+a-b1-c2,-m+n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-a+b1+c2] Gamma[-b2+c2] Pochhammer[1+a-b1+b2-c2,-m+n])+((x/(1-y))^-m (1-y)^-a (x/(-1+y))^-n (y/(-1+y))^n Gamma[-a+b1] Gamma[c1] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-a,((-1+y)/x)^a] Pochhammer[a,m+n] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n])+((x/(1-y))^-m (1-y)^-a (-((-1+y)/y))^(-b2+c2) (y/(-1+y))^-n Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] If[1+Re[x]-Re[y]>0,(x/(-1+y))^-b1,((-1+y)/x)^b1] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[1+b1-c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m+n]);


Subscript[Roc, 33][x_,y_]:=1/Abs[y/(-1+y)]<1&&Abs[y/x]<1&&Abs[y/x] (1+Abs[y/(-1+y)])<Abs[y/(-1+y)]&&1/Abs[x/(1-y)]<1&&1/Abs[x/(1-y)]+Abs[y/x]<1;


(*ACs obtained using ET3 are numbered as 34-44 *)


Subscript[AC, 34][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (x/(-1+x+y))^m (y/(-1+x+y))^n Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Pochhammer[c1,m] Pochhammer[c2,n]);


Subscript[Roc, 34][x_,y_]:=Abs[x/(-1+x+y)]+Abs[y/(-1+x+y)]<1&&Abs[x/(-1+x+y)]<1;


Subscript[AC, 35][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((x/(-1+x))^m (1-x-y)^-a ((-1+x)/(-1+x+y))^n Gamma[a-b2] Gamma[c2] If[-1+Re[x]-Re[y]>0,((-1+x)/(-1+x+y))^(-a+b2),((-1+x+y)/(-1+x))^(a-b2)] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m] Pochhammer[-a+c2,-m+n])/(m! n! Gamma[a] Gamma[-b2+c2] Pochhammer[1-a+b2,-m+n] Pochhammer[c1,m])+((1-x-y)^-a ((-1+x)/(-1+x+y))^n (x/(-1+x+y))^m Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m]);


Subscript[Roc, 35][x_,y_]:=Abs[x/(-1+x)]<1&&Abs[(-1+x)/(-1+x+y)]<1;


Subscript[AC, 36][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (y/(-1+y))^n ((-1+y)/(-1+x+y))^m Gamma[a-b1] Gamma[c1] If[-1-Re[x]+Re[y]>0,((-1+y)/(-1+x+y))^(-a+b1),((-1+x+y)/(-1+y))^(a-b1)] Pochhammer[b1,m] Pochhammer[1+a-c1,n] Pochhammer[-a+c1,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[-b1+c1] Pochhammer[1-a+b1,m-n] Pochhammer[c2,n])+((1-x-y)^-a ((-1+y)/(-1+x+y))^m (y/(-1+x+y))^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[1+a-c1,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n]);


Subscript[Roc, 36][x_,y_]:=Abs[y/(-1+y)]<1&&Abs[(-1+y)/(-1+x+y)]<1;


Subscript[AC, 37][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((-((-1+x)/x))^(a-b2) ((-1+x)/x)^m (1-x-y)^-a (x/(-1+x+y))^n Gamma[a-b2] Gamma[c1] Gamma[-a-b1+b2+c1] Gamma[c2] If[-1+Re[x]-Re[y]>0,((-1+x)/(-1+x+y))^(-a+b2),((-1+x+y)/(-1+x))^(a-b2)] Pochhammer[a-b2,m-n] Pochhammer[b2,n] Pochhammer[1+a-b2-c1,m-n] Pochhammer[-b2+c2,m])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-a+b2+c1] Gamma[-b2+c2] Pochhammer[1+a+b1-b2-c1,m-n] Pochhammer[-b2+c2,m-n])+((1-x-y)^-a ((-1+x)/(-1+x+y))^n (x/(-1+x+y))^m Gamma[-a+b2] Gamma[c2] Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[b2] Gamma[-a+c2] Pochhammer[1+a-b2,m+n] Pochhammer[c1,m])+((-((-1+x)/x))^(-b1+c1) ((-1+x)/x)^m (1-x-y)^-a ((-1+x)/(-1+x+y))^n Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] If[-1+Re[x]-Re[y]>0,((-1+x)/(-1+x+y))^(-a+b2),((-1+x+y)/(-1+x))^(a-b2)] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[-a-b1+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m+n] Pochhammer[-a-b1+c1+c2,m]);


Subscript[Roc, 37][x_,y_]:=Abs[(-1+x)/x]<1&&Abs[x/(-1+x+y)]<1;


Subscript[AC, 38][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (-((-1+y)/y))^(a-b1) ((-1+y)/y)^n (y/(-1+x+y))^m Gamma[a-b1] Gamma[c1] Gamma[c2] Gamma[-a+b1-b2+c2] If[-1-Re[x]+Re[y]>0,((-1+y)/(-1+x+y))^(-a+b1),((-1+x+y)/(-1+y))^(a-b1)] Pochhammer[a-b1,-m+n] Pochhammer[b1,m] Pochhammer[-b1+c1,n] Pochhammer[1+a-b1-c2,-m+n])/(m! n! Gamma[a] Gamma[-b1+c1] Gamma[-a+b1+c2] Gamma[-b2+c2] Pochhammer[-b1+c1,-m+n] Pochhammer[1+a-b1+b2-c2,-m+n])+((1-x-y)^-a ((-1+y)/(-1+x+y))^m (y/(-1+x+y))^n Gamma[-a+b1] Gamma[c1] Pochhammer[a,m+n] Pochhammer[1+a-c1,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[b1] Gamma[-a+c1] Pochhammer[1+a-b1,m+n] Pochhammer[c2,n])+((1-x-y)^-a (-((-1+y)/y))^(-b2+c2) ((-1+y)/y)^n ((-1+y)/(-1+x+y))^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] If[-1-Re[x]+Re[y]>0,((-1+y)/(-1+x+y))^(-a+b1),((-1+x+y)/(-1+y))^(a-b1)] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[-b2+c2,n] Pochhammer[-a-b2+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m+n] Pochhammer[-a-b2+c1+c2,n]);


Subscript[Roc, 38][x_,y_]:=Abs[(-1+y)/y]<1&&Abs[y/(-1+x+y)]<1;


Subscript[AC, 39][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=(((-1+x)/x)^n (1-x-y)^-a ((-1+x+y)/x)^m Gamma[c1] Gamma[-a-b1+c1] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^-a,(-((-1+x+y)/x))^a] Pochhammer[a,m+n] Pochhammer[b2,m] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,m+n])+((1-x-y)^-a ((-1+x)/(-1+x+y))^n ((-1+x+y)/x)^m Gamma[a+b1-c1] Gamma[a+b1-b2-c1] Gamma[c1] Gamma[1-a-b1+b2+c1] Gamma[1+a+b1-c1-c2] Gamma[c2] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] Pochhammer[1-b1,m] Pochhammer[-b1+c1,m] Pochhammer[-a-b1+b2+c1,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Gamma[1+a+b1-b2-c1] Gamma[a+b1-c1-c2] Gamma[1-a-b1+c1+c2] Pochhammer[1-a-b1+c1,m-n] Pochhammer[-a-b1+c1+c2,m])+(((-1+x)/x)^m (1-x-y)^-a ((-1+x)/(-1+x+y))^n Gamma[a+b1-b2-c1] Gamma[c1] Gamma[c2] If[-1+Re[x]-Re[y]>0,((-1+x)/(-1+x+y))^(-a-b1+b2+c1),((-1+x+y)/(-1+x))^(a+b1-b2-c1)] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] Pochhammer[1-b1,m] Pochhammer[b2,n] Pochhammer[-b1+c1,m] Pochhammer[-a-b1+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Pochhammer[1-a-b1+b2+c1,m+n] Pochhammer[-a-b1+c1+c2,m]);


Subscript[Roc, 39][x_,y_]:=Abs[(-1+x+y)/x]<1&&Abs[(-1+x)/(-1+x+y)]<1;


Subscript[AC, 40][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a ((-1+y)/y)^m ((-1+x+y)/y)^n Gamma[c2] Gamma[-a-b2+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^-a,(-((-1+x+y)/y))^a] Pochhammer[a,m+n] Pochhammer[b1,n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m+n] Pochhammer[1+a+b2-c2,m+n])+((1-x-y)^-a ((-1+y)/(-1+x+y))^m ((-1+x+y)/y)^n Gamma[c1] Gamma[a+b2-c2] Gamma[a-b1+b2-c2] Gamma[1+a+b2-c1-c2] Gamma[c2] Gamma[1-a+b1-b2+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] Pochhammer[1-b2,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n] Pochhammer[-a+b1-b2+c2,-m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Gamma[1+a-b1+b2-c2] Gamma[a+b2-c1-c2] Gamma[1-a-b2+c1+c2] Pochhammer[1-a-b2+c2,-m+n] Pochhammer[-a-b2+c1+c2,n])+((1-x-y)^-a ((-1+y)/y)^n ((-1+y)/(-1+x+y))^m Gamma[c1] Gamma[a-b1+b2-c2] Gamma[c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] If[-1-Re[x]+Re[y]>0,((-1+y)/(-1+x+y))^(-a+b1-b2+c2),((-1+x+y)/(-1+y))^(a-b1+b2-c2)] Pochhammer[b1,m] Pochhammer[1-b2,n] Pochhammer[-b2+c2,n] Pochhammer[-a-b2+c1+c2,m+n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Pochhammer[1-a+b1-b2+c2,m+n] Pochhammer[-a-b2+c1+c2,n]);


Subscript[Roc, 40][x_,y_]:=Abs[(-1+x+y)/y]<1&&Abs[(-1+y)/(-1+x+y)]<1;


Subscript[AC, 41][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (-(x/y))^m ((-1+x+y)/y)^n Gamma[c2] Gamma[-a-b2+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^-a,(-((-1+x+y)/y))^a] Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((1-x-y)^-a (x/(-1+x+y))^m (-((-1+x+y)/y))^n Gamma[a+b2-c2] Gamma[c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] Pochhammer[1-b2,n] Pochhammer[-b1+c1,m] Pochhammer[a+b2-c2,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Pochhammer[c1,m]);


Subscript[Roc, 41][x_,y_]:=Abs[x/(-1+x+y)]<1&&Abs[(-1+x+y)/y]<1&&(1+Abs[x/(-1+x+y)]) Abs[(-1+x+y)/y]<1;


Subscript[AC, 42][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (-(y/x))^n ((-1+x+y)/x)^m Gamma[c1] Gamma[-a-b1+c1] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^-a,(-((-1+x+y)/x))^a] Pochhammer[a,m+n] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n])+((1-x-y)^-a (y/(-1+x+y))^n (-((-1+x+y)/x))^m Gamma[a+b1-c1] Gamma[c1] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] Pochhammer[1-b1,m] Pochhammer[a+b1-c1,-m+n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Pochhammer[c2,n]);


Subscript[Roc, 42][x_,y_]:=Abs[y/(-1+x+y)]<1&&Abs[(-1+x+y)/x]<1&&(1+Abs[y/(-1+x+y)]) Abs[(-1+x+y)/x]<1;


Subscript[AC, 43][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (x/(-1+x+y))^m (-(y/(-1+x+y)))^-m (y/(-1+x+y))^-n Gamma[c2] Gamma[-a-b2+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^-a,(-((-1+x+y)/y))^a] Pochhammer[a,m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a-c2,m+n])/(m! n! Gamma[-a+c2] Gamma[-b2+c2] Pochhammer[c1,m] Pochhammer[1+a+b2-c2,m+n])+((1-x-y)^-a (x/(-1+x+y))^(-m+n) (y/(-1+x+y))^-n Gamma[c1] Gamma[a+b2-c2] Gamma[c2] Gamma[-a-b1-b2+c1+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(-a-b2+c2),(-((-1+x+y)/x))^(a+b2-c2)] Pochhammer[1-b2,n] Pochhammer[a+b2-c2,m-n] Pochhammer[1+a+b2-c1-c2,m-n] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b2] Gamma[-b1+c1] Gamma[-a-b2+c1+c2] Pochhammer[1+a+b1+b2-c1-c2,m-n])+((1-x-y)^-a (x/(-1+x+y))^-m (y/(-1+x+y))^-n Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] Pochhammer[1-b1,m] Pochhammer[1-b2,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m+n]);


Subscript[Roc, 43][x_,y_]:=1/Abs[x/(-1+x+y)]<1&&Abs[x/y]<1&&Abs[x/y] (1+Abs[x/(-1+x+y)])<Abs[x/(-1+x+y)]&&1/Abs[y/(-1+x+y)]<1&&Abs[x/y]+1/Abs[y/(-1+x+y)]<1;


Subscript[AC, 44][a_,b1_,b2_,c1_,c2_,x_,y_,m_,n_]:=((1-x-y)^-a (x/(-1+x+y))^-m (y/(-1+x+y))^(m-n) Gamma[a+b1-c1] Gamma[c1] Gamma[c2] Gamma[-a-b1-b2+c1+c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(-a-b1+c1),(-((-1+x+y)/y))^(a+b1-c1)] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] Pochhammer[1-b1,m] Pochhammer[a+b1-c1,-m+n] Pochhammer[-b1+c1,m] Pochhammer[1+a+b1-c1-c2,-m+n])/(m! n! Gamma[a] Gamma[b1] Gamma[-b2+c2] Gamma[-a-b1+c1+c2] Pochhammer[1+a+b1+b2-c1-c2,-m+n])+((1-x-y)^-a (-(x/(-1+x+y)))^-n (x/(-1+x+y))^-m (y/(-1+x+y))^n Gamma[c1] Gamma[-a-b1+c1] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^-a,(-((-1+x+y)/x))^a] Pochhammer[a,m+n] Pochhammer[1+a-c1,m+n] Pochhammer[-b2+c2,n])/(m! n! Gamma[-a+c1] Gamma[-b1+c1] Pochhammer[1+a+b1-c1,m+n] Pochhammer[c2,n])+((1-x-y)^-a (x/(-1+x+y))^-m (y/(-1+x+y))^-n Gamma[c1] Gamma[a+b1+b2-c1-c2] Gamma[c2] If[-1+Re[x]-Re[y]>0,(-(y/(-1+x+y)))^(b2-c2),(-((-1+x+y)/y))^(-b2+c2)] If[-1-Re[x]+Re[y]>0,(-(x/(-1+x+y)))^(b1-c1),(-((-1+x+y)/x))^(-b1+c1)] Pochhammer[1-b1,m] Pochhammer[1-b2,n] Pochhammer[-b1+c1,m] Pochhammer[-b2+c2,n])/(m! n! Gamma[a] Gamma[b1] Gamma[b2] Pochhammer[1-a-b1-b2+c1+c2,m+n]);


Subscript[Roc, 44][x_,y_]:=1/Abs[y/(-1+x+y)]<1&&Abs[y/x]<1&&Abs[y/x] (1+Abs[y/(-1+x+y)])<Abs[y/(-1+x+y)]&&1/Abs[x/(-1+x+y)]<1&&Abs[y/x]+1/Abs[x/(-1+x+y)]<1;


End[]


EndPackage[]
