clear[standardBase];
makeSimpleRootSystem[A,rank_Integer]:=standardBase @@@ Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank},{j,1,rank+1}];
makeSimpleRootSystem[B,rank_Integer]:=standardBase @@@ Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],1]];
makeSimpleRootSystem[C,rank_Integer]:=standardBase @@@ Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],2]];
makeSimpleRootSystem[D,rank_Integer]:=standardBase @@@ Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Append[Table[0,{rank-2}],1],1]];
Unprotect[Dot];
Dot[x_standardBase,y_standardBase]:=standardBase @@ Dot[List@@x,List@@y];
Unprotect[Plus];
Plus[x_standardBase,y_standardBase]:=standardBase @@ Plus[List@@x,List@@y];
(* Unprotect[Power]; *)
(* Power[x_standardBase,y_Integer]:=standardBase @@ Power[List@@x,y]; *)
Unprotect[Times];
Times[x_,y_standardBase]:=standardBase @@ Times[x,List@@y];
reflection[x_standardBase]:=Function[y, y-2*(x.y)/(x.x)*x];
coroot[x_standardBase]:=2*x/(x.x);
cartanMatrix[{x__standardBase}]:=Transpose[Outer[Dot,{x},coroot/@{x}]];
clear[weylGroupElement];
revApply[x_,f_]:=f[x];
weylGroupElement[x__Integer][{y__standardBase}]:=Function[z,Fold[revApply,z,reflection /@ Part[{y},{x}]]];
fundamentalWeights[{simpleRoots__standardBase}]:=Plus@@(Inverse[cartanMatrix[{simpleRoots}]]*{simpleRoots});
rho[{simpleRoots__standardBase}]:=Plus@@fundamentalWeights[{simpleRoots}];
toFundamentalChamber[{simpleRoots__standardBase}][vec_standardBase]:=
    First[NestWhile[Function[v,
		       reflection[Scan[If[#.v<0,Return[#]]&,{simpleRoots}]][v]],
	      vec,
	      Head[#]=!=reflection[Null]&]];
orbit[{simpleRoots__standardBase}][{weights__standardBase}]:=
    Module[{total={}},
	   NestWhileList[
	       Function[x,
			Union[Flatten[Map[Function[y,
						   Map[reflection[#][y]&,Cases[{simpleRoots},z_ /; z.y>0]]],x]]]],
	       {weights},
	       Function[y,
			Module[{t=y=!={}},total=Union[total,y];t]
		       ]]];
orbit[{simpleRoots__standardBase}][weight_standardBase]:=orbit[{simpleRoots}][{toFundamentalChamber[{simpleRoots}][weight]}];
positiveRoots[{simpleRoots__standardBase}]:=Map[-#&,Flatten[orbit[{simpleRoots}][Map[-#&,{simpleRoots}]]]];
weightSystem[{simpleRoots__standardBase}][higestWeight_standardBase]:=Module[{minusPosRoots=-positiveRoots[{simpleRoots}]},
									     NestWhileList[Function[x,Complement[
										 Cases[Flatten[Outer[Plus,minusPosRoots,x]],y_/;And@@(#.y>=0&/@{simpleRoots})]
										 ,x]],{higestWeight},#=!={}&]]

freudenthalMultiplicities[{simpleRoots__standardBase}][highestWeight_standardBase]:=
    Module[{weights=Rest[Flatten[weightSystem[{simpleRoots}][highestWeight]]],
	    posroots=positiveRoots[{simpleRoots}],
	    rh=rho[{simpleRoots}],
	    mults,
	    toFC=toFundamentalChamber[{simpleRoots}]},
	   mults[highestWeight]=1;
	   Print[weights];
	   Scan[Function[v,
			 mults[v]=
			 2/((highestWeight+rh).(highestWeight+rh)-(v+rh).(v+rh))*
			 Plus@@
			     Map[Function[r,
					  Plus@@Map[mults[toFC[#[[1]]]]*#[[2]]&,
						    Rest[NestWhileList[{#[[1]]+r,#[[2]]+r.r}&,
								       {v,v.r},
								       IntegerQ[mults[toFC[#[[1]]+r]]]&]]]]
				 ,posroots];
			 Print[mults[v]];
			],
		weights];
	   mults]



freudenthal[{simpleRoots__standardBase}][highestWeight_standardBase]:=
    

Out[135]= {}

NestWhileList[#+1&,1,#+1<1&]

Out[143]= {1}

Out[117]= {1, 2, 3}

Out[116]= {1}

rho[b2]

                       3  1
Out[173]= standardBase[-, -]
                       2  2

fundamentalWeights[b2]

                                            1  1
Out[172]= {standardBase[1, 0], standardBase[-, -]}
                                            2  2

                                            1  1
Out[169]= {standardBase[1, 0], standardBase[-, -]}
                                            2  2

rho[b2]
                                            1  1
Out[168]= {standardBase[1, 0], standardBase[-, -]}
                                            2  2

                              
Out[167]= {standardBase[1, -1], standardBase[0, 1]}

Out[166]= standardBase[1, 0]

mts=freudenthalMultiplicities[b2][standardBase[2,0]]

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
2
-
3
1
4
-
5

Out[11]= mults$91

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
0
0
0

Out[4]= mults$81

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{{standardBase[1, 0], 1}}
{{standardBase[1, 0], 0}}
{{standardBase[1, 0], 1}}
{{standardBase[1, 0], 1}, {standardBase[2, 0], 2}}
2
-
3
{{standardBase[1, 1], 0}, {standardBase[2, 0], 2}}
{{standardBase[1, 1], 1}}
{{standardBase[1, 1], 2}}
{{standardBase[1, 1], 1}}
1
{{standardBase[0, 0], 0}, {standardBase[1, -1], 2}}
{{standardBase[0, 0], 0}}
{{standardBase[0, 0], 0}, {standardBase[1, 1], 2}}
{{standardBase[0, 0], 0}}
4
-
5

Out[174]= mults$1121

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{{standardBase[1, 0], 1}}
{{standardBase[1, 0], 0}}
{{standardBase[1, 0], 1}}
{{standardBase[1, 0], 1}, {standardBase[2, 0], 2}}
4
-
5
{{standardBase[1, 1], 0}, {standardBase[2, 0], 2}}
{{standardBase[1, 1], 1}}
{{standardBase[1, 1], 2}}
{{standardBase[1, 1], 1}}
1
{{standardBase[0, 0], 0}, {standardBase[1, -1], 2}}
{{standardBase[0, 0], 0}}
{{standardBase[0, 0], 0}, {standardBase[1, 1], 2}}
{{standardBase[0, 0], 0}}
1

Out[163]= mults$1097

mts[standardBase[0,0]]

Out[9]= 0

Out[8]= 0

Out[7]= 1

Out[6]= 0

Out[5]= 0

          2
Out[175]= -
          3


          4
Out[161]= -
          5

Out[160]= mults$1075[standardBase[0, 1]]

Out[159]= 1

Out[158]= mults$1075[standardBase[2, 1]]

Out[157]= 1

Out[156]= 1

Out[155]= {HoldPattern[mts] :> mults$1075}

Out[154]= {}

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[1, 0], 1}
mults$1075[standardBase[2, 1]]
False
{standardBase[1, 0], 0}
mults$1075[standardBase[1, 1]]
False
{standardBase[1, 0], 1}
mults$1075[standardBase[2, 1]]
False
{standardBase[1, 0], 1}
1
True
{standardBase[2, 0], 2}
mults$1075[standardBase[3, 0]]
False
4
-
5
{standardBase[1, 1], 0}
1
True
{standardBase[2, 0], 2}
mults$1075[standardBase[3, 1]]
False
{standardBase[1, 1], 1}
mults$1075[standardBase[2, 1]]
False
{standardBase[1, 1], 2}
mults$1075[standardBase[2, 2]]
False
{standardBase[1, 1], 1}
mults$1075[standardBase[2, 1]]
False
1
{standardBase[0, 0], 0}
1
True
{standardBase[1, -1], 2}
mults$1075[standardBase[2, 2]]
False
{standardBase[0, 0], 0}
4
-
5
False
{standardBase[0, 0], 0}
1
True
{standardBase[1, 1], 2}
mults$1075[standardBase[2, 2]]
False
{standardBase[0, 0], 0}
4
-
5
False
1

Out[153]= mults$1075

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[1, 0], 1}
mults$1053[standardBase[2, 1]]
False
{standardBase[1, 0], 0}
mults$1053[standardBase[1, 1]]
False
{standardBase[1, 0], 1}
mults$1053[standardBase[2, 1]]
False
{standardBase[1, 0], 1}
1
True
{standardBase[2, 0], 2}
mults$1053[standardBase[3, 0]]
False
          4
{0, 0, 0, -}
          5
{standardBase[1, 1], 0}
1
True
{standardBase[2, 0], 2}
mults$1053[standardBase[3, 1]]
False
{standardBase[1, 1], 1}
mults$1053[standardBase[2, 1]]
False
{standardBase[1, 1], 2}
mults$1053[standardBase[2, 2]]
False
{standardBase[1, 1], 1}
mults$1053[standardBase[2, 1]]
False
{1, 0, 0, 0}
{standardBase[0, 0], 0}
{1, 0, 0, 0}
False
{standardBase[0, 0], 0}
          4
{0, 0, 0, -}
          5
False
{standardBase[0, 0], 0}
{1, 0, 0, 0}
False
{standardBase[0, 0], 0}
          4
{0, 0, 0, -}
          5
False
{0, 0, 0, 0}

Out[151]= mults$1053

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
2 positiveRoots
---------------
       5
positiveRoots
-------------
      2
positiveRoots
-------------
      4

Out[149]= mults$1043

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[1, 0], 1}
mults$1027[standardBase[2, 1]]
False
{standardBase[1, 0], 0}
mults$1027[standardBase[1, 1]]
False
{0, 0}
{standardBase[1, 1], 0}
1
True
{standardBase[2, 0], 2}
mults$1027[standardBase[3, 1]]
False
{standardBase[1, 1], 1}
mults$1027[standardBase[2, 1]]
False
{1, 0}
{standardBase[0, 0], 0}
{1, 0}
False
{standardBase[0, 0], 0}
{0, 0}
False
{0, 0}

Out[147]= mults$1027

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[1, 0], 1}
mults$1011[standardBase[2, 1]]
{standardBase[1, 0], 0}
mults$1011[standardBase[1, 1]]
{0, 0}
{standardBase[1, 1], 0}
1
{standardBase[2, 0], 2}
mults$1011[standardBase[3, 1]]
{standardBase[1, 1], 1}
mults$1011[standardBase[2, 1]]
{1, 0}
{standardBase[0, 0], 0}
{1, 0}
{standardBase[0, 0], 0}
{0, 0}
{0, 0}

Out[145]= mults$1011

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[1, 0], 1}
False
{standardBase[1, 0], 0}
False
{0, 0}
{standardBase[1, 1], 0}
True
{standardBase[2, 0], 2}
False
{standardBase[1, 1], 1}
False
{1, 0}
{standardBase[0, 0], 0}
False
{standardBase[0, 0], 0}
False
{0, 0}

Out[142]= mults$995

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[1, 0], 1}
False
{standardBase[1, 0], 0}
False
{standardBase[1, 1], 0}
True
{standardBase[2, 0], 2}
False
{standardBase[1, 1], 1}
False
{standardBase[0, 0], 0}
False
{standardBase[0, 0], 0}
False

Out[140]= mults$979

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
standardBase[1, 0]
standardBase[1, -1]
{standardBase[2, -1], 3}
False
standardBase[1, 0]
standardBase[0, 1]
{standardBase[1, 1], 1}
False
standardBase[1, 1]
standardBase[1, -1]
{standardBase[2, 0], 2}
False
standardBase[1, 1]
standardBase[0, 1]
{standardBase[1, 2], 2}
False
standardBase[0, 0]
standardBase[1, -1]
{standardBase[1, -1], 2}
False
standardBase[0, 0]
standardBase[0, 1]
{standardBase[0, 1], 1}
True
{standardBase[0, 2], 2}
False

Out[134]= mults$963

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
standardBase[2, -1]
{standardBase[2, -1], 3}
False
standardBase[1, 1]
{standardBase[1, 1], 1}
False
standardBase[2, 0]
{standardBase[2, 0], 2}
False
standardBase[1, 2]
{standardBase[1, 2], 2}
False
standardBase[1, -1]
{standardBase[1, -1], 2}
False
standardBase[0, 1]
{standardBase[0, 1], 1}
True
{standardBase[0, 2], 2}
False

Out[132]= mults$947

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[2, -1], 3}
False
{standardBase[1, 1], 1}
False
{standardBase[2, 0], 2}
False
{standardBase[1, 2], 2}
False
{standardBase[1, -1], 2}
False
{standardBase[0, 1], 1}
True
{standardBase[0, 2], 2}
False

Out[130]= mults$931

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{standardBase[2, -1], 3}
{standardBase[1, 1], 1}
{standardBase[2, 0], 2}
{standardBase[1, 2], 2}
{standardBase[1, -1], 2}
{standardBase[0, 1], 1}
{standardBase[0, 2], 2}

Out[128]= mults$915

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
{{standardBase[2, -1], 3}}
{{standardBase[1, 1], 1}}
{{standardBase[2, 0], 2}}
{{standardBase[1, 2], 2}}
{{standardBase[1, -1], 2}}
{{standardBase[0, 1], 1}, {standardBase[0, 2], 2}}

Out[126]= mults$899

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
 6 mults$889[standardBase[2, 1]]  2 mults$889[standardBase[1, 1]]
{-------------------------------, -------------------------------}
                5                                5
{1, mults$889[standardBase[2, 1]]}
  1  mults$889[standardBase[2, 1]]
{{-, -----------------------------}, 
  2                2
 
         6 mults$889[standardBase[2, 1]]
     2 + -------------------------------
                        5
>   {-----------------------------------, 
                      4
 
             2 mults$889[standardBase[2, 1]]
         2 + -------------------------------
      3                     5
>    {-, -----------------------------------}}}
      5                   4

Out[124]= mults$889

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
 6 mults$879[standardBase[2, -1]]  2 mults$879[standardBase[1, 1]]
{--------------------------------, -------------------------------}
                5                                 5
{1, mults$879[standardBase[1, 2]]}
 mults$879[standardBase[1, -1]]
{------------------------------, 
               2
 
    mults$879[standardBase[0, 1]] + 2 mults$879[standardBase[0, 2]]
>   ---------------------------------------------------------------}
                                   4

Out[122]= mults$879

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}
 6 mults$869[standardBase[2, -1]]  2 mults$869[standardBase[1, 1]]
{--------------------------------, -------------------------------}
                5                                 5
{1, mults$869[standardBase[1, 2]]}
 mults$869[standardBase[1, -1]]
{------------------------------, 
               2
 
    mults$869[standardBase[0, 1]] + 2 mults$869[standardBase[0, 2]]
>   ---------------------------------------------------------------}
                                   4

Out[115]= mults$869


mts[standardBase[0,0]]

           mults$859[standardBase[1, -1]]
Out[112]= {------------------------------, 
                         2
 
     mults$859[standardBase[0, 1]] + 2 mults$859[standardBase[0, 2]]
>    ---------------------------------------------------------------}
                                    4

           6 mults$859[standardBase[2, -1]]
Out[111]= {--------------------------------, 
                          5
 
      2  2 mults$859[standardBase[1, 2]]
>    {-, -------------------------------}}
      5                 5

           2 mults$849[standardBase[2, -1]]
Out[108]= {--------------------------------, 0}
                          5

           2 mults$839[standardBase[2, -1]]
Out[105]= {--------------------------------, 0}
                          5

Out[104]= {0, 0}

Out[103]= {}

{standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}

Out[102]= mults$839

Module[{r=standardBase[0,1],v=standardBase[1,1],mults,toFC=toFundamentalChamber[b2]},
       mults[standardBase[2,1]]=1;
       NestWhileList[{#[[1]]+r,#[[2]]+r.r}&,
		     {v+r,v.r},
		     Module[{t=mults[toFC[First[#]]]},Print[t];IntegerQ[t]]&
		    ]]

                                                  1
mults$836[standardBase[3, 1]]

Out[100]= {{standardBase[1, 2], 1}, {standardBase[1, 3], 2}}

                                             mults$834[standardBase[2, 1]]

Out[99]= {{standardBase[1, 2], 1}}

                                             mults$832[standardBase[1, 1]]

Out[98]= {{standardBase[1, 1], 1}}

weightSystem[b2][standardBase[-5,0]]

Out[61]= {{standardBase[-5, 0]}, {}}

Out[60]= {{standardBase[5, 0]}, {standardBase[4, 0], standardBase[4, 1]}, 
 
>    {standardBase[3, 0], standardBase[3, 1], standardBase[3, 2]}, 
 
>    {standardBase[2, 0], standardBase[2, 1], standardBase[2, 2]}, 
 
>    {standardBase[1, 0], standardBase[1, 1]}, {standardBase[0, 0]}, {}}


Timing[Length[positiveRoots[makeSimpleRootSystem[B,64]]]]

Flatten[orbit[{simpleRoots}][simpleRoots]],x_ /; Dot[x,rho[{simpleRoots}]]>=0]


Map[-#&,Sequence[1,2,3]]
Out[54]= Map[-#1 & , 1, 2, 3]

Length[Flatten[orbit[makeSimpleRootSystem[B,8]][standardBase[-1,-1/2,0,0,0,0,0,0]]]]

rho[b2]

                      3  1
Out[98]= standardBase[-, -]
                      2  2

b2

Out[67]= {standardBase[1, -1], standardBase[0, 1]}

Union[Flatten[Outer[If[#1.#2>0,reflection[#1][#2],Sequence[]]&,
      b2,
      b2]]]


Print[expr] prints expr as output. 

toFundamentalChamber[b2][standardBase[0,-1]]

Out[63]= standardBase[1, 0]

Out[43]= standardBase[1, 0]

Out[42]= standardBase[1, 1]

Out[38]= Sequence[standardBase[1, 1]]

Out[36]= 

First[reflection[Null][standardBase[1, 1]]]

Out[40]= standardBase[1, 1]

Out[39]= standardBase[1, 1]

Out[32]= standardBase[-1, -1]

Function::flpar: 
   Parameter specification standardBase[-1, -1] in 
    Function[standardBase[-1, -1], 
     reflection[Scan[If[#1 . standardBase[-1, -1] < 0, Return[#1]] & , 
        {standardBase[1, -1], standardBase[0, 1]}]][standardBase[-1, -1]]]
     should be a symbol or a list of symbols.

Out[28]= standardBase[-1, -1]

Head[expr] gives the head of expr. 

Head[reflection[Null][1,2,3]]

Out[24]= reflection[Null]

v=standardBase[-1,-1]

Out[17]= standardBase[-1, -1]

reflection[Scan[If[#.v<0,Return[#]]&,b2]][v]

Out[30]= standardBase[-1, 1]

Out[29]= reflection[Null][standardBase[-1, -1]]

Out[22]= reflection[Null]





Timing[Length[Map[-#&,Flatten[orbit[makeSimpleRootSystem[B,64]][Sequence @@ Map[-#&,makeSimpleRootSystem[B,64]]]]]]]

Out[15]= {8.15651, 4096}

Out[12]= 1024

Out[7]= {standardBase[1, -1], standardBase[0, 1], standardBase[1, 1], 
 
>    standardBase[1, 0]}

Out[6]= {{standardBase[-1, 1], standardBase[0, -1]}, 
 
>    {standardBase[-1, -1], standardBase[-1, 0]}, {}}

b2

Out[5]= {standardBase[1, -1], standardBase[0, 1]}

Out[4]= {{standardBase[1, -1], standardBase[0, 1]}, 
 
>    {standardBase[-1, 1], standardBase[0, -1]}, 
 
>    {standardBase[-1, -1], standardBase[-1, 0]}, {}}


Out[8]= $Aborted

positiveRoots[{simpleRoots__standardBase}]:=Cases[Flatten[orbit[{simpleRoots}][simpleRoots]],x_ /; Dot[x,rho[{simpleRoots}]]>=0]

weightSystem[{simpleRoots__standardBase}][higestWeight_standardBase]:=Module[{minusPosRoots=-positiveRoots[{simpleRoots}]},
									     NestWhileList[Function[x,Complement[
										 Cases[Flatten[Outer[Plus,minusPosRoots,x]],y_/;And@@(#.y>=0&/@{simpleRoots})]
										 ,x]],{higestWeight},#=!={}&]]

D18=makeSimpleRootSystem[D,18];


Length[positiveRoots[D18]]

[Calculating...]

Out[192]= 56

          
Out[191]= 56

minusPosRoots=-positiveRoots[b2];
Map[Function[x,And@@Map[#.x>=0&,b2]],Flatten[Outer[Plus,minusPosRoots,{standardBase[2,0]}]]]

          
Out[183]= {True, False, True, False}

          
Out[182]= {{True, True}, {True, False}, {True, True}, {True, False}}

          
Out[181]= {{True, True}, {True, False}, {True, True}, {True, False}}

          
Out[176]= {{standardBase[1, -1] . {standardBase[1, 1]} >= 0, 
 
>     standardBase[0, 1] . {standardBase[1, 1]} >= 0}, 
 
>    {standardBase[1, -1] . {standardBase[2, -1]} >= 0, 
 
>     standardBase[0, 1] . {standardBase[2, -1]} >= 0}, 
 
>    {standardBase[1, -1] . {standardBase[1, 0]} >= 0, 
 
>     standardBase[0, 1] . {standardBase[1, 0]} >= 0}, 
 
>    {standardBase[1, -1] . {standardBase[1, -1]} >= 0, 
 
>     standardBase[0, 1] . {standardBase[1, -1]} >= 0}}

          
Out[175]= {}

          
Out[174]= {}

          
Out[173]= {}

          
Out[172]= {}

          
Out[171]= {{standardBase[1, 1]}, {standardBase[2, -1]}, {standardBase[1, 0]}, 
 
>    {standardBase[1, -1]}}

Out[170]= {standardBase[-1, 1], standardBase[0, -1], standardBase[-1, 0], 
 
>    standardBase[-1, -1]}

Flatten[weightSystem[b2][standardBase[4,2]]]

Out[189]= {standardBase[4, 2], standardBase[3, 1], standardBase[3, 2], 
 
>    standardBase[3, 3], standardBase[4, 1], standardBase[2, 0], 
 
>    standardBase[2, 1], standardBase[2, 2], standardBase[3, 0], 
 
>    standardBase[4, 0], standardBase[1, 0], standardBase[1, 1], 
 
>    standardBase[3, 1], standardBase[0, 0], standardBase[2, 0], 
 
>    standardBase[2, 1], standardBase[2, 2], standardBase[3, 0], 
 
>    standardBase[1, 0], standardBase[1, 1], standardBase[0, 0]}

Out[188]= {{standardBase[4, 2]}, 
 
>    {standardBase[3, 1], standardBase[3, 2], standardBase[3, 3], 
 
>     standardBase[4, 1]}, {standardBase[2, 0], standardBase[2, 1], 
 
>     standardBase[2, 2], standardBase[3, 0], standardBase[4, 0]}, 
 
>    {standardBase[1, 0], standardBase[1, 1], standardBase[3, 1]}, 
 
>    {standardBase[0, 0], standardBase[2, 0], standardBase[2, 1], 
 
>     standardBase[2, 2], standardBase[3, 0]}, 
 
>    {standardBase[1, 0], standardBase[1, 1]}, {standardBase[0, 0]}, {}}

Out[187]= {{standardBase[2, 2]}, {standardBase[1, 1], standardBase[2, 1]}, 
 
>    {standardBase[0, 0], standardBase[1, 0], standardBase[2, 0]}, 
 
>    {standardBase[1, 1]}, {standardBase[0, 0], standardBase[1, 0]}, {}}

Out[186]= {{standardBase[2, 4]}, {}}

Out[185]= {{standardBase[2, 0]}, {standardBase[1, 0], standardBase[1, 1]}, 
 
>    {standardBase[0, 0]}, {}}

Out[180]= {{standardBase[2, 0]}, {}}

Out[178]= {standardBase[-1, 1], standardBase[0, -1], standardBase[-1, 0], 
 
>    standardBase[-1, -1]}

Out[169]= {standardBase[-1, 1], standardBase[0, -1], standardBase[-1, 0], 
 
>    standardBase[-1, -1]}

Module::argrx: Module called with 3 arguments; 2 arguments are expected.

Out[167]= Module[{minusPosRoots$ = 
 
>      -positiveRoots[{standardBase[1, -1], standardBase[0, 1]}]}, 
 
>    NestWhileList[Function[x$, 
 
>      Complement[Cases[Outer[Plus, minusPosRoots$, x$], 
 
>        And[x$_ /; 
 
>          (#1 . x$ >= 0 & ) {standardBase[1, -1], standardBase[0, 1]}]], x$]]
 
>      , {standardBase[2, 0]}, #1 =!= {} & ], minusPosRoots$]

Out[165]= {{standardBase[2, 0]}, {}}

Out[164]= {{standardBase[1, 0]}, {}}

?Minus

-x is the arithmetic negation of x. 

positiveRoots[makeSimpleRootSystem[B,4]]

Out[156]= {standardBase[1, -1, 0, 0], standardBase[0, 1, -1, 0], 
 
>    standardBase[0, 0, 1, -1], standardBase[0, 0, 0, 1], 
 
>    standardBase[0, 0, 1, 0], standardBase[0, 0, 1, 1], 
 
>    standardBase[0, 1, 0, -1], standardBase[1, 0, -1, 0], 
 
>    standardBase[0, 1, 0, 0], standardBase[0, 1, 0, 1], 
 
>    standardBase[1, 0, 0, -1], standardBase[0, 1, 1, 0], 
 
>    standardBase[1, 0, 0, 0], standardBase[1, 0, 0, 1], 
 
>    standardBase[1, 0, 1, 0], standardBase[1, 1, 0, 0]}

Out[155]= {standardBase[1, -1], standardBase[0, 1], standardBase[1, 0], 
 
>    standardBase[1, 1]}

Apply::level: Level specification standardBase[0, 1]
     is not of the form n, {n}, or {m, n}.

Out[153]= {}

Apply::level: Level specification standardBase[0, 1]
     is not of the form n, {n}, or {m, n}.

Out[151]= {}

Apply::level: Level specification standardBase[0, 1]
     is not of the form n, {n}, or {m, n}.

Out[149]= {}

Cases[Flatten[orbit[b2][Sequence@@(-#&/@b2)]],x_ /; Dot[x,rho[b2]]>=0]

Out[147]= {standardBase[0, 1], standardBase[1, -1], standardBase[1, 0], 
 
>    standardBase[1, 1]}

Out[145]= {}

Out[144]= {}

Out[143]= {}

Out[142]= {{standardBase[-1, 1], standardBase[0, -1]}, 
 
>    {standardBase[-1, -1], standardBase[-1, 0], standardBase[0, 1], 
 
>     standardBase[1, -1]}, {standardBase[1, 0], standardBase[1, 1]}, {}}


rho[b2]

                       3  1
Out[141]= standardBase[-, -]
                       2  2

Out[140]= {}

?Cases

Cases[{e , e , ...}, pattern] gives a list of the e
        1   2                                      i
     that match the pattern. Cases[{e , ...}, pattern -> rhs]
                                     1
      gives a list of the values of rhs
       corresponding to the e  that match the pattern. 
                             i
        Cases[expr, pattern, levelspec]
         gives a list of all parts of expr
          on levels specified by levelspec
           that match the pattern. Cases[expr, pattern -> rhs, levelspec]

            gives the values of rhs
             that match the pattern. Cases[expr, pattern, levelspec, n]

              gives the first n parts in expr that match the pattern. 

Out[137]= Function[]

Out[136]= Filter[#1 . rho[b2] > 0 & , 
 
>    {{standardBase[-1, 1], standardBase[0, -1]}, 
 
>     {standardBase[-1, -1], standardBase[-1, 0], standardBase[0, 1], 
 
>      standardBase[1, -1]}, {standardBase[1, 0], standardBase[1, 1]}, {}}]



Out[135]= {{standardBase[-1, 1], standardBase[0, -1]}, 
 
>    {standardBase[-1, -1], standardBase[-1, 0], standardBase[0, 1], 
 
>     standardBase[1, -1]}, {standardBase[1, 0], standardBase[1, 1]}, {}}


b2

Sequence@@(-#&/@b2)

Out[134]= Sequence[standardBase[-1, 1], standardBase[0, -1]]

Out[132]= {standardBase[-1, 1], standardBase[0, -1]}

Out[131]= orbit[{standardBase[1, -1], standardBase[0, 1]}][{standardBase[-1, 
 
>      1], standardBase[0, -1]}]

Out[130]= orbit[{standardBase[1, -1], standardBase[0, 1]}][{standardBase[-1, 
 
>      1], standardBase[0, -1]}]

positiveRoots[{simpleRoots__standardBase}]:=Flatten[orbit

calcMults[hw_standardBase,{simpleRoots__standardBase}]:=Module[{mults[hw]:=1},
							       
multiplicity[x_standardBase]:=0


orbit[b2][standardBase[1,0]]

Out[65]= {{standardBase[1, 0]}, {}}

Out[10]= {{standardBase[1, 0]}, {}}

Out[2]= orbit[b2][standardBase[1, 0]]

Out[129]= {{standardBase[1, 0]}, {standardBase[0, 1]}, {standardBase[0, -1]}, 
 
>    {standardBase[-1, 0]}, {}}

??First

First[expr] gives the first element in expr. 

Attributes[First] = {Protected}



Information::notfound: Symbol Tail not found.

Last[expr] gives the last element in expr. 

First[expr] gives the first element in expr. 

          
                         3  1                  1  3                3    1
Out[123]= {{standardBase[-, -]}, {standardBase[-, -], standardBase[-, -(-)]}, 
                         2  2                  2  2                2    2
 
                     1   3                1    3
>    {standardBase[-(-), -], standardBase[-, -(-)]}, 
                     2   2                2    2
 
                     3   1                  1     3
>    {standardBase[-(-), -], standardBase[-(-), -(-)]}, 
                     2   2                  2     2
 
                     3     1
>    {standardBase[-(-), -(-)]}, {}}
                     2     2

          
                         3  1                  1  3                3    1
Out[122]= {{standardBase[-, -]}, {standardBase[-, -], standardBase[-, -(-)]}, 
                         2  2                  2  2                2    2
 
                     1   3                1    3
>    {standardBase[-(-), -], standardBase[-, -(-)]}, 
                     2   2                2    2
 
                     3   1                  1     3
>    {standardBase[-(-), -], standardBase[-(-), -(-)]}, 
                     2   2                  2     2
 
                     3     1
>    {standardBase[-(-), -(-)]}, {}}
                     2     2

?SetMinus

Information::notfound: Symbol SetMinus not found.

          
                         3  1                  1  3                3    1
Out[120]= {{standardBase[-, -]}, {standardBase[-, -], standardBase[-, -(-)]}, 
                         2  2                  2  2                2    2
 
                     1   3                1    3                 3  1
>    {standardBase[-(-), -], standardBase[-, -(-)], standardBase[-, -]}, 
                     2   2                2    2                 2  2
 
                     3   1                  1     3                 1  3
>    {standardBase[-(-), -], standardBase[-(-), -(-)], standardBase[-, -], 
                     2   2                  2     2                 2  2
 
                   3    1                     3     1
>     standardBase[-, -(-)]}, {standardBase[-(-), -(-)], 
                   2    2                     2     2
 
                     1   3                1    3                 3  1
>     standardBase[-(-), -], standardBase[-, -(-)], standardBase[-, -]}, 
                     2   2                2    2                 2  2
 
                     3   1                  1     3                 1  3
>    {standardBase[-(-), -], standardBase[-(-), -(-)], standardBase[-, -], 
                     2   2                  2     2                 2  2
 
                   3    1
>     standardBase[-, -(-)]}}
                   2    2

          
                         3  1                  1  3                3    1
Out[119]= {{standardBase[-, -]}, {standardBase[-, -], standardBase[-, -(-)]}, 
                         2  2                  2  2                2    2
 
                     1   3                1    3                 3  1
>    {standardBase[-(-), -], standardBase[-, -(-)], standardBase[-, -]}, 
                     2   2                2    2                 2  2
 
                     3   1                  1     3                 1  3
>    {standardBase[-(-), -], standardBase[-(-), -(-)], standardBase[-, -], 
                     2   2                  2     2                 2  2
 
                   3    1
>     standardBase[-, -(-)]}}
                   2    2

Outer[Apply,reflection/@b2,{rho[b2]}]

            3                                        3
Out[118]= {{- + standardBase[-standardBase[1, -1] . (-), 
            2                                        2
 
                               3
>       standardBase[1, -1] . (-)]}, 
                               2
 
      3                                            3
>    {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}
      2                                            2

Outer[#1[#2]&,reflection/@b2,{rho[b2]}]

                         1  3                  3    1
Out[117]= {{standardBase[-, -]}, {standardBase[-, -(-)]}}
                         2  2                  2    2

?Apply

Apply[f, expr] or f @@ expr replaces the head of expr
      by f. Apply[f, expr, levelspec]

        replaces heads in parts of expr specified by levelspec. 

                         1  3                  3    1
Out[115]= {{standardBase[-, -]}, {standardBase[-, -(-)]}}
                         2  2                  2    2

                          1  3                    3    1
Out[114]= {{{standardBase[-, -]}}, {{standardBase[-, -(-)]}}}
                          2  2                    2    2

Out[113]= {{{Function[y$, y$ - 
 
         2 standardBase[1, -1] . y$ standardBase[1, -1]                3  1
>        ----------------------------------------------], standardBase[-, -]}}
           standardBase[1, -1] . standardBase[1, -1]                   2  2
 
                            2 standardBase[0, 1] . y$ standardBase[0, 1]
>     , {{Function[y$, y$ - --------------------------------------------], 
                              standardBase[0, 1] . standardBase[0, 1]
 
                    3  1
>      standardBase[-, -]}}}
                    2  2

            3                                        3
Out[111]= {{- + standardBase[-standardBase[1, -1] . (-), 
            2                                        2
 
                               3
>       standardBase[1, -1] . (-)]}, 
                               2
 
      3                                            3
>    {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}
      2                                            2

                             2 standardBase[1, -1] . y$ standardBase[1, -1]
Out[110]= {Function[y$, y$ - ----------------------------------------------], 
                               standardBase[1, -1] . standardBase[1, -1]
 
                       2 standardBase[0, 1] . y$ standardBase[0, 1]
>    Function[y$, y$ - --------------------------------------------]}
                         standardBase[0, 1] . standardBase[0, 1]
          
                         3  1
Out[109]= {{standardBase[-, -]}, 
                         2  2
 
       3                                            3
>    {{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
       2                                            2
 
       3                                        3
>     {- + standardBase[-standardBase[1, -1] . (-), 
       2                                        2
 
                                3
>        standardBase[1, -1] . (-)]}}, 
                                2
 
        3                                            3
>    {{{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
        2                                            2
 
        3                                            3
>      {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}, 
        2                                            2
 
        3                                        3
>     {{- + standardBase[-standardBase[1, -1] . (-), 
        2                                        2
 
                                 3
>         standardBase[1, -1] . (-)]}, 
                                 2
 
        3                                        3
>      {- + standardBase[-standardBase[1, -1] . (-), 
        2                                        2
 
                                 3
>         standardBase[1, -1] . (-)]}}}, 
                                 2
 
         3                                            3
>    {{{{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
         2                                            2
 
         3                                            3
>       {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}, 
         2                                            2
 
         3                                            3
>      {{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
         2                                            2
 
         3                                            3
>       {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}}, 
         2                                            2
 
         3                                        3
>     {{{- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}, 
                                  2
 
         3                                        3
>       {- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}}, 
                                  2
 
         3                                        3
>      {{- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}, 
                                  2
 
         3                                        3
>       {- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}}}}}
                                  2

          
Intersection::normal: 
   Nonatomic expression expected at position 2 in Intersection[{}, x].

Union::normal: Nonatomic expression expected at position 2 in Union[{}, x].

Intersection::normal: 
   Nonatomic expression expected at position 2 in 
    Intersection[Union[{}, x], x].

Union::normal: Nonatomic expression expected at position 2 in Union[{}, x, x].

Intersection::normal: 
   Nonatomic expression expected at position 2 in 
    Intersection[Union[{}, x, x], x].

General::stop: Further output of Intersection::normal
     will be suppressed during this calculation.

Union::normal: Nonatomic expression expected at position 2 in 
    Union[{}, x, x, x].

General::stop: Further output of Union::normal
     will be suppressed during this calculation.

                         3  1
Out[108]= {{standardBase[-, -]}, 
                         2  2
 
       3                                            3
>    {{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
       2                                            2
 
       3                                        3
>     {- + standardBase[-standardBase[1, -1] . (-), 
       2                                        2
 
                                3
>        standardBase[1, -1] . (-)]}}, 
                                2
 
        3                                            3
>    {{{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
        2                                            2
 
        3                                            3
>      {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}, 
        2                                            2
 
        3                                        3
>     {{- + standardBase[-standardBase[1, -1] . (-), 
        2                                        2
 
                                 3
>         standardBase[1, -1] . (-)]}, 
                                 2
 
        3                                        3
>      {- + standardBase[-standardBase[1, -1] . (-), 
        2                                        2
 
                                 3
>         standardBase[1, -1] . (-)]}}}, 
                                 2
 
         3                                            3
>    {{{{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
         2                                            2
 
         3                                            3
>       {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}, 
         2                                            2
 
         3                                            3
>      {{- + standardBase[0, -2 standardBase[0, 1] . (-)]}, 
         2                                            2
 
         3                                            3
>       {- + standardBase[0, -2 standardBase[0, 1] . (-)]}}}, 
         2                                            2
 
         3                                        3
>     {{{- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}, 
                                  2
 
         3                                        3
>       {- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}}, 
                                  2
 
         3                                        3
>      {{- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}, 
                                  2
 
         3                                        3
>       {- + standardBase[-standardBase[1, -1] . (-), 
         2                                        2
 
                                  3
>          standardBase[1, -1] . (-)]}}}}}
                                  2

          
Module::argrx: Module called with 4 arguments; 2 arguments are expected.

Out[107]= Module[{total = {}}, NestWhileList[Function[x, 
 
>      Union[Outer[Apply, reflection /@ b2, x]]], {rho[b2]}, 
 
>     Function[y, Module[{t = Intersection[total, x] =!= x}, 
 
>       total = Union[total, x]; t]]], 1, 3]

Out[106]= $Aborted

Union[{1,2,3}]=!=Union[{2,1,3}]

Out[105]= False

Out[104]= True

Out[103]= False

?Intersection

Intersection[list , list , ...] gives a sorted list of the elements common to
                 1      2
     all the list . 
                 i

?Union

Union[list , list , ...] gives a sorted list of all the distinct elements that
          1      2
     appear in any of the list . Union[list]
                              i
      gives a sorted version of a list, in which all duplicated elements have
      been dropped. 

?NestWhile

?NestWhileList

NestWhileList[f, expr, test] generates a list of the results of applying f
     repeatedly, starting with expr
     , and continuing until applying test
       to the result no longer yields True. NestWhileList[f, expr, test, m]

        supplies the most recent m
         results as arguments for test
          at each step. NestWhileList[f, expr, test, All]

           supplies all results so far as arguments for test
            at each step. NestWhileList[f, expr, test, m, max]

             applies f at most max times. 

NestWhile[f, expr, test] starts with expr
    , then repeatedly applies f until applying test
       to the result no longer yields True. NestWhile[f, expr, test, m]

        supplies the most recent m
         results as arguments for test
          at each step. NestWhile[f, expr, test, All]

           supplies all results so far as arguments for test
            at each step. NestWhile[f, expr, test, m, max]

             applies f at most max
               times. NestWhile[f, expr, test, m, max, n]

                applies f an extra n
                  times. NestWhile[f, expr, test, m, max, -n]

                   returns the result found when f
                    had been applied n fewer times. 

                      3  1
Out[95]= standardBase[-, -]
                      2  2

                      3  1
Out[94]= standardBase[-, -]
                      2  2

                      3  1
Out[92]= standardBase[-, -]
                      2  2

Out[91]= rho[{standardBase[1, -1], standardBase[0, 1]}]

fundamentalWeights[Sequence@@(makeSimpleRootSystem[D,5])]

Out[89]= {standardBase[1, 0, 0, 0, 0], standardBase[1, 1, 0, 0, 0], 
 
                                               1  1  1  1    1
>    standardBase[1, 1, 1, 0, 0], standardBase[-, -, -, -, -(-)], 
                                               2  2  2  2    2
 
                  1  1  1  1  1
>    standardBase[-, -, -, -, -]}
                  2  2  2  2  2

                                           1  1
Out[88]= {standardBase[1, 0], standardBase[-, -]}
                                           2  2

?Sequence

Sequence[expr , expr , ...] represents a sequence of arguments to be spliced
             1      2
    automatically into any function. 

Out[86]= fundamentalWeights[{standardBase[1, -1], standardBase[0, 1]}]


Plus@@(Inverse[cartanMatrix[b2]]*b2)

b2

Out[84]= {standardBase[1, -1], standardBase[0, 1]}

                                           1  1
Out[83]= {standardBase[1, 0], standardBase[-, -]}
                                           2  2

                                               3
Out[82]= {standardBase[2, -2], standardBase[0, -]}
                                               2

                       3    3
Out[81]= {standardBase[-, -(-)], standardBase[0, 2]}
                       2    2

Out[80]= {{standardBase[1, -1], standardBase[0, 1]}, 
 
                   1    1
>    {standardBase[-, -(-)], standardBase[0, 1]}}
                   2    2

                                             1    1
Out[74]= {{standardBase[1, -1], standardBase[-, -(-)]}, 
                                             2    2
 
>    {standardBase[0, 1], standardBase[0, 1]}}

              1
Out[73]= {{1, -}, {1, 1}}
              2

FullForm[b2.b2]

{{1,2},{3,4}}.{{1,2},{3,4}}

b2.b2

Out[72]= standardBase[1, 2]

{1,2}^2

Out[70]= {1, 4}

                           2                      2
Out[69]= standardBase[0, 1]  + standardBase[1, -1]

Out[68]= {{7, 10}, {15, 22}}

Out[67]= 5

??Power

x^y gives x to the power y. 

Attributes[Power] = {Listable, NumericFunction, OneIdentity, Protected}
 
Power /: Default[Power, 2] := 1

x^y gives x to the power y. 

Out[64]//FullForm= 
 
>   Plus[Power[standardBase[0, 1], 2], Power[standardBase[1, -1], 2]]

                           2                      2
Out[63]= standardBase[0, 1]  + standardBase[1, -1]

weylGroupElement[x__Integer][{y__standardBase}]:=Part[{y},{x}]

weylGroupElement[1,2,1,1,2][b2][standardBase[1,1]]

Out[58]= standardBase[1, 1]

we:=weylGroupElement[1,2,1]

Length[we]

Out[56]= 3

we[b2]

Out[55]= Function[z$, Fold[revApply, z$, 
 
>     reflection /@ {standardBase[1, -1], standardBase[0, 1]}[[{1, 2, 1}]]]]

Out[54]= Function[z$, Fold[revApply, z$, 
 
>     reflection /@ {standardBase[1, -1], standardBase[0, 1]}[[{1, 2, 1}]]]]

reflection[b2[[1]]][standardBase[1,1]]

Out[51]= standardBase[1, 1]

Out[50]= {standardBase[1, -1], standardBase[0, 1]}

Out[49]= standardBase[1, 1]

Out[48]= standardBase[1, -1]

Out[47]= standardBase[-1, 1]

Out[45]= standardBase[1, -1][standardBase[0, 1][standardBase[1, -1][
 
>      standardBase[1, 1]]]]

Out[44]= Function[z$, Fold[revApply, z$, 
 
>     {standardBase[1, -1], standardBase[0, 1]}[[{1, 2, 1}]]]]

Out[42]= {standardBase[1, -1], standardBase[0, 1], standardBase[1, -1]}

Part::pspec: Part specification standardBase[0, 1]
     is neither an integer nor a list of integers.

Out[40]= standardBase[1, -1][[standardBase[0, 1],{1, 2, 1}]]

?Part

expr[[i]] or Part[expr, i] gives the i
     th
        part of expr. expr[[-i]]

         counts from the end. expr[[i, j, ...]]

          or Part[expr, i, j, ...]
           is equivalent to expr[[i]][[j]] ...
           . expr[[{i , i , ...}]]
                     1   2
             gives a list of the parts i
                                        1
             , i , ... of expr. expr[[m;;n]]
                2
                  gives parts m through n
                   .expr[[m;;n;;s]] gives parts m through n in steps of s.


Out[37]= Part

Part::pspec: Part specification standardBase[0, 1]
     is neither an integer nor a list of integers.

Out[36]= standardBase[1, -1][[standardBase[0, 1],1,2,1]]

Part::pspec: Part specification standardBase[1, -1]
     is neither an integer nor a list of integers.

Out[34]= standardBase[0, 1][standardBase[1, -1][1[2[1[standardBase[1, 1]]]]]]

Part::pspec: Part specification standardBase[0, 1]
     is neither an integer nor a list of integers.

Out[32]= 1[2[1[standardBase[0, 1][standardBase[1, -1][standardBase[1, 1]]]]]]

Out[30]= weylGroupElement[1, 2, 1][{standardBase[1, -1], standardBase[0, 1]}][
 
>    standardBase[1, 1]]


b2

Out[29]= {standardBase[1, -1], standardBase[0, 1]}


Part[{1,2,3},{1,2,1,1,2,3}]

Out[27]= {1, 2, 1, 1, 2, 3}

Length[weylGroupElement[1,2,3]]



Apply[{reflection[standardBase[1,0]],reflection[standardBase[0,1]]},standardBase[1,1]]



Fold[revApply,standardBase[1,1],{reflection[standardBase[1,0]],reflection[standardBase[0,1]]}]

?Part

Composition[Apply,Reverse][x,f]

?Reverse

Reverse[expr] reverses the order of the elements in expr
    . Reverse[expr, n] reverses elements at level n

       in expr.Reverse[expr, {n , n , ...}]
                               1   2
         reverses elements at levels n , n , ... in expr.
                                      1   2

Reverse::normal: Nonatomic expression expected at position 1 in Reverse[x, f].

Apply::argtu: Apply called with 1 argument; 2 or 3 arguments are expected.

Out[25]= Apply[Reverse[x, f]]

Out[24]= Composition[Apply, Reverse]

?Composition

Composition[f , f , f , ...] represents a composition of the functions f
             1   2   3                                                  1
    , f , f , ... . 
       2   3

System`Compose

Attributes[Compose] = {Protected}

expr[[i]] or Part[expr, i] gives the i
     th
        part of expr. expr[[-i]]

         counts from the end. expr[[i, j, ...]]

          or Part[expr, i, j, ...]
           is equivalent to expr[[i]][[j]] ...
           . expr[[{i , i , ...}]]
                     1   2
             gives a list of the parts i
                                        1
             , i , ... of expr. expr[[m;;n]]
                2
                  gives parts m through n
                   .expr[[m;;n;;s]] gives parts m through n in steps of s.


Out[20]= standardBase[-1, -1]

Out[19]= {standardBase[1, 1], standardBase[-1, 1], standardBase[-1, -1]}

?@

\[FormalA]        \[FormalQ]        \[FormalCapitalG] \[FormalCapitalW]
\[FormalB]        \[FormalR]        \[FormalCapitalH] \[FormalCapitalX]
\[FormalC]        \[FormalS]        \[FormalCapitalI] \[FormalCapitalY]
\[FormalD]        \[FormalT]        \[FormalCapitalJ] \[FormalCapitalZ]
\[FormalE]        \[FormalU]        \[FormalCapitalK] b2
\[FormalF]        \[FormalV]        \[FormalCapitalL] clear
\[FormalG]        \[FormalW]        \[FormalCapitalM] cm
\[FormalH]        \[FormalX]        \[FormalCapitalN] coroot
\[FormalI]        \[FormalY]        \[FormalCapitalO] i
\[FormalJ]        \[FormalZ]        \[FormalCapitalP] j
\[FormalK]        \[FormalCapitalA] \[FormalCapitalQ] rank
\[FormalL]        \[FormalCapitalB] \[FormalCapitalR] reflection
\[FormalM]        \[FormalCapitalC] \[FormalCapitalS] x
\[FormalN]        \[FormalCapitalD] \[FormalCapitalT] y
\[FormalO]        \[FormalCapitalE] \[FormalCapitalU] y$
\[FormalP]        \[FormalCapitalF] \[FormalCapitalV]

b2         clear      cm         coroot     rank       reflection y$

                            2 standardBase[1, 0] . y$ standardBase[1, 0]
Out[15]= {Function[y$, y$ - --------------------------------------------], 
                              standardBase[1, 0] . standardBase[1, 0]
 
                        2 standardBase[0, 1] . y$ standardBase[0, 1]
>     Function[y$, y$ - --------------------------------------------]}[1, 1]
                          standardBase[0, 1] . standardBase[0, 1]

Apply[f, expr] or f @@ expr replaces the head of expr
      by f. Apply[f, expr, levelspec]

        replaces heads in parts of expr specified by levelspec. 

FoldList[f, x, {a, b, ...}] gives {x, f[x, a], f[f[x, a], b], ...}. 

Fold[f, x, list] gives the last element of FoldList[f, x, list]. 

Out[11]= 3

b2=makeSimpleRootSystem[B,2]

Out[2]= {standardBase[1, -1], standardBase[0, 1]}

Out[2]= {standardBase[1, -1], standardBase[0, 1]}

Out[3]= {standardBase[1, -1], standardBase[0, 1]}

Out[64]= {standardBase[1, -1], standardBase[0, 1]}

Out[9]= {standardBase[1, -1], standardBase[0, 1]}

Out[3]= makeSimpleRootSystem[B, 2]

Out[2]= {standardBase[1, -1], standardBase[0, 1]}

cm=cartanMatrix[b2]

Out[8]= {{2, -1}, {-2, 2}}

Inverse[cm].cm

Out[9]= {{1, 0}, {0, 1}}

Out[7]= Inverse[cm] . cm

             1
Out[6]= {{1, -}, {1, 1}}
             2

Out[5]= MatrixInverse[{{2, -1}, {-2, 2}}]

          1          1   1
Out[4]= {{-, -1}, {-(-), -}}
          2          2   2

Out[3]= {{2, -1}, {-2, 2}}

Mathematica 8.0 for Linux x86 (32-bit)
Copyright 1988-2010 Wolfram Research, Inc.

weights[
coroot[b2[[2]]]

coroot/@b2


b2.b2

                           2                      2
Out[61]= standardBase[0, 1]  + standardBase[1, -1]

Out[60]= {standardBase[1, -1], standardBase[0, 1]}

Map[f, expr] or f/@expr applies f
      to each element on the first level in expr
      . Map[f, expr, levelspec] applies f

         to parts of expr specified by levelspec. 


coroot/@makeSimpleRootSystem[D,5]

Out[113]= {standardBase[1, -1, 0, 0, 0, 0], standardBase[0, 1, -1, 0, 0, 0], 
 
>    standardBase[0, 0, 1, -1, 0, 0], standardBase[0, 0, 0, 1, -1, 0], 
 
>    standardBase[0, 0, 0, 0, 1, 1]}

Out[111]= {standardBase[1, -1, 0, 0, 0, 0], standardBase[0, 1, -1, 0, 0, 0], 
 
>    standardBase[0, 0, 1, -1, 0, 0], standardBase[0, 0, 0, 1, -1, 0], 
 
>    standardBase[0, 0, 0, 0, 1, 1]}

cartanMatrix[makeSimpleRootSystem[D,5]]//MatrixForm

Out[125]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   -1

                      0    0    -1   2    0

                      0    0    -1   0    2

???Append

Information::notfound: Symbol ?Append not found.

Append[expr, elem] gives expr with elem appended. 

Attributes[Append] = {Protected}

Append[expr, elem] gives expr with elem appended. 

Out[120]//MatrixForm= 2    -1   0    0

                      -1   2    -1   0

                      0    -1   2    -1

                      0    0    -1   2

Out[118]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -1

                      0    0    0    -1   2

Out[117]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -2

                      0    0    0    -1   2

Out[115]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -1

                      0    0    0    -2   2

Out[114]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -2

                      0    0    0    -1   2

Out[112]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -1

                      0    0    0    -1   2

Out[110]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -1

                      0    0    0    -1   2

Out[109]//MatrixForm= 2    -1   0    0    0

                      -1   2    -1   0    0

                      0    -1   2    -1   0

                      0    0    -1   2    -1

                      0    0    0    -1   2

?//

Postfix[f[expr]] prints with f[expr]
     given in default postfix form: expr//f
     . Postfix[f[expr], h] prints as exprh. 


Out[107]//MatrixForm= 2    -1   0

                      -1   2    -1

                      0    -1   2

Out[106]= {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}

Out[105]= {{2, -2}, {-1, 2}}

Map::level: Level specification standardBase[0, 1]
     is not of the form n, {n}, or {m, n}.

Outer::heads: Heads Map and standardBase at positions 4 and 2
     are expected to be the same.

Out[103]= Outer[Dot, standardBase[1, -1], standardBase[0, 1], 
 
>    Map[coroot, standardBase[1, -1], standardBase[0, 1]]]

Out[101]= cartanMatrix[{standardBase[1, -1], standardBase[0, 1]}]

Outer[f, list , list , ...] gives the generalized outer product of the list
             1      2                                                      i
    , forming all possible combinations of the lowest\[Hyphen]level elements
      in each of them, and feeding them as arguments to f
     . Outer[f, list , list , ..., n]
                    1      2
       treats as separate elements only sublists at level n
        in the list . Outer[f, list , list , ..., n , n , ...]
                   i               1      2        1   2
          treats as separate elements only sublists at level n
                                                              i
           in the corresponding list . 
                                    i

Out[98]= {standardBase[1, -1], standardBase[0, 2]}

Out[97]= standardBase[0, 2]

b2[[1]]

Out[96]= standardBase[1, -1]

Out[95]= standardBase[1, -1]

FullForm[reflection[b2[[2]]][b2[[1]]]]

Out[93]//FullForm= standardBase[1, 1]

Out[90]//FullForm= Plus[Times[2, standardBase[0, 1]], standardBase[1, -1]]

Out[89]= 2 standardBase[0, 1] + standardBase[1, -1]

Out[87]= 2 standardBase[0, 1] + standardBase[1, -1]

Out[85]= 2 standardBase[0, 1] + standardBase[1, -1]

Out[84]= standardBase[1, 0]

Out[81]= standardBase[0, 1] + standardBase[1, -1]

Function::flpar: 
   Parameter specification y_standardBase in 
    Function[y_standardBase, y - 
      2 standardBase[1, -1] . y standardBase[1, -1]
      ---------------------------------------------] should be a symbol or a
        standardBase[1, -1] . standardBase[1, -1]
     list of symbols.

Function::flpar: 
   Parameter specification y_standardBase in 
    Function[y_standardBase, y - 
      2 standardBase[1, -1] . y standardBase[1, -1]
      ---------------------------------------------] should be a symbol or a
        standardBase[1, -1] . standardBase[1, -1]
     list of symbols.

Out[77]= Function[y_standardBase, 
 
          2 standardBase[1, -1] . y standardBase[1, -1]
>     y - ---------------------------------------------][standardBase[0, 1]]
            standardBase[1, -1] . standardBase[1, -1]

?Function

Function[body] or body &
       is a pure function. The formal parameters are # (or #1), #2, etc. Funct

       ion[x, body] is a pure function with a single formal parameter x
      . Function[{x , x , ...}, body]
                   1   2
        is a pure function with a list of formal parameters. 

b2=makeSimpleRootSystem[B,2]

Union[b2,b2]

Out[10]= {standardBase[0, 1], standardBase[1, -1]}

Out[9]= {standardBase[2, -2], standardBase[0, 2]}

Out[8]= {standardBase[1, -1], standardBase[0, 1]}

b2[[1]].b2[[2]]

Out[74]= -1

List@@b2[[1]]

Out[72]= {1, -1}

Out[71]= {standardBase[1, -1]}

Out[70]= standardBase[{1}, {-1}]

Out[69]= {{standardBase[1, -1]}, {standardBase[0, 1]}}

Out[68]= standardBase[List, List]

Out[67]= {standardBase[1, -1], standardBase[0, 1]}

         
SetDelayed::write: 
   Tag Dot in (x_standardBase) . (y_standardBase) is Protected.

Out[64]= $Failed

SetDelayed::write: 
   Tag Dot in (x_standardBase) . (y_standardBase) is Protected.

Out[63]= $Failed        

ExpandNCM[(h : NonCommutativeMultiply)[a___, b_Plus, c___]] :=  Distribute[h[a, b, c], Plus, h, Plus, ExpandNCM[h[##]] &]

ExpandNCM[a_] := ExpandAll[a]

ExpandNCM[(a + b) ** (a + b) ** (a + b)]

keys = DownValues[#,Sort->False][[All,1,1,1]]&;

a[1]=5;
a[3]=7;
a[2]=4

                  
Out[57]= 4

a /@ keys[a]

Out[58]= {5, 7, 4}

Out[56]= {5, 7}

Out[55]= {1, 3}
         
Out[54]= 7


         
Set::write: Tag Times in 5 a[3] is Protected.

Out[53]= 7

LatticeReduce[{{1,-1},{1,0},{1,1},{0,1}}]

h[1,2,3].h[3,4,5]

??Dot

a . b . c or Dot[a, b, c] gives products of vectors, matrices and tensors. 

Attributes[Dot] = {Flat, OneIdentity, Protected}

a . b . c or Dot[a, b, c] gives products of vectors, matrices and tensors. 

Out[49]= h[1, 2, 3] . h[3, 4, 5]



Get::noopen: Cannot open LLLalgorithm.m.

Out[48]= $Failed

Out[47]= {{1, 0}, {0, 1}}

LatticeReduce[{v , v , ...}] gives a reduced basis for the set of vectors v . 
                1   2                                                      i

Out[45]= a ** a ** a + a ** a ** b + a ** b ** a + a ** b ** b + 
 
>    b ** a ** a + b ** a ** b + b ** b ** a + b ** b ** b

{1,2}ExpandNCM[(a - b) ** (a + b)]

Out[44]= {a ** a + a ** b + (-b) ** a + (-b) ** b, 
 
>    2 (a ** a + a ** b + (-b) ** a + (-b) ** b)}

Append[{1,2,3},{4,5}]

Out[33]= {1, 2, 3, {4, 5}}

RuleDelayed::rhs: 
   Pattern rank_Integer appears on the right-hand side of rule 
    makeSimpleRootSystem[A, rank_Integer] :> 
     (standardBase @@ 
         Table[If[i == j, 1, If[i == j - 1, -1, 0]], {i, 1, rank}, 
          {j, 1, rank + 1}] makeSimpleRootSystem[B, rank_Integer] := <<1>>).

         
RuleDelayed::rhs: 
   Pattern rank_Integer appears on the right-hand side of rule 
    makeSimpleRootSystem[A, rank_Integer] :> 
     (standardBase @@ 
         Table[If[i == j, 1, If[i == j - 1, -1, 0]], {i, 1, rank}, 
          {j, 1, rank + 1}] makeSimpleRootSystem[B, rank_Integer] := <<1>>).


tt=makeSimpleRootSystem[D,5]

Out[62]= {standardBase[1, -1, 0, 0, 0, 0], standardBase[0, 1, -1, 0, 0, 0], 
 
>    standardBase[0, 0, 1, -1, 0, 0], standardBase[0, 0, 0, 1, -1, 0], 
 
>    standardBase[0, 0, 0, 0, 1, 1]}

Out[60]= {standardBase[1, -1, 0, 0, 0, 0], standardBase[0, 1, -1, 0, 0, 0], 
 
>    standardBase[0, 0, 1, -1, 0, 0], standardBase[0, 0, 0, 1, -1, 0], 
 
>    standardBase[0, 0, 0, 0, 1, -1]}

Out[41]= standardBase[{1, -1, 0, 0, 0, 0}, {0, 1, -1, 0, 0, 0}, 
 
>    {0, 0, 1, -1, 0, 0}, {0, 0, 0, 1, -1, 0}, {0, 0, 0, 0, 1, 1}]

Module::argrx: Module called with 3 arguments; 2 arguments are expected.

Out[39]= Module[{tmp = makeSimpleRootSystem[A, 5]}, tmp[[-1,-1]] := 1, tmp]

tt[[-1,-1]]

Out[37]= 2

Out[36]= standardBase[{0, 0, 0, 0, 2}, {0, 0, 0, 0, 2}]

Out[35]= standardBase[{1, -1, 0, 0, 0}, {0, 1, -1, 0, 0}, {0, 0, 1, -1, 0}, 
 
>     {0, 0, 0, 1, -1}, {0, 0, 0, 0, 2}][-1]

Out[34]= standardBase[{1, -1, 0, 0, 0}, {0, 1, -1, 0, 0}, {0, 0, 1, -1, 0}, 
 
>    {0, 0, 0, 1, -1}, {0, 0, 0, 0, 2}]

Out[32]= standardBase[{1, -1, 0, 0, 0}, {0, 1, -1, 0, 0}, {0, 0, 1, -1, 0}, 
 
>    {0, 0, 0, 1, -1}, {0, 0, 0, 0, 2}]

Out[30]= standardBase[{1, -1, 0, 0, 0}, {0, 1, -1, 0, 0}, {0, 0, 1, -1, 0}, 
 
>    {0, 0, 0, 1, -1}, {0, 0, 0, 0, 1}]

Out[29]= standardBase[{1, -1}, {0, 1}]

Table::itform: Argument 2 - 1 at position 2
     does not have the correct form for an iterator.

Table::itform: Argument 2 - 1 at position 2
     does not have the correct form for an iterator.

Out[27]= standardBase[{1, -1}, Table[0, 2 - 1, 1]]

Out[22]= makeSimpleRootSystem[B, 2]

Out[20]= standardBase[{1, -1, 0}, {0, 1, -1}]


?@@

Global`rank


t@@{1,2,3}

Out[17]= t[1, 2, 3]

Out[16]= t[{1, 2, 3}]

Out[15]= t

Out[14]= {1, 2, 3}[t]

Out[13]= t

Information::nomatch: No symbol matching //@ found.

\[FormalA]        \[FormalO]        \[FormalCapitalC] \[FormalCapitalQ]
\[FormalB]        \[FormalP]        \[FormalCapitalD] \[FormalCapitalR]
\[FormalC]        \[FormalQ]        \[FormalCapitalE] \[FormalCapitalS]
\[FormalD]        \[FormalR]        \[FormalCapitalF] \[FormalCapitalT]
\[FormalE]        \[FormalS]        \[FormalCapitalG] \[FormalCapitalU]
\[FormalF]        \[FormalT]        \[FormalCapitalH] \[FormalCapitalV]
\[FormalG]        \[FormalU]        \[FormalCapitalI] \[FormalCapitalW]
\[FormalH]        \[FormalV]        \[FormalCapitalJ] \[FormalCapitalX]
\[FormalI]        \[FormalW]        \[FormalCapitalK] \[FormalCapitalY]
\[FormalJ]        \[FormalX]        \[FormalCapitalL] \[FormalCapitalZ]
\[FormalK]        \[FormalY]        \[FormalCapitalM] i
\[FormalL]        \[FormalZ]        \[FormalCapitalN] j
\[FormalM]        \[FormalCapitalA] \[FormalCapitalO] rank
\[FormalN]        \[FormalCapitalB] \[FormalCapitalP]

Information::nomatch: No symbol matching /@ found.

Out[9]= standardBase

?If

If[condition, t, f] gives t if condition
      evaluates to True, and f if it evaluates to False. 

       If[condition, t, f, u] gives u
         if condition evaluates to neither True nor False. 


?Table

Table[expr, {i   }] generates a list of i
              max                        max
     copies of expr. Table[expr, {i, i   }]
                                      max
       generates a list of the values of expr
        when i runs from 1 to i
                               max
          . Table[expr, {i, i   , i   }]
                             min   max
            starts with i = i
                             min
            . Table[expr, {i, i   , i   , di}]
                               min   max
              uses steps di. Table[expr, {i, {i , i , ...}}]
                                               1   2
                uses the successive values i
                                            1
                , i , ....Table[expr, {i, i   , i   }, {j, j   , j   }, ...]
                   2                       min   max        min   max
                    gives a nested list. The list associated with i
                     is outermost. 

?A

Global`A


?Integer

        
Out[4]= A?Integer

Integer is the head used for integers. 

Information::notfound: Symbol PositiveInteger not found.

Mathematica 8.0 for Linux x86 (32-bit)
Copyright 1988-2010 Wolfram Research, Inc.


clear[weight[standardBase]]

makeWeight[coords_standardBase]:=Module[{weight},weight[standardBase]=coords;weight]

w=makeWeight[standardBase[1,2,3]]

w2:=makeWeight[standardBase[1,2,3]]

FullForm[Unevaluated[weight$90]]

Out[58]//FullForm= Unevaluated[weight$90]

FullForm[Symbol["aaa$123"]]



Out[61]= weight?SymbolQ

Out[59]//FullForm= aaa$123

Out[57]//FullForm= aaa

Out[56]= aaa

Symbol::string: String expected at position 1 in Symbol[aaa].

Out[55]= Symbol[aaa]

Out[54]//FullForm= weight$90

Out[53]= Symbol

Out[52]= weight$90[rootSystem]

Out[51]= {HoldPattern[weight$90[standardBase]] :> standardBase[1, 2, 3]}

Out[48]= {HoldPattern[weight$89[rootSystem]] :> {weight$89}, 
 
>    HoldPattern[weight$89[standardBase]] :> standardBase[1, 2, 3]}

Out[47]= {}

Out[46]= {HoldPattern[w] :> weight$89}

Out[45]= {}

Out[42]= weight$89

Out[40]= w$88

w[rootSystem]={w}

Out[43]= {weight$89}

w[rootSystem][[1]][rootSystem][[1]][standardBase]

Out[44]= standardBase[1, 2, 3]


Out[15]= standardBase[1, 2, 3]

Out[14]= weight$82

Out[13]= {weight$82}

w[standardBase]

Out[12]= standardBase[1, 2, 3]

makeWeight[rs__weight,coords_standardBase]:=Module[{weight},weight[rootSystem]=rs;weight[standardBase]=coords;weight]

makeWeight[weight,standardBase[1,2,3]]

Out[3]= makeWeight[rs, standardBase[1, 2, 3]]

Out[1]= clear[weight[standardBase]]

weight[rootSystem]