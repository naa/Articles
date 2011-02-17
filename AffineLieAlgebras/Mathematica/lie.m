ExpandNCM[(h : NonCommutativeMultiply)[a___, b_Plus, c___]] :=  Distribute[h[a, b, c], Plus, h, Plus, ExpandNCM[h[##]] &];
ExpandNCM[a_] := ExpandAll[a];
ExpandNCM[(a + b) ** (a + b) ** (a + b)];
keys = DownValues[#,Sort->False][[All,1,1,1]]&;
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
    NestWhileList[
	Function[x,
		 Union[Flatten[Map[Function[y,
					    Map[reflection[#][y]&,Cases[{simpleRoots},z_ /; z.y>0]]],x]]]],
	{weights},
	#=!={}&];
orbit[{simpleRoots__standardBase}][weight_standardBase]:=orbit[{simpleRoots}][{toFundamentalChamber[{simpleRoots}][weight]}];
positiveRoots[{simpleRoots__standardBase}]:=Map[-#&,Flatten[orbit[{simpleRoots}][Map[-#&,{simpleRoots}]]]];
weightSystem[{simpleRoots__standardBase}][higestWeight_standardBase]:=Module[{minusPosRoots=-positiveRoots[{simpleRoots}]},
									     NestWhileList[Function[x,Complement[
										 Cases[Flatten[Outer[Plus,minusPosRoots,x]],y_/;And@@(#.y>=0&/@{simpleRoots})]
										 ,x]],{higestWeight},#=!={}&]];
freudenthalMultiplicities[{simpleRoots__standardBase}][highestWeight_standardBase]:=
    Module[{rh=rho[{simpleRoots}],weights,mults,c,insideQ,
	    posroots=positiveRoots[{simpleRoots}],
	    toFC=toFundamentalChamber[{simpleRoots}]},
	   weights=SortBy[ Rest[Flatten[weightSystem[{simpleRoots}][highestWeight]]], -#.rh&];
	   c:=(#+rh).(#+rh)&;
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 2/(c[highestWeight]-c[v])*
			 Plus@@
			     Map[Function[r,
					  Plus@@Map[mults[toFC[#[[1]]]]*#[[2]]&,
						    Rest[NestWhileList[({#[[1]]+r,#[[2]]+r.r})&,
								       {v,v.r},
								       insideQ[#[[1]]+r]&]]]]
				 ,posroots]],
		weights];
	   mults];
orbitWithEps[{simpleRoots__standardBase}][weight_standardBase]:=Flatten[Most[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[{simpleRoots}][weight]]],1];
racahMultiplicities[{simpleRoots__standardBase}][highestWeight_standardBase]:=
    Module[{rh=rho[{simpleRoots}],weights,mults,c,insideQ,
	    fan,
	    toFC=toFundamentalChamber[{simpleRoots}]},
	   fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[{simpleRoots}][rh]]];
	   weights=Sort[ Rest[Flatten[weightSystem[{simpleRoots}][highestWeight]]], #1.rh>#2.rh&];
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 Plus@@(fan /. {x_standardBase,e_Integer}:> If[insideQ[v+x],-e*mults[toFC[v+x]],0])],
		weights];
	   mults]


freudenthalMultiplicities[{simpleRoots__standardBase}][highestWeight_standardBase,posroots]:=
    Module[{rh=1/2*Plus@@posroots,weights,mults,c,insideQ,
	    toFC=toFundamentalChamber[{simpleRoots}]},
	   weights=Sort[ Rest[Flatten[weightSystem[{simpleRoots}][highestWeight]]], #1.rh>#2.rh&];
	   c:=(#+rh).(#+rh)&;
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 2/(c[highestWeight]-c[v])*
			 Plus@@
			     Map[Function[r,
					  Plus@@Map[mults[toFC[#[[1]]]]*#[[2]]&,
						    Rest[NestWhileList[({#[[1]]+r,#[[2]]+r.r})&,
								       {v,v.r},
								       insideQ[#[[1]]+r]&]]]]
				 ,posroots]],
		weights];
	   mults];
orbitWithEps[{simpleRoots__standardBase}][weight_standardBase]:=Flatten[Drop[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[{simpleRoots}][weight]],-1],1];
racahMultiplicities[{simpleRoots__standardBase}][highestWeight_standardBase]:=
    Module[{rh=rho[{simpleRoots}],weights,mults,c,insideQ,
	    fan,
	    toFC=toFundamentalChamber[{simpleRoots}]},
	   fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[{simpleRoots}][rh]]];
	   weights=Sort[ Rest[Flatten[weightSystem[{simpleRoots}][highestWeight]]], #1.rh>#2.rh&];
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 Plus@@(fan /. {x_standardBase,e_Integer}:> If[insideQ[v+x],-e*mults[toFC[v+x]],0])],
		weights];
	   mults]


{{standardBase[1,0],2},{standardBase[0,1],2}} /. {x_standardBase,e_Integer} /; x[[1]]>0 :> e*mults[toFC[v+x]]

Out[26]= {2 mults[toFC[v + standardBase[1, 0]]], {standardBase[0, 1], 2}}

Out[25]= {2 mults[toFC[v + standardBase[1, 0]]], 
 
>    2 mults[toFC[v + standardBase[0, 1]]]}

b8=makeSimpleRootSystem[B,8];

wg=fundamentalWeights[b8];

Out[97]= weights[{standardBase[1, -1, 0, 0, 0, 0, 0, 0], 
 
>     standardBase[0, 1, -1, 0, 0, 0, 0, 0], 
 
>     standardBase[0, 0, 1, -1, 0, 0, 0, 0], 
 
>     standardBase[0, 0, 0, 1, -1, 0, 0, 0], 
 
>     standardBase[0, 0, 0, 0, 1, -1, 0, 0], 
 
>     standardBase[0, 0, 0, 0, 0, 1, -1, 0], 
 
>     standardBase[0, 0, 0, 0, 0, 0, 1, -1], 
 
>     standardBase[0, 0, 0, 0, 0, 0, 0, 1]}]

freudenthalMultiplicities[b8][wg[[-1]]*2]//Timing

Out[4]= {0.776049, mults$89}

Out[4][[2]][standardBase[0,0,0,0,0,0,0,1]]

Out[6]= mults$89[standardBase[0, 0, 0, 0, 0, 0, 0, 1]]

Out[5]= 0

racahMultiplicities[b8][wg[[-1]]*2]//Timing



[Calculating...]

Out[4][[2]][standardBase[0, 0, 0, 0, 0, 0, 0, 0]]

Out[7]= 0

Out[5]= 70

Out[4]= {1.22408, mults$89}

mts4=freudenthalMultiplicities[b8][wg[[-1]]*2]

[Calculating...]

wg

Out[103]= mults$266

mts4[standardBase[0, 0, 0, 0, 0, 0, 0, 0]]

Out[104]= 792

Out[102]= 792

Out[101]= {standardBase[5, 1, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[2, 2, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[5, 1, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[4, 2, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[2, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 3, 2, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 2, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[3, 3, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[6, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 2, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[7, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 1, 1, 1, 1, 1, 1], 
 
>    standardBase[2, 2, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 1, 1, 0, 0, 0], 
 
>    standardBase[10, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[1, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[8, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[8, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 2, 2, 0, 0, 0, 0], 
 
>    standardBase[4, 3, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 1, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[5, 2, 2, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 3, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 1, 1, 1, 1, 1, 0], 
 
>    standardBase[2, 1, 1, 1, 1, 1, 1, 0], 
 
>    standardBase[5, 3, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 2, 2, 0, 0, 0, 0], 
 
>    standardBase[1, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 2, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 1, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[4, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[6, 3, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 3, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[2, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 1, 1, 1, 1, 1, 1, 1], 
 
>    standardBase[1, 1, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[8, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 2, 1, 1, 0, 0, 0], 
 
>    standardBase[6, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 3, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 4, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 3, 3, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 2, 1, 1, 0, 0], 
 
>    standardBase[1, 1, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[6, 3, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 2, 2, 1, 0, 0, 0], 
 
>    standardBase[5, 4, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 1, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[3, 3, 3, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 3, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[1, 1, 1, 1, 1, 1, 1, 1], 
 
>    standardBase[3, 2, 2, 1, 1, 0, 0, 0], 
 
>    standardBase[2, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 3, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 2, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[4, 2, 2, 1, 0, 0, 0, 0], 
 
>    standardBase[6, 1, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 2, 1, 0, 0, 0], 
 
>    standardBase[6, 2, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[6, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[6, 2, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 1, 1, 1, 1, 1, 1], 
 
>    standardBase[7, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[5, 5, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[8, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[1, 1, 1, 1, 1, 1, 1, 0], 
 
>    standardBase[5, 3, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 1, 1, 1, 1, 1, 0], 
 
>    standardBase[6, 4, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[7, 3, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 2, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 4, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[7, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[1, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 2, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 2, 2, 0, 0, 0], 
 
>    standardBase[3, 2, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 1, 1, 1, 0, 0, 0], 
 
>    standardBase[6, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 1, 1, 1, 1, 0], 
 
>    standardBase[4, 2, 2, 2, 0, 0, 0, 0], 
 
>    standardBase[5, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 1, 1, 1, 1, 1, 0], 
 
>    standardBase[4, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 2, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[5, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 1, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[5, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[7, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 1, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[4, 2, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[1, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 1, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[7, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 3, 2, 1, 0, 0, 0, 0], 
 
>    standardBase[4, 4, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 2, 2, 1, 1, 0, 0, 0], 
 
>    standardBase[3, 3, 3, 1, 0, 0, 0, 0], 
 
>    standardBase[6, 2, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 1, 1, 1, 1, 1, 1, 0], 
 
>    standardBase[4, 4, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[5, 2, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[9, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 1, 1, 1, 1, 0, 0], 
 
>    standardBase[4, 3, 2, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 1, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 1, 1, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 2, 0, 0, 0, 0], 
 
>    standardBase[2, 2, 2, 1, 1, 1, 0, 0], 
 
>    standardBase[3, 2, 2, 1, 1, 1, 0, 0], 
 
>    standardBase[9, 1, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[0, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[7, 2, 1, 0, 0, 0, 0, 0], 
 
>    standardBase[6, 0, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[3, 2, 0, 0, 0, 0, 0, 0], 
 
>    standardBase[4, 4, 1, 0, 0, 0, 0, 0]}

keys[mts2]

mts2=freudenthalMultiplicities[b2][standardBase[100,0]]

Timing[freudenthalMultiplicities[b2][standardBase[100,0]]]

Out[91]= {74.3446, mults$250}

mts2[standardBase[0,0]]

Out[90]= 51

Out[89]= mults$248

mts2[standardBase[0,0]]

Out[88]= 51

Out[87]= mults$246

Timing[racahMultiplicities[b2][standardBase[100,0]]]

Out[92]= {4.53228, mults$252}

Timing[racahMultiplicities[b2][standardBase[200,0]]]


Out[85]= {18.1491, mults$242}

Timing[freudenthalMultiplicities[b2][standardBase[200,0]]]

Out[86]= $Aborted

mts2[standardBase[0,0]]

Out[84]= 101

Out[83]= mults$240

standardBase[1, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[1, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{0, -1, 0, 0, 0, 0, 0}
1
standardBase[0, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 0, 0, 0, 0, 0}
2

Out[70]= mults$226

standardBase[1, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
-1
standardBase[1, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{0, 1, 0, 0, 0, 0, 0}
1
standardBase[0, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{1, -1, 0, 0, 0, 0, 0}
0

Out[68]= mults$224

mts2[standardBase[1,0]]

Out[66]= 1

Out[65]= 0

standardBase[1, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
standardBase[1, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{0, 1, 0, 0, 0, 0, 0}
standardBase[0, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{1, -1, 0, 0, 0, 0, 0}

Out[64]= mults$222

standardBase[1, 1]
               1  3                      3    1
{{standardBase[-, -], -1}, {standardBase[-, -(-)], -1}, 
               2  2                      2    2
 
                    1   3                     1    3
>   {standardBase[-(-), -], 1}, {standardBase[-, -(-)], 1}, 
                    2   2                     2    2
 
                    3   1                        1     3
>   {standardBase[-(-), -], -1}, {standardBase[-(-), -(-)], -1}, 
                    2   2                        2     2
 
                    3     1
>   {standardBase[-(-), -(-)], 1}}
                    2     2
{0, 0, 0, 0, 0, 0, 0}
standardBase[1, 0]
               1  3                      3    1
{{standardBase[-, -], -1}, {standardBase[-, -(-)], -1}, 
               2  2                      2    2
 
                    1   3                     1    3
>   {standardBase[-(-), -], 1}, {standardBase[-, -(-)], 1}, 
                    2   2                     2    2
 
                    3   1                        1     3
>   {standardBase[-(-), -], -1}, {standardBase[-(-), -(-)], -1}, 
                    2   2                        2     2
 
                    3     1
>   {standardBase[-(-), -(-)], 1}}
                    2     2
{0, 0, 0, 0, 0, 0, 0}
standardBase[0, 0]
               1  3                      3    1
{{standardBase[-, -], -1}, {standardBase[-, -(-)], -1}, 
               2  2                      2    2
 
                    1   3                     1    3
>   {standardBase[-(-), -], 1}, {standardBase[-, -(-)], 1}, 
                    2   2                     2    2
 
                    3   1                        1     3
>   {standardBase[-(-), -], -1}, {standardBase[-(-), -(-)], -1}, 
                    2   2                        2     2
 
                    3     1
>   {standardBase[-(-), -(-)], 1}}
                    2     2
{0, 0, 0, 0, 0, 0, 0}

Out[62]= mults$220

standardBase[1, 1]
                1  3                      3    1
{{{standardBase[-, -], -1}, {standardBase[-, -(-)], -1}}, 
                2  2                      2    2
 
                     1   3                     1    3
>   {{standardBase[-(-), -], 1}, {standardBase[-, -(-)], 1}}, 
                     2   2                     2    2
 
                     3   1                        1     3
>   {{standardBase[-(-), -], -1}, {standardBase[-(-), -(-)], -1}}, 
                     2   2                        2     2
 
                     3     1
>   {{standardBase[-(-), -(-)], 1}}}
                     2     2
{{0, 0}, {0, 0}, {0, 0}, {0}}

Thread::tdlen: Objects of unequal length in {0, 0} + {0, 0} + {0, 0} + {0}
     cannot be combined.

Thread::tdlen: Objects of unequal length in {0} + {0, 0} cannot be combined.
standardBase[1, 0]
                1  3                      3    1
{{{standardBase[-, -], -1}, {standardBase[-, -(-)], -1}}, 
                2  2                      2    2
 
                     1   3                     1    3
>   {{standardBase[-(-), -], 1}, {standardBase[-, -(-)], 1}}, 
                     2   2                     2    2
 
                     3   1                        1     3
>   {{standardBase[-(-), -], -1}, {standardBase[-(-), -(-)], -1}}, 
                     2   2                        2     2
 
                     3     1
>   {{standardBase[-(-), -(-)], 1}}}
                     2     2
{{0, 0}, {0, 0}, {0, 0}, {0}}

Thread::tdlen: Objects of unequal length in {0, 0} + {0, 0} + {0, 0} + {0}
     cannot be combined.

General::stop: Further output of Thread::tdlen
     will be suppressed during this calculation.
standardBase[0, 0]
                1  3                      3    1
{{{standardBase[-, -], -1}, {standardBase[-, -(-)], -1}}, 
                2  2                      2    2
 
                     1   3                     1    3
>   {{standardBase[-(-), -], 1}, {standardBase[-, -(-)], 1}}, 
                     2   2                     2    2
 
                     3   1                        1     3
>   {{standardBase[-(-), -], -1}, {standardBase[-(-), -(-)], -1}}, 
                     2   2                        2     2
 
                     3     1
>   {{standardBase[-(-), -(-)], 1}}}
                     2     2
{{0, 0}, {0, 0}, {0, 0}, {0}}

Out[60]= mults$218

standardBase[1, 1]
{{0, 0}, {0, 0}, {0, 0}, {0}}

Thread::tdlen: Objects of unequal length in {0, 0} + {0, 0} + {0, 0} + {0}
     cannot be combined.

Thread::tdlen: Objects of unequal length in {0} + {0, 0} cannot be combined.
standardBase[1, 0]
{{0, 0}, {0, 0}, {0, 0}, {0}}

Thread::tdlen: Objects of unequal length in {0, 0} + {0, 0} + {0, 0} + {0}
     cannot be combined.

General::stop: Further output of Thread::tdlen
     will be suppressed during this calculation.
standardBase[0, 0]
{{0, 0}, {0, 0}, {0, 0}, {0}}

Out[58]= mults$216

Thread::tdlen: Objects of unequal length in {0, 0} + {0, 0} + {0, 0} + {0}
     cannot be combined.

Thread::tdlen: Objects of unequal length in {0} + {0, 0} cannot be combined.

Thread::tdlen: Objects of unequal length in {0, 0} + {0, 0} + {0, 0} + {0}
     cannot be combined.

General::stop: Further output of Thread::tdlen
     will be suppressed during this calculation.

Out[52]= mults$210

mts2[standardBase[1,0]]

Out[50]= 0

Out[49]= 0

             3  5
standardBase[-, -]
             2  2
             5  1
standardBase[-, -]
             2  2
             1  5
standardBase[-, -]
             2  2
             3    1
standardBase[-, -(-)]
             2    2
               1   3
standardBase[-(-), -]
               2   2
             1    1
standardBase[-, -(-)]
             2    2
               1   1
standardBase[-(-), -]
               2   2
             3  3
standardBase[-, -]
             2  2
             5    1
standardBase[-, -(-)]
             2    2
             1  3
standardBase[-, -]
             2  2
             3    3
standardBase[-, -(-)]
             2    2
               1   1
standardBase[-(-), -]
               2   2
             1    3
standardBase[-, -(-)]
             2    2
               1     1
standardBase[-(-), -(-)]
               2     2
             1  3
standardBase[-, -]
             2  2
             3    1
standardBase[-, -(-)]
             2    2
               1   3
standardBase[-(-), -]
               2   2
             1    3
standardBase[-, -(-)]
             2    2
               3   1
standardBase[-(-), -]
               2   2
               1     3
standardBase[-(-), -(-)]
               2     2
               3     1
standardBase[-(-), -(-)]
               2     2

Out[48]= mults$208

mts2[1,0]

Scan::argtu: Scan called with 1 argument; 2 or 3 arguments are expected.

Out[17]= mults$99[1, 0]

MapIndexed[g,{{1,2},{3,4}},1]

Out[5]= {g[{1, 2}, {1}], g[{3, 4}, {2}]}

Out[4]= {g[{g[1, {1, 1}], g[2, {1, 2}]}, {1}], 
 
>    g[{g[3, {2, 1}], g[4, {2, 2}]}, {2}]}

Out[3]= {g[{g[1, {1, 1}], g[2, {1, 2}]}, {1}], 
 
>    g[{g[3, {2, 1}], g[4, {2, 2}]}, {2}]}

Out[2]= {{{{1, {1, 1}}, {2, {1, 2}}}, {1}}, {{{3, {2, 1}}, {4, {2, 2}}}, {2}}}

Out[1]= {{{1, 2}, {1}}, {{3, 4}, {2}}}

orbit[b2][rho[b2]]

                        3  1                  1  3                3    1
Out[36]= {{standardBase[-, -]}, {standardBase[-, -], standardBase[-, -(-)]}, 
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

                        3  1
Out[34]= {{standardBase[-, -]}}
                        2  2

b2

Out[33]= {standardBase[1, -1], standardBase[0, 1]}

                        3  1
Out[32]= {{standardBase[-, -]}}
                        2  2

                           
Syntax::sntxf: "b2 \                          " cannot be followed by 
                                        3  1 \
                Out[32]= {{standardBase[-
    ", -]}} \".

                        3  1                  1  3                3    1
Out[76]= {{standardBase[-, -]}, {standardBase[-, -], standardBase[-, -(-)]}, 
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

freudenthal[{simpleRoots__standardBase}][highestWeight_standardBase]:=
    

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

rho[b2]

                      3  1
Out[47]= standardBase[-, -]
                      2  2

Sort[ Rest[Flatten[weightSystem[b2][standardBase[2,0]]]], #1.rho[b2]>#2.rho[b2]&]

Out[48]= {standardBase[1, 1], standardBase[1, 0], standardBase[0, 0]}

mts=freudenthalMultiplicities[b2][standardBase[44,0]]

Timing[racahMultiplicities[b2][standardBase[54,0]]]


mts3:=racahMultiplicities[b2][standardBase[54,0]]

mts3[standardBase[1,0]]

Out[77]= 27

Out[76]= 28

Out[73]= {1.41609, mults$230}

standardBase[53, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[53, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{0, -1, 0, 0, 0, 0, 0}
1
standardBase[52, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[52, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[51, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[52, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 0, 0, 0, 0, 0}
2
standardBase[51, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[50, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[51, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[50, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[51, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -2, 1, 1, -1, 0, 0}
2
standardBase[49, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[50, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[49, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[50, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[48, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[49, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[50, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 1, 1, -1, -1, 1}
3
standardBase[48, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[49, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[47, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[48, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[49, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[47, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[48, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[49, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -3, 2, 2, -2, -1, 1}
3
standardBase[46, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[47, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[48, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[46, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[47, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[48, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[45, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[46, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[47, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[48, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 2, 2, -2, -2, 2}
4
standardBase[45, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[46, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[47, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[44, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[45, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[46, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[47, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[44, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[45, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[46, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[47, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -4, 3, 3, -3, -2, 2}
4
standardBase[43, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[44, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[45, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[46, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[43, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[44, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[45, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[46, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[42, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[43, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[44, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[45, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[46, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 3, 3, -3, -3, 3}
5
standardBase[42, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[43, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[44, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[45, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[41, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[42, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[43, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[44, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[45, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[41, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[42, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[43, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[44, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[45, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -5, 4, 4, -4, -3, 3}
5
standardBase[40, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[41, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[42, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[43, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[44, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[40, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[41, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[42, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[43, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[44, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[39, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[40, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[41, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[42, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[43, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[44, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 4, 4, -4, -4, 4}
6
standardBase[39, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[40, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[41, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[42, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[43, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[38, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[39, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[40, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[41, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[42, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[43, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[38, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[39, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[40, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[41, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[42, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[43, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -6, 5, 5, -5, -4, 4}
6
standardBase[37, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[38, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[39, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[40, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[41, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[42, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[37, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[38, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[39, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[40, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[41, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[42, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[36, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[37, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[38, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[39, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[40, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[41, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[42, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 5, 5, -5, -5, 5}
7
standardBase[36, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[37, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[38, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[39, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[40, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[41, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[35, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[36, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[37, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[38, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[39, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[40, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[41, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[35, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[36, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[37, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[38, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[39, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[40, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[41, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -7, 6, 6, -6, -5, 5}
7
standardBase[34, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[35, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[36, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[37, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[38, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[39, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[40, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[34, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[35, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[36, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[37, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[38, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[39, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[40, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[33, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[34, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[35, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[36, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[37, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[38, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[39, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[40, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 6, 6, -6, -6, 6}
8
standardBase[33, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[34, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[35, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[36, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[37, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[38, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[39, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[32, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[33, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[34, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[35, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[36, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[37, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[38, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[39, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[32, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[33, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[34, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[35, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[36, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[37, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[38, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[39, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -8, 7, 7, -7, -6, 6}
8
standardBase[31, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[32, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[33, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[34, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[35, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[36, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[37, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[38, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[31, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[32, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[33, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[34, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[35, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[36, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[37, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[38, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[30, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[31, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[32, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[33, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[34, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[35, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[36, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[37, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[38, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 7, 7, -7, -7, 7}
9
standardBase[30, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[31, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[32, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[33, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[34, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[35, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[36, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[37, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[29, 25]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[30, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[31, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[32, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[33, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[34, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[35, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[36, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[37, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[29, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[30, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[31, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[32, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[33, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[34, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[35, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[36, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[37, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -9, 8, 8, -8, -7, 7}
9
standardBase[28, 26]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[29, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[30, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[31, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[32, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[33, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[34, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[35, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[36, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[28, 25]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[29, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[30, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[31, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[32, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[33, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[34, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[35, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[36, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[27, 27]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, 0, 0, 0, 0, 0, 0}
1
standardBase[28, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[29, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[30, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[31, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[32, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[33, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[34, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[35, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[36, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 8, 8, -8, -8, 8}
10
standardBase[27, 26]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-1, -1, 1, 0, 0, 0, 0}
1
standardBase[28, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[29, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[30, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[31, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[32, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[33, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[34, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[35, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[27, 25]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[28, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[29, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[30, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[31, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[32, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[33, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[34, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[35, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[27, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[28, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[29, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[30, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[31, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[32, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[33, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[34, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[35, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -10, 9, 9, -9, -8, 8}
10
standardBase[26, 26]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -1, 1, 0, 0, 0, 0}
2
standardBase[27, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[28, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[29, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[30, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[31, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[32, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[33, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[34, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[26, 25]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-2, -2, 2, 1, -1, 0, 0}
2
standardBase[27, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[28, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[29, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[30, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[31, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[32, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[33, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[34, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[26, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[27, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[28, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[29, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[30, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[31, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[32, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[33, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[34, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 9, 9, -9, -9, 9}
11
standardBase[26, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[27, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[28, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[29, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[30, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[31, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[32, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[33, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[25, 25]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -2, 2, 1, -1, -1, 1}
3
standardBase[26, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[27, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[28, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[29, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[30, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[31, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[32, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[33, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[25, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-3, -3, 3, 2, -2, -1, 1}
3
standardBase[26, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[27, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[28, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[29, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[30, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[31, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[32, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[33, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -11, 10, 10, -10, -9, 9}
11
standardBase[25, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[26, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[27, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[28, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[29, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[30, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[31, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[32, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[25, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[26, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[27, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[28, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[29, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[30, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[31, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[32, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[24, 24]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -3, 3, 2, -2, -2, 2}
4
standardBase[25, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[26, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[27, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[28, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[29, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[30, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[31, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[32, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 10, 10, -10, -10, 10}
12
standardBase[24, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-4, -4, 4, 3, -3, -2, 2}
4
standardBase[25, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[26, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[27, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[28, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[29, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[30, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[31, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[24, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[25, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[26, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[27, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[28, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[29, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[30, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[31, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[24, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[25, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[26, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[27, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[28, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[29, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[30, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[31, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -12, 11, 11, -11, -10, 10}
12
standardBase[23, 23]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -4, 4, 3, -3, -3, 3}
5
standardBase[24, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[25, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[26, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[27, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[28, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[29, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[30, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[23, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-5, -5, 5, 4, -4, -3, 3}
5
standardBase[24, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[25, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[26, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[27, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[28, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[29, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[30, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[23, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[24, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[25, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[26, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[27, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[28, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[29, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[30, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 11, 11, -11, -11, 11}
13
standardBase[23, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[24, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[25, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[26, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[27, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[28, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[29, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[22, 22]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -5, 5, 4, -4, -4, 4}
6
standardBase[23, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[24, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[25, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[26, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[27, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[28, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[29, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[22, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-6, -6, 6, 5, -5, -4, 4}
6
standardBase[23, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[24, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[25, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[26, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[27, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[28, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[29, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -13, 12, 12, -12, -11, 11}
13
standardBase[22, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[23, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[24, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[25, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[26, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[27, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[28, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[22, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[23, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[24, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[25, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[26, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[27, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[28, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[21, 21]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -6, 6, 5, -5, -5, 5}
7
standardBase[22, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[23, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[24, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[25, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[26, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[27, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[28, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 12, 12, -12, -12, 12}
14
standardBase[21, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-7, -7, 7, 6, -6, -5, 5}
7
standardBase[22, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[23, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[24, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[25, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[26, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[27, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[21, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[22, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[23, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[24, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[25, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[26, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[27, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[21, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[22, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[23, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[24, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[25, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[26, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[27, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -14, 13, 13, -13, -12, 12}
14
standardBase[20, 20]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -7, 7, 6, -6, -6, 6}
8
standardBase[21, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[22, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[23, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[24, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[25, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[26, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[20, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-8, -8, 8, 7, -7, -6, 6}
8
standardBase[21, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[22, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[23, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[24, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[25, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[26, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[20, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[21, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[22, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[23, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[24, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[25, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[26, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 13, 13, -13, -13, 13}
15
standardBase[20, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[21, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[22, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[23, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[24, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[25, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[19, 19]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -8, 8, 7, -7, -7, 7}
9
standardBase[20, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[21, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[22, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[23, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[24, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[25, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[19, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-9, -9, 9, 8, -8, -7, 7}
9
standardBase[20, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[21, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[22, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[23, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[24, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[25, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -15, 14, 14, -14, -13, 13}
15
standardBase[19, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[20, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[21, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[22, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[23, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[24, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[19, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[20, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[21, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[22, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[23, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[24, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[18, 18]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -9, 9, 8, -8, -8, 8}
10
standardBase[19, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[20, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[21, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[22, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[23, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[24, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 14, 14, -14, -14, 14}
16
standardBase[18, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-10, -10, 10, 9, -9, -8, 8}
10
standardBase[19, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[20, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[21, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[22, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[23, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[18, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[19, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[20, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[21, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[22, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[23, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[18, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[19, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[20, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[21, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[22, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[23, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -16, 15, 15, -15, -14, 14}
16
standardBase[17, 17]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -10, 10, 9, -9, -9, 9}
11
standardBase[18, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[19, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[20, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[21, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[22, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[17, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-11, -11, 11, 10, -10, -9, 9}
11
standardBase[18, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[19, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[20, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[21, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[22, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[17, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[18, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[19, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[20, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[21, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[22, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 15, 15, -15, -15, 15}
17
standardBase[17, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[18, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[19, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[20, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[21, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[16, 16]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -11, 11, 10, -10, -10, 10}
12
standardBase[17, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[18, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[19, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[20, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[21, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[16, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-12, -12, 12, 11, -11, -10, 10}
12
standardBase[17, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[18, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[19, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[20, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[21, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -17, 16, 16, -16, -15, 15}
17
standardBase[16, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[17, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[18, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[19, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[20, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[16, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[17, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[18, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[19, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[20, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[15, 15]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -12, 12, 11, -11, -11, 11}
13
standardBase[16, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[17, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[18, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[19, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[20, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 16, 16, -16, -16, 16}
18
standardBase[15, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-13, -13, 13, 12, -12, -11, 11}
13
standardBase[16, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[17, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[18, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[19, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[15, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[16, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[17, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[18, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[19, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[15, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[16, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[17, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[18, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[19, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -18, 17, 17, -17, -16, 16}
18
standardBase[14, 14]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -13, 13, 12, -12, -12, 12}
14
standardBase[15, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[16, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[17, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[18, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[14, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-14, -14, 14, 13, -13, -12, 12}
14
standardBase[15, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[16, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[17, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[18, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[14, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[15, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[16, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[17, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[18, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 17, 17, -17, -17, 17}
19
standardBase[14, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[15, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[16, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[17, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[13, 13]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -14, 14, 13, -13, -13, 13}
15
standardBase[14, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[15, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[16, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[17, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[13, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-15, -15, 15, 14, -14, -13, 13}
15
standardBase[14, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[15, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[16, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[17, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -19, 18, 18, -18, -17, 17}
19
standardBase[13, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[14, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[15, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[16, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[13, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[14, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[15, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[16, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[12, 12]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -15, 15, 14, -14, -14, 14}
16
standardBase[13, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[14, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[15, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[16, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 18, 18, -18, -18, 18}
20
standardBase[12, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-16, -16, 16, 15, -15, -14, 14}
16
standardBase[13, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[14, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[15, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[12, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[13, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[14, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[15, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[12, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[13, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[14, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[15, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -20, 19, 19, -19, -18, 18}
20
standardBase[11, 11]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -16, 16, 15, -15, -15, 15}
17
standardBase[12, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[13, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[14, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[11, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-17, -17, 17, 16, -16, -15, 15}
17
standardBase[12, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[13, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[14, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[11, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[12, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[13, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[14, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 19, 19, -19, -19, 19}
21
standardBase[11, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[12, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[13, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[10, 10]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -17, 17, 16, -16, -16, 16}
18
standardBase[11, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[12, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[13, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[10, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-18, -18, 18, 17, -17, -16, 16}
18
standardBase[11, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[12, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[13, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -21, 20, 20, -20, -19, 19}
21
standardBase[10, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[11, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[12, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[10, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[11, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[12, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 21, 20, -20, -19, 19}
21
standardBase[9, 9]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -18, 18, 17, -17, -17, 17}
19
standardBase[10, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[11, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[12, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 20, 20, -20, -20, 20}
22
standardBase[9, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-19, -19, 19, 18, -18, -17, 17}
19
standardBase[10, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[11, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 21, 20, -20, -19, 19}
21
standardBase[9, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[10, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[11, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -21, 21, 20, -20, -20, 20}
22
standardBase[9, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[10, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 21, 20, -20, -19, 19}
21
standardBase[11, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -22, 21, 21, -21, -20, 20}
22
standardBase[8, 8]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -19, 19, 18, -18, -18, 18}
20
standardBase[9, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[10, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -21, 21, 20, -20, -20, 20}
22
standardBase[8, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-20, -20, 20, 19, -19, -18, 18}
20
standardBase[9, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 21, 20, -20, -19, 19}
21
standardBase[10, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -22, 22, 21, -21, -20, 20}
22
standardBase[8, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[9, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -21, 21, 20, -20, -20, 20}
22
standardBase[10, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -22, 21, 21, -21, -21, 21}
23
standardBase[8, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 21, 20, -20, -19, 19}
21
standardBase[9, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -22, 22, 21, -21, -20, 20}
22
standardBase[7, 7]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -20, 20, 19, -19, -19, 19}
21
standardBase[8, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -21, 21, 20, -20, -20, 20}
22
standardBase[9, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -22, 22, 21, -21, -21, 21}
23
standardBase[7, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-21, -21, 21, 20, -20, -19, 19}
21
standardBase[8, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -22, 22, 21, -21, -20, 20}
22
standardBase[9, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -23, 22, 22, -22, -21, 21}
23
standardBase[7, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -21, 21, 20, -20, -20, 20}
22
standardBase[8, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -22, 22, 21, -21, -21, 21}
23
standardBase[7, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -22, 22, 21, -21, -20, 20}
22
standardBase[8, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -23, 23, 22, -22, -21, 21}
23
standardBase[6, 6]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -21, 21, 20, -20, -20, 20}
22
standardBase[7, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -22, 22, 21, -21, -21, 21}
23
standardBase[8, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -23, 22, 22, -22, -22, 22}
24
standardBase[6, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-22, -22, 22, 21, -21, -20, 20}
22
standardBase[7, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -23, 23, 22, -22, -21, 21}
23
standardBase[6, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -22, 22, 21, -21, -21, 21}
23
standardBase[7, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -23, 23, 22, -22, -22, 22}
24
standardBase[6, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -23, 23, 22, -22, -21, 21}
23
standardBase[7, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -24, 23, 23, -23, -22, 22}
24
standardBase[5, 5]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -22, 22, 21, -21, -21, 21}
23
standardBase[6, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -23, 23, 22, -22, -22, 22}
24
standardBase[5, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-23, -23, 23, 22, -22, -21, 21}
23
standardBase[6, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -24, 24, 23, -23, -22, 22}
24
standardBase[5, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -23, 23, 22, -22, -22, 22}
24
standardBase[6, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -24, 23, 23, -23, -23, 23}
25
standardBase[5, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -24, 24, 23, -23, -22, 22}
24
standardBase[4, 4]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -23, 23, 22, -22, -22, 22}
24
standardBase[5, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -24, 24, 23, -23, -23, 23}
25
standardBase[4, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -24, 24, 23, -23, -22, 22}
24
standardBase[5, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-24, -25, 24, 24, -24, -23, 23}
25
standardBase[4, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -24, 24, 23, -23, -23, 23}
25
standardBase[4, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -25, 25, 24, -24, -23, 23}
25
standardBase[3, 3]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -24, 24, 23, -23, -23, 23}
25
standardBase[4, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -25, 24, 24, -24, -24, 24}
26
standardBase[3, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -25, 25, 24, -24, -23, 23}
25
standardBase[3, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-26, -25, 25, 24, -24, -24, 24}
26
standardBase[3, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-25, -26, 25, 25, -25, -24, 24}
26
standardBase[2, 2]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-26, -25, 25, 24, -24, -24, 24}
26
standardBase[2, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-26, -26, 26, 25, -25, -24, 24}
26
standardBase[2, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-26, -26, 25, 25, -25, -25, 25}
27
standardBase[1, 1]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-27, -26, 26, 25, -25, -25, 25}
27
standardBase[1, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-26, -27, 26, 26, -26, -25, 25}
27
standardBase[0, 0]
{{standardBase[1, -1], -1}, {standardBase[0, 1], -1}, 
 
>   {standardBase[2, -1], 1}, {standardBase[1, 2], 1}, 
 
>   {standardBase[3, 0], -1}, {standardBase[2, 2], -1}, 
 
>   {standardBase[3, 1], 1}}
{-27, -27, 26, 26, -26, -26, 26}
28

Out[71]= {2.68417, mults$228}

Out[75]= {12.5328, mults$564}

mts[standardBase[0,0]]

Out[72]= 3

1
1
1
1
2
2
2
3

                                                                        
Out[71]= 72

Out[70]= mults$544

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

Out[9]= {standardBase[1, -1], standardBase[0, 1]}

Out[7]= makeSimpleRootSystem[B, 2]

Out[6]= {standardBase[1, -1], standardBase[0, 1]}

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


a[0]=1;
a[1]=2

keys[a]

Out[5]= {0, 1}

        
Out[3]= 2

Out[2]= 1

Set::noval: Symbol a in part assignment does not have an immediate value.

Out[1]= 1

hashtable:=Module[{table},table[5]=1;table]

t1=hashtable

Out[16]= table$97

t1[5]

t2=hashtable

t2[5]=7

t1[5]

Out[22]= 1

Out[21]= 7

Out[20]= 7

Out[19]= 1

Out[18]= table$98

Out[17]= 1

Out[15]= table$96[][5]

Out[14]= table$96[]

Out[13]= 1

Module::argrx: Module called with 3 arguments; 2 arguments are expected.

Out[11]= Module[{table}, table[5] = 1, table][5]

Out[9]= table$92[5]

Out[8]= table$91[]

b2=makeSimpleRootSystem[B,2]

makeSimpleAffineRootSystem[A,rank_Integer]:=Module[{roots=makeSimpleRootSystem[A,rank],affineRoots,affineRoot}, 
						   affineRoot[grade]=1;
						   affineRoot[finite]=-Plus@@roots;
						   affineRoots[0]=affineRoot;
						   MapIndexed[affine

b1[1]=1;b1[2]=2;b1[3]=3

DownValues[b1]

Out[31]= {HoldPattern[b1[1]] :> 1, HoldPattern[b1[2]] :> 2, 
 
>    HoldPattern[b1[3]] :> 3}

Out[30]= b1

Out[29]= 3

b2[[1]]

Out[26]= standardBase[1, -1]

Out[25]= List

Out[24]= {standardBase[1, -1], standardBase[0, 1]}[0]

Out[23]= {standardBase[1, -1], standardBase[0, 1]}