(* 

   finiteWeight[dimension]=dimension_Integer
   finiteWeight[standardBase]=List[x__Integer]
   
   affineWeight[finitePart]=x_finiteWeight
   affineWeight[grade]=grade_Integer
   affineWeight[level]=level_Integer
   

   Finite Weight: finiteWeight[dimension][<id>][dimension]
                                              [standardBase]
   Affine Weight: affineWeight[dimension][<id>][finitePart]
   Affine Weight: affineWeight[dimension][<id>][level]
   Affine Weight: affineWeight[dimension][<id>][grade]
   
   *)

Clear[finiteWeight,affineWeight,dimension,standardBase,simpleRootBase,level,grade,Expect];

Expect[ description_, val_, expr_ ] := 
If[
    val != expr,
    Throw[
        StringJoin[ "GOT UNEXPECTED VALUE ", ToString[expr], 
        " INSTEAD OF ", ToString[val] ]
        , "assertion exception"
    ]
]

(* finiteWeight[dimension,coordinates] *)

finiteWeight/:x_finiteWeight[dimension]:=x[[1]];
finiteWeight/:x_finiteWeight[standardBase]:=x[[2]];

makeFiniteWeight[{coordinates__?NumberQ}]:=finiteWeight @@ {Length[{coordinates}],{coordinates}}

Expect["Dimension equals to length",3,makeFiniteWeight[{1,2,3}][dimension]]

finiteWeight/:x_finiteWeight . y_finiteWeight/;x[dimension]==y[dimension]:=x[standardBase].y[standardBase]

Expect["Scalar product for vectors from weigth space of finite-dimensional Lie algebras",10,makeFiniteWeight[{1,2,3}].makeFiniteWeight[{3,2,1}]]


Expect["Scalar product for vectors from different spaces are left unevaluated",True,
       MatchQ[makeFiniteWeight[{1,2,3}].makeFiniteWeight[{3,2,1,2}],x_finiteWeight . y_finiteWeight]]

finiteWeight/:x_finiteWeight+y_finiteWeight/;x[dimension]==y[dimension]:=makeFiniteWeight[x[standardBase]+y[standardBase]]

Expect["Plus for finite-dimensional weights",{2,4,6},(makeFiniteWeight[{1,2,3}]+makeFiniteWeight[{1,2,3}])[standardBase]]

Expect["Plus product for vectors from different spaces are left unevaluated",True,
       MatchQ[makeFiniteWeight[{1,2,3}]+makeFiniteWeight[{3,2,1,2}],x_finiteWeight + y_finiteWeight]]

finiteWeight/:x_finiteWeight==y_finiteWeight:=x[standardBase]==y[standardBase]

Expect["Equal for finite weights compares standard base representations", False,makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{1,3,2}]]

Expect["Equal for finite weights compares standard base representations", True,makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{1,2,3}]]

Expect["Equal for finite weights compares standard base representations", False,makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{1,2,3,4}]]

finiteWeight/:x_?NumberQ*y_finiteWeight:=makeFiniteWeight[x*y[standardBase]]

Expect["Multiplication by scalar", True,makeFiniteWeight[{1,2,3}]*2==makeFiniteWeight[{2,4,6}]]

Expect["Multiplication by scalar", True,2*makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{2,4,6}]]

makeAffineWeight[fw_finiteWeight,lev_?NumberQ,gr_?NumberQ]:=affineWeight[fw[dimension],fw,lev,gr]

affineWeight/:x_affineWeight[dimension]:=x[[1]];
affineWeight/:x_affineWeight[finitePart]:=x[[2]];
affineWeight/:x_affineWeight[level]:=x[[3]];
affineWeight/:x_affineWeight[grade]:=x[[4]];

Expect["Affine weight has the same real dimension as the finite-dimensional part, since we hold level and grade separetely",
       True,makeAffineWeight[makeFiniteWeight[{1,2,3,4,5}],1,2][dimension]==Length[{1,2,3,4,5}]]

makeAffineWeight[{fp__?NumberQ},lev_?NumberQ,gr_?NumberQ]:=makeAffineWeight[makeFiniteWeight[{fp}],lev,gr]

Expect["Shortened constructor",
       True,makeAffineWeight[{1,2,3,4,5},1,2][dimension]==Length[{1,2,3,4,5}]]

affineWeight/:x_affineWeight==y_affineWeight:=And[x[finitePart]==y[finitePart],
						      x[level]==y[level],
						      x[grade]==y[grade]]

Expect["Equal for affine weights compares finite parts, levels and grades", True,makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]==makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]]

Expect["Equal for affine weights compares finite parts, levels and grades", False,makeAffineWeight[makeFiniteWeight[{1,2,3}],2,1]==makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]]

affineWeight/:x_affineWeight+y_affineWeight/;x[dimension]==y[dimension]:=
    makeAffineWeight[x[finitePart]+y[finitePart],
		     x[level]+y[level],
		     x[grade]+y[grade]]


Expect["Plus for affine weights",{2,4,6},(makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]+ makeAffineWeight[makeFiniteWeight[{1,2,3}],3,1])[finitePart][standardBase]]

Expect["Plus for affine weights",4,(makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]+ makeAffineWeight[makeFiniteWeight[{1,2,3}],3,1])[level]]

Expect["We compare dimensions of vectors before sum calculation, the expression is left unevaluated in case of dimension mismatch ",
       True,MatchQ[makeAffineWeight[makeFiniteWeight[{1,2}],1,2] + makeAffineWeight[ makeFiniteWeight[{3,2,1}],2,1], x_affineWeight + y_affineWeight]]

affineWeight/:x_affineWeight.y_affineWeight/;x[dimension]==y[dimension]:= 
    x[finitePart].y[finitePart] + 
    x[level]* y[grade] + 
    x[grade]* y[level]

Expect["Scalar product for vectors from weigth space of affine Lie algebras",20,
       makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]. 
       makeAffineWeight[makeFiniteWeight[{3,2,1}],3,4]]

Expect["We compare dimensions of vectors before product calculation, the expression is left unevaluated in case of dimension mismatch ",
       True,MatchQ[makeAffineWeight[makeFiniteWeight[{1,2}],1,2]. makeAffineWeight[ makeFiniteWeight[{3,2,1}],2,1], x_affineWeight . y_affineWeight]]

affineWeight/:x_?NumberQ*y_affineWeight:=makeAffineWeight[x*y[finitePart],x*y[level],x*y[grade]]

Expect["Multiplication by scalar", True,makeAffineWeight[makeFiniteWeight[{1,2,3}],1,2]*2==makeAffineWeight[makeFiniteWeight[{2,4,6}],2,4]]


ExpandNCM[(h : NonCommutativeMultiply)[a___, b_Plus, c___]] :=  Distribute[h[a, b, c], Plus, h, Plus, ExpandNCM[h[##]] &];
ExpandNCM[a_] := ExpandAll[a];
ExpandNCM[(a + b) ** (a + b) ** (a + b)];
keys = DownValues[#,Sort->False][[All,1,1,1]]&;

makeFiniteRootSystem[{roots__finiteWeight}]:=finiteRootSystem[Length[{roots}],{roots}]

finiteRootSystem/:x_finiteRootSystem[rank]:=x[[1]];
finiteRootSystem/:x_finiteRootSystem[simpleRoots]:=x[[2]];
finiteRootSystem/:x_finiteRootSystem[simpleRoot][n_Integer]:=x[[2]][[n]];

makeFiniteRootSystem[{roots__List}]:=makeFiniteRootSystem[makeFiniteWeight/@{roots}]

makeSimpleRootSystem[A,1]:=makeFiniteRootSystem[{{1}}];
makeSimpleRootSystem[A,r_Integer]:=makeFiniteRootSystem[makeFiniteWeight /@ Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,r},{j,1,r+1}]];
makeSimpleRootSystem[B,rank_Integer]:=makeFiniteRootSystem[Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],1]]];
makeSimpleRootSystem[C,rank_Integer]:=makeFiniteRootSystem[Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],2]]];
makeSimpleRootSystem[D,rank_Integer]:=makeFiniteRootSystem[Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Append[Table[0,{rank-2}],1],1]]];

Expect["B2:",True,makeSimpleRootSystem[B,2][simpleRoots][[1]]==makeFiniteWeight[{1,-1}]]

Expect["B2: rank",2,makeSimpleRootSystem[B,2][rank]]

reflection[x_finiteWeight]:=Function[y, y-2*(x.y)/(x.x)*x];
reflection[x_affineWeight]:=Function[y, y-2*(x.y)/(x.x)*x];

Expect["Reflection for finite weights", True,reflection[makeFiniteWeight[{1,0}]][makeFiniteWeight[{1,1}]]==makeFiniteWeight[{-1,1}]]

coroot[x_finiteWeight]:=2*x/(x.x);

Expect["Co root of [1,0]", True, coroot[makeFiniteWeight[{1,0}]]==makeFiniteWeight[{2,0}]]

cartanMatrix[rs_finiteRootSystem]:=Transpose[Outer[Dot,rs[simpleRoots],coroot/@rs[simpleRoots]]];

Expect["Cartan matrix of B2",{{2, -1}, {-2, 2}},cartanMatrix[makeSimpleRootSystem[B,2]]]

clear[weylGroupElement];

revApply[x_,f_]:=f[x];

(* Print[makeSimpleRootSystem[B,2][simpleRoot][2]] *)

rootSystemQ[rs_]:=MatchQ[rs,x_finiteRootSystem|x_affineRootSystem]

Expect["Predicate for finite and affine root systems",True,rootSystemQ[makeSimpleRootSystem[B,2]]]

weylGroupElement[x__Integer][rs_?rootSystemQ]:=Function[z,Fold[revApply,z,reflection /@ rs[simpleRoot]/@{x}]];

Expect["Weyl reflection s1 s2 s1 in algebra B2",True,
       weylGroupElement[1,2,1][makeSimpleRootSystem[B,2]][makeFiniteWeight[{1,0}]]==
       makeFiniteWeight[{-1,0}]]

fundamentalWeights[rs_finiteRootSystem]:=Plus@@(Inverse[cartanMatrix[rs]]*rs[simpleRoots]);

(*
Print/@fundamentalWeights[makeSimpleRootSystem[B,2]]

finiteWeight[$166]{1, 0}
                   1  1
finiteWeight[$167]{-, -}
                   2  2

Out[86]= {Null, Null}

Out[85]= {finiteWeight[$159], finiteWeight[$160]}
*)


rho[rs_finiteRootSystem]:=Plus@@fundamentalWeights[rs];

Print[rho[makeSimpleRootSystem[B,2]]]

                 3  1
finiteWeight[2, {-, -}]
                 2  2

Expect["Weyl vector for B2",True,makeFiniteWeight[{3/2,1/2}]==rho[makeSimpleRootSystem[B,2]]]

toFundamentalChamber[rs_finiteRootSystem][vec_finiteWeight]:=
    First[NestWhile[Function[v,
		       reflection[Scan[If[#.v<0,Return[#]]&,rs[simpleRoots]]][v]],
	      vec,
	      Head[#]=!=reflection[Null]&]]

Expect["To fundamental chamber",True,makeFiniteWeight[{1,1/2}]==toFundamentalChamber[makeSimpleRootSystem[B,2]][makeFiniteWeight[{-1,1/2}]]]

orbit[rs_finiteRootSystem][{weights__finiteWeight}]:=
    NestWhileList[
	Function[x,
		 Union[Flatten[Map[Function[y,
					    Map[reflection[#][y]&,Cases[rs[simpleRoots],z_ /; z.y>0]]],x]],SameTest->(#1==#2&)]],
	{weights},
	#=!={}&];
orbit[rs_finiteRootSystem][weight_finiteWeight]:=orbit[rs][{toFundamentalChamber[rs][weight]}];

positiveRoots[rs_finiteRootSystem]:=Map[-#&,Flatten[orbit[rs][Map[-#&,rs[simpleRoots]]]]]

orbit[makeSimpleRootSystem[B,2]][rho[makeSimpleRootSystem[B,2]]]

Print/@positiveRoots[makeSimpleRootSystem[B,2]]

Map[Print,orbit[makeSimpleRootSystem[B, 2]][rho[makeSimpleRootSystem[B, 2]]],2]

SubValues[finiteWeight]

weightSystem[rs_finiteRootSystem][higestWeight_finiteWeight]:=Module[{minusPosRoots=-positiveRoots[rs]},
									     NestWhileList[Function[x,Complement[
										 Cases[Flatten[Outer[Plus,minusPosRoots,x]],y_/;And@@(#.y>=0&/@rs[simpleRoots])]
										 ,x]],{higestWeight},#=!={}&]];

weightSystem[makeSimpleRootSystem[B,2]][makeFiniteWeight[{2,2}]]

Out[88]= {{finiteWeight[2, {2, 2}]}, 
 
>    {finiteWeight[2, {1, 1}], finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {0, 0}], finiteWeight[2, {1, 0}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {1, 1}]}, 
 
>    {finiteWeight[2, {0, 0}], finiteWeight[2, {1, 0}]}, {}}

freudenthalMultiplicities[rs_finiteRootSystem][highestWeight_finiteWeight]:=
    Module[{rh=rho[rs],weights,mults,c,insideQ,
	    posroots=positiveRoots[rs],
	    toFC=toFundamentalChamber[rs]},
	   weights=SortBy[ Rest[Flatten[weightSystem[rs][highestWeight]]], -#.rh&];
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

freudenthalMultiplicities[makeSimpleRootSystem[B,2]][makeFiniteWeight[{5,5}]]

Out[94]= mults$98

Out[90]= mults$96

Out[94][finiteWeight[2, {0, 0}]]

Out[96]= 6

Out[95]= 1

Out[93]= 1

Out[92]= 1

Out[91]= {finiteWeight[2, {1, 0}], finiteWeight[2, {0, 0}]}

orbitWithEps[rs_finiteRootSystem][weight_finiteWeight]:=Flatten[Most[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[rs][weight]]],1];

racahMultiplicities[rs_finiteRootSystem][highestWeight_finiteWeight]:=
    Module[{rh=rho[rs],weights,mults,c,insideQ,
	    fan,
	    toFC=toFundamentalChamber[rs]},
	   fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[rs][rh]]];
	   weights=Sort[ Rest[Flatten[weightSystem[rs][highestWeight]]], #1.rh>#2.rh&];
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 Plus@@(fan /. {x_finiteWeight,e_Integer}:> If[insideQ[v+x],-e*mults[toFC[v+x]],0])],
		weights];
	   mults]

racahMultiplicities[makeSimpleRootSystem[B,2]][makeFiniteWeight[{5,5}]]

Out[111]= mults$102

Out[99]= mults$100

(Out[111][#]&) /@ keys[Out[111]]

Out[112]= {4, 2, 2, 6, 1, 2, 3, 1, 4, 1, 5, 4, 1, 2, 3, 3, 1, 5, 1, 3, 2}

Out[105]= {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

Out[104]= {finiteWeight[2, {5, 3}], finiteWeight[2, {4, 1}], 
 
>    finiteWeight[2, {0, 0}], finiteWeight[2, {5, 4}], 
 
>    finiteWeight[2, {2, 2}], finiteWeight[2, {1, 0}], 
 
>    finiteWeight[2, {3, 2}], finiteWeight[2, {4, 0}], 
 
>    finiteWeight[2, {4, 4}], finiteWeight[2, {5, 5}], 
 
>    finiteWeight[2, {5, 2}], finiteWeight[2, {3, 0}], 
 
>    finiteWeight[2, {4, 3}], finiteWeight[2, {3, 1}], 
 
>    finiteWeight[2, {1, 1}], finiteWeight[2, {3, 3}], 
 
>    finiteWeight[2, {5, 0}], finiteWeight[2, {2, 0}], 
 
>    finiteWeight[2, {5, 1}], finiteWeight[2, {2, 1}], 
 
>    finiteWeight[2, {4, 2}]}

Out[103]= 0

Out[102]= mults$100[finiteWeight[2, {4, 5}]]

Out[101]= 1

Out[100]= 0


b2=makeSimpleRootSystem[B,2];
rh=rho[b2];
rs=b2;

fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[rs][rh]]];

fan

Out[108]= orbitWithEps[finiteRootSystem[2, 
 
>      {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}]][]

makeUntwistedAffineRootSystem[fs_finiteRootSystem]:=Module[{rs=affineRootSystem[Unique[]],ar},
							   rs[finiteRootSystem]=fs;
							   ar=makeAffineWeight[#,1,0]&/@fs[simpleRoots];
							   rs[simpleRoots]=Append[ar,makeAffineWeight[-highestRoot[fs],1,1]];
							   rs[simpleRoot][0]=rs[[-1]];
							   rs[simpleRoot][n_Integer]:=rs[simpleRoots][[n]];
							   rs]



