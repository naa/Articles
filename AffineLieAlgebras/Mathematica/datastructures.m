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

Expect::"usage"="Expect[comment,value,expression] is used for unit tests. If expression != value it throws exception.";
Expect[ description_, val_, expr_ ] := 
If[
    val != expr,
    Throw[
        StringJoin[ ToString[description],": GOT UNEXPECTED VALUE ", ToString[expr], 
        " INSTEAD OF ", ToString[val] ]
        , "assertion exception"
    ]
]

standardBase::"usage"=
    "standardBase[coordinates] represents weight vector in standard (Bourbaki) coordinates.\
    E.g. simple roots for B2 are standardBase[1,-1] and standardBase[0,1]"

(* finiteWeight[dimension,coordinates] *)

finiteWeight::"usage"=
    "finiteWeight[dimension_?NumberQ,coordinates_standardBase] represents\ 
    vector in weight space of finite-dimensional Lie algebra.\n 
    finiteWeight[dimension] returns dimension of the space, where weight vector is embedded\ 
    (i.e. for sl_n it is n+1).\n 
     finiteWeight[standardBase] returns standard base coordinates of weight of finite-dimensional Lie algebra";

finiteWeight/:x_finiteWeight[dimension]:=x[[1]];
finiteWeight/:x_finiteWeight[standardBase]:=x[[2]];

makeFiniteWeight::"usage"=
    "makeFiniteWeight[{coordinates__?NumberQ}] creates finiteWeight with given coordinates in standard base";

makeFiniteWeight[{coordinates__?NumberQ}]:=finiteWeight @@ {Length[{coordinates}],{coordinates}}

Expect["Dimension equals to length",3,makeFiniteWeight[{1,2,3}][dimension]]

Dot::"usage"=Dot::"usage" <> "\n It is defined for weights of finite and affine Lie algebras";

finiteWeight/:x_finiteWeight . y_finiteWeight/;x[dimension]==y[dimension]:=x[standardBase].y[standardBase]

Expect["Scalar product for vectors from weigth space of finite-dimensional Lie algebras",
       10,makeFiniteWeight[{1,2,3}].makeFiniteWeight[{3,2,1}]]

Expect["Scalar product for vectors from different spaces are left unevaluated",True,
       MatchQ[makeFiniteWeight[{1,2,3}].
	      makeFiniteWeight[{3,2,1,2}],x_finiteWeight . y_finiteWeight]]

Plus::"usage"=Plus::"usage" <> "\n It is defined for weights of finite and affine Lie algebras";

finiteWeight/:x_finiteWeight+y_finiteWeight/;x[dimension]==y[dimension]:=
    makeFiniteWeight[x[standardBase]+y[standardBase]]

Expect["Plus for finite-dimensional weights",{2,4,6},(makeFiniteWeight[{1,2,3}]+makeFiniteWeight[{1,2,3}])[standardBase]]

Expect["Plus product for vectors from different spaces are left unevaluated",True,
       MatchQ[makeFiniteWeight[{1,2,3}]+makeFiniteWeight[{3,2,1,2}],x_finiteWeight + y_finiteWeight]]

Equal::"usage"=Equal::"usage" <> "\n It is defined for weights of finite and affine Lie algebras";

finiteWeight/:x_finiteWeight==y_finiteWeight:=x[standardBase]==y[standardBase]

Expect["Equal for finite weights compares standard base representations", 
       False,
       makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{1,3,2}]]

Expect["Equal for finite weights compares standard base representations", 
       True,
       makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{1,2,3}]]

Expect["Equal for finite weights compares standard base representations", 
       False,
       makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{1,2,3,4}]]


Times::"usage"=Times::"usage"<>"\n\n Mutliplication by numbers is defined for the elements of weight space of affine and finite-dimensional Lie algebras.\ 
    \n For example makeFiniteWeight[{1,2,3}]*2==makeFiniteWeight[{2,4,6}]]";

finiteWeight/:0*y_finiteWeight:=makeFiniteWeight[0*y[standardBase]];

finiteWeight/:x_?NumberQ*y_finiteWeight:=makeFiniteWeight[x*y[standardBase]];

Expect["Multiplication by scalar", True,makeFiniteWeight[{1,2,3}]*2==makeFiniteWeight[{2,4,6}]]

Expect["Multiplication by scalar", True,2*makeFiniteWeight[{1,2,3}]==makeFiniteWeight[{2,4,6}]]


affineWeight::"usage"=
    "affineWeight[dimension_?NumberQ,fw_finiteWeigt,level_?NumberQ, grade_?NumberQ] represents\ 
    vector in weight space of affine Lie algebra.\n 
    affineWeight[dimension] returns dimension of the space of real roots, where finite weight vector is embedded\ 
    (i.e. for sl_n it is n+1).\n 
    affineWeight[finitePart] returns finite part of weight as finiteWeight structure\n
    affineWeight[level] returns level of affine weight\n
    affineWeight[grade] returns grade of affine weight";
affineWeight/:x_affineWeight[dimension]:=x[[1]];
affineWeight/:x_affineWeight[finitePart]:=x[[2]];
affineWeight/:x_affineWeight[level]:=x[[3]];
affineWeight/:x_affineWeight[grade]:=x[[4]];

makeAffineWeight::"usage"= 
    "makeAffineWeight[fw_finiteWeight,level_?NumberQ,grade_?NumberQ] creates affine weight with the given finite part fw, level and grade\n
    makeAffineWeight[{fw__?NumberQ},level_?NumberQ,grade_?NumberQ] creates affine weight with the given finite part fw, level and grade
    ";
makeAffineWeight[fp_finiteWeight,lev_?NumberQ,gr_?NumberQ]:=affineWeight[fp[dimension],fp,lev,gr];
makeAffineWeight[{fp__?NumberQ},lev_?NumberQ,gr_?NumberQ]:=makeAffineWeight[makeFiniteWeight[{fp}],lev,gr]


Expect["Affine weight has the same real dimension as the finite-dimensional part, since we hold level and grade separetely",
       True,makeAffineWeight[makeFiniteWeight[{1,2,3,4,5}],1,2][dimension]==Length[{1,2,3,4,5}]]

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

keys::"usage"="keys[hashtable] gives all the keys in hashtable";
keys = DownValues[#,Sort->False][[All,1,1,1]]&;


finiteRootSystem::"usage"=
    "finiteRootSystem[rank_Integer,{roots_finiteWeight}] represents root system of finite-dimensional Lie algebra.\n
    finiteRootSystem[rank] returns rank of the root system\n
    finiteRootSystem[simpleRoots] returns unsorted list of simple roots in the root system.\n
    finiteRootSystem[simpleRoot][n_Integer] returns n'th simple root"

makeFiniteRootSystem::"usage"=
    "makeFiniteRootSystem[{roots__finiteWeight}] creates finiteRootSystem structure from the List of simple roots\n
    makeFiniteRootSystem[{roots__finiteWeight}] creates  finiteRootSystem structure from the List of simple roots\n
    given as the lists of coordinates";

makeFiniteRootSystem[{roots__finiteWeight}]:=finiteRootSystem[Length[{roots}],{roots}]

finiteRootSystem/:x_finiteRootSystem[rank]:=x[[1]];
finiteRootSystem/:x_finiteRootSystem[simpleRoots]:=x[[2]];
finiteRootSystem/:x_finiteRootSystem[simpleRoot][n_Integer]:=x[[2]][[n]];

makeFiniteRootSystem[{roots__List}]:=makeFiniteRootSystem[makeFiniteWeight/@{roots}]

makeSimpleRootSystem::"usage"=
    "makeSimpleRootSystem[series, rank] creates root system of the algebra from given series with given rank.\n
    Exceptional Lie algebras are not supported yet";

makeSimpleRootSystem[A,1]:=makeFiniteRootSystem[{{1}}];
makeSimpleRootSystem[A,r_Integer]:=makeFiniteRootSystem[makeFiniteWeight /@ Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,r},{j,1,r+1}]];
makeSimpleRootSystem[B,rank_Integer]:=makeFiniteRootSystem[Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],1]]];
makeSimpleRootSystem[C,rank_Integer]:=makeFiniteRootSystem[Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],2]]];
makeSimpleRootSystem[D,rank_Integer]:=makeFiniteRootSystem[Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Append[Table[0,{rank-2}],1],1]]];

Expect["B2:",True,makeSimpleRootSystem[B,2][simpleRoots][[1]]==makeFiniteWeight[{1,-1}]]

Expect["B2: rank",2,makeSimpleRootSystem[B,2][rank]]

checkGrade[rs_?rootSystemQ][w_finiteWeight]=True;
checkGrade[rs_?rootSystemQ][w_affineWeight]=(Abs[w[grade]]<rs[gradeLimit]) /. rs[gradeLimit]->10;



weightQ::"usage"=
    "Weight predicate"

weightQ[x_finiteWeight]=True;
weightQ[x_affineWeight]=True;
weightQ[_]=False;

reflection::"usage"=
    "reflection[x_finiteWeight] represents reflection in the hyperplane, orthogonal to the given root x\n
    reflection[x_affineWeight] represents affine reflection in the hyperplane, orthogonal to the given root x\n
    reflection[x][y] reflects y in the hyperplane, orthogonal to x"

reflection[x_?weightQ]:=Function[y, y-2*(x.y)/(x.x)*x];

Expect["Reflection for finite weights", True,reflection[makeFiniteWeight[{1,0}]][makeFiniteWeight[{1,1}]]==makeFiniteWeight[{-1,1}]]

coroot::"usage"="returns coroot of given root"
    
coroot[x_?weightQ]:=2*x/(x.x);

Expect["Co root of [1,0]", True, coroot[makeFiniteWeight[{1,0}]]==makeFiniteWeight[{2,0}]]

Expect["Co root of affine [1,0]", True, coroot[makeAffineWeight[{1,0},1,0]]==makeAffineWeight[{2,0},2,0]]

cartanMatrix::"usage"=
    "cartanMatrix[rs_?rootSystemQ] returns Cartan matrix of root system rs";
cartanMatrix[rs_?rootSystemQ]:=Transpose[Outer[Dot,rs[simpleRoots],coroot/@rs[simpleRoots]]];

Expect["Cartan matrix of B2",{{2, -1}, {-2, 2}},cartanMatrix[makeSimpleRootSystem[B,2]]]

revApply[x_,f_]:=f[x];

(* Print[makeSimpleRootSystem[B,2][simpleRoot][2]] *)

rootSystemQ::"usage"=
    "Predicate for root system data structures";
rootSystemQ[rs_]:=MatchQ[rs,x_finiteRootSystem|x_affineRootSystem]

Expect["Predicate for finite and affine root systems",True,rootSystemQ[makeSimpleRootSystem[B,2]]]

weylGroupElement::"usage"=
    "weylGroupElement[x__Integer][rs_?rootSystemQ] represents element of Weyl group composed of basic reflecetions x \n
    of root system rs. Example: weylGroupElement[1,2,1][makeSimpleRootSystem[B,2]] = s1 s2 s1";
clear[weylGroupElement];
weylGroupElement[x__Integer][rs_?rootSystemQ]:=Function[z,Fold[revApply,z,reflection /@ rs[simpleRoot]/@{x}]];

Expect["Weyl reflection s1 s2 s1 in algebra B2",True,
       weylGroupElement[1,2,1][makeSimpleRootSystem[B,2]][makeFiniteWeight[{1,0}]]==
       makeFiniteWeight[{-1,0}]]

fundamentalWeights::"usage"=
    "fundamentalWeights[rs_?rootSystemQ] calculates fundamental weights of given root system rs";
fundamentalWeights[rs_finiteRootSystem]:=Plus@@(Inverse[cartanMatrix[rs]]*rs[simpleRoots]);

Expect["Fundamental weights for B2", True, {makeFiniteWeight[{1,0}],makeFiniteWeight[{1/2,1/2}]}==fundamentalWeights[makeSimpleRootSystem[B,2]]]

Expect["Fundamental weights for A1", True, {makeFiniteWeight[{1/2}]}==fundamentalWeights[makeSimpleRootSystem[A,1]]]

Expect["Fundamental weights for A2", True, {makeFiniteWeight[{2/3,-1/3,-1/3}],makeFiniteWeight[{1/3,1/3,-2/3}]}==fundamentalWeights[makeSimpleRootSystem[A,2]]]

rho::"usage"=
    "rho[rs_?rootSystemQ] - Weyl vector of root system rs (sum of fundamental weights)";
rho[rs_?rootSystemQ]:=Plus@@fundamentalWeights[rs];

Expect["Weyl vector for B2",True,makeFiniteWeight[{3/2,1/2}]==rho[makeSimpleRootSystem[B,2]]]

toFundamentalChamber::"usage"=
    "toFundamentalChamber[rs_?rootSystemQ][vec_?weightQ] acts on a weight vec by simple reflections of rs till it gets to main Weyl chamber";
toFundamentalChamber[rs_?rootSystemQ][vec_?weightQ]:=
    First[NestWhile[Function[v,
		       reflection[Scan[If[#.v<0,Return[#]]&,rs[simpleRoots]]][v]],
	      vec,
	      Head[#]=!=reflection[Null]&]]

Expect["To fundamental chamber",True,makeFiniteWeight[{1,1/2}]==toFundamentalChamber[makeSimpleRootSystem[B,2]][makeFiniteWeight[{-1,1/2}]]]

Expect["We can use this and other functions for mapping",True,
       Map[toFundamentalChamber[makeSimpleRootSystem[B,2]],{makeFiniteWeight[{-1,-1}],makeFiniteWeight[{-2,-1}]}]==
       {makeFiniteWeight[{1, 1}], makeFiniteWeight[{2, 1}]}]



partialOrbit::"usage"=
    "partialOrbit[rs_finiteRootSystem][{weights__finiteWeight}] constructs Weyl partial orbit of given set of weights. \n
    It starts with the list of weights and reflects them away from main Weyl chamber.\n
    For w in fundamental chamber it gives the whole orbit.";
partialOrbit[rs_?rootSystemQ][{weights__?weightQ}]:=
    Most[NestWhileList[
	Function[x,
		 Union[
		     Flatten[
			 Map[
			     Function[y,
				      Map[reflection[#][y]&,
					  Cases[rs[simpleRoots],
						z_ /; And[z.y>0,
							  checkGrade[rs][reflection[z][y]]]
					       ]
					 ]
				     ],
			     x]
			    ],
		     SameTest->(#1==#2&)]
		],
	{weights},
	#=!={}&]];

Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["Partial orbit test", True, partialOrbit[b2][{rho[b2]}]==
	      {{makeFiniteWeight[{3/2,1/2}]},{makeFiniteWeight[{1/2,3/2}],makeFiniteWeight[{3/2,-1/2}]},
	       {makeFiniteWeight[{-1/2,3/2}],makeFiniteWeight[{1/2,-3/2}]},
	       {makeFiniteWeight[{-3/2,1/2}],makeFiniteWeight[{-1/2,-3/2}]},
	       {makeFiniteWeight[{-3/2,-1/2}]}}]]

orbit::"usage"=
    "orbit[rs_?rootSystemQ][wg_?weightQ] returns the weight of Weyl orbit, containing given weight wg. \n
    Orbit is given as the list of lists starting with the weight in dominant Weyl chamber\n
    orbit[rs_?rootSystemQ][{wg_?weightQ}] works for a list of weights";

orbit[rs_?rootSystemQ][weight_?weightQ]:=partialOrbit[rs][{toFundamentalChamber[rs][weight]}];
orbit[rs_?rootSystemQ][{weights__?weightQ}]:=partialOrbit[rs][toFundamentalChamber[rs] /@ {weights}];

Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["orbit is equivalent to partial orbit",True, partialOrbit[b2][{rho[b2]}]== orbit[b2][rho[b2]]]]


positiveRoots::"usage"=
    "positiveRoots[rs_?rootSystemQ] returns positive roots of root system rs";
(* Ugly hack *)
positiveRoots[rs_affineRootSystem]:=Join[Map[-#&,Flatten[partialOrbit[rs][Map[-#&,rs[simpleRoots]]]]],
					 Join@@Table[NestWhileList[#+makeAffineWeight[zeroWeight[rs[finiteRootSystem]],0,1]&,zeroWeight[rs],checkGrade[rs][#]&],{rs[rank]}]];
positiveroots[rs_?rootSystemQ]:=Map[-#&,Flatten[partialOrbit[rs][Map[-#&,rs[simpleRoots]]]]]


Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["Positive roots of B2",True, positiveRoots[b2]=={makeFiniteWeight[{1,-1}],makeFiniteWeight[{0,1}],
							       makeFiniteWeight[{1,0}],makeFiniteWeight[{1,1}]}]]

(*orbit[makeSimpleRootSystem[B,2]][rho[makeSimpleRootSystem[B,2]]]

Print/@positiveRoots[makeSimpleRootSystem[B,2]]

Map[Print,orbit[makeSimpleRootSystem[B, 2]][rho[makeSimpleRootSystem[B, 2]]],2]

SubValues[finiteWeight]
*)

weightSystem::"usage"=
    "weightSystem[rs_?rootSystemQ][higestWeight_?weightQ] returns the set of dominant weights in the highest weight module. \n
    The list is split in pieces by number of root substractions";
weightSystem[rs_?rootSystemQ][higestWeight_?weightQ]:=Module[{minusPosRoots=-positiveRoots[rs]},
							     Most[NestWhileList[Function[x,Complement[
										 Cases[Flatten[Outer[Plus,minusPosRoots,x]],y_/;
										       And[checkGrade[rs][y],And@@(#.y>=0&/@rs[simpleRoots])]]
										 ,x]],{higestWeight},#=!={}&]]];


Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["Weights of [2,1] module of B2",True, weightSystem[b2][makeFiniteWeight[{2,1}]]==
	      {{makeFiniteWeight[{2,1}]},
	      {makeFiniteWeight[{1,0}],makeFiniteWeight[{1,1}],makeFiniteWeight[{2,0}]},
	      {makeFiniteWeight[{0,0}]}}]]

freudenthalMultiplicities::"usage"=
    "freudenthalMultiplicities[rs_finiteRootSystem][highestWeight_finiteWeight] returns hashtable with the multiplicities of 
    weights in the highest weight module";
freudenthalMultiplicities[rs_?rootSystemQ][highestWeight_?weightQ]:=
    Module[{rh=rho[rs],weights,mults,c,insideQ,
	    posroots=positiveRoots[rs],
	    toFC=toFundamentalChamber[rs]},
	   weights=SortBy[ Rest[Flatten[weightSystem[rs][highestWeight]]], -#.rh&];
	   c:=(#+rh).(#+rh)&;
	   mults[highestWeight]=1;
	   insideQ:=And[checkGrade[rs][#],IntegerQ[mults[toFC[#]]]]&;
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

Module[{b2=makeSimpleRootSystem[B,2],fm},
       fm=freudenthalMultiplicities[b2][makeFiniteWeight[{2,1}]];
       Expect["Mutliplicites of [2,1] B2 representation",{1, 1, 2, 3, 3},
	      fm[#]&/@keys[fm]]]

(* freudenthalMultiplicities[makeSimpleRootSystem[B,2]][makeFiniteWeight[{5,5}]] *)
orbitWithEps::"usage"="returns orbit with the determinants of Weyl reflections";
orbitWithEps[rs_?rootSystemQ][weight_?weightQ]:=Flatten[Most[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[rs][weight]]],1];

racahMultiplicities[rs_?rootSystemQ][highestWeight_?weightQ]:=
    Module[{rh=rho[rs],weights,mults,c,insideQ,
	    fan,
	    toFC=toFundamentalChamber[rs]},
	   fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[rs][rh]]];
	   weights=Sort[ Rest[Flatten[weightSystem[rs][highestWeight]]], #1.rh>#2.rh&];
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 Plus@@(fan /. {x_?weightQ,e_Integer}:> If[insideQ[v+x],-e*mults[toFC[v+x]],0])],
		weights];
	   mults]

(* racahMultiplicities[makeSimpleRootSystem[B,2]][makeFiniteWeight[{5,5}]] *)

(* b2=makeSimpleRootSystem[B,2];
rh=rho[b2];
rs=b2; 

fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[rs][rh]]]; *)

higestRoot[rs_finiteRootSystem]:=toFundamentalChamber[rs][rs[simpleRoot][Ordering[(#.#&)/@rs[simpleRoots],-1][[1]]]]

Expect["Higest root for B2",makeFiniteWeight[{1, 1}],higestRoot[makeSimpleRootSystem[B,2]]]

makeAffineExtension[fs_finiteRootSystem]:=affineRootSystem[fs[rank],fs,makeAffineWeight[-higestRoot[fs],0,1],(makeAffineWeight[#,0,0]&)/@fs[simpleRoots]]

Append[makeAffineExtension[makeSimpleRootSystem[B,2]][simpleRoots],makeAffineExtension[makeSimpleRootSystem[B,2]][[3]]]

(* makeAffineExtension[makeSimpleRootSystem[B,2]][simpleRoot][0] *)

affineRootSystem/:rs_affineRootSystem[rank]:=rs[[1]]

affineRootSystem/:rs_affineRootSystem[finiteRootSystem]:=rs[[2]]

affineRootSystem/:rs_affineRootSystem[realRoots]:=rs[[4]]

affineRootSystem/:rs_affineRootSystem[imaginaryRoot]:=rs[[3]]

affineRootSystem/:rs_affineRootSystem[simpleRoots]:=Prepend[rs[[4]],rs[[3]]]

affineRootSystem/:rs_affineRootSystem[simpleRoot][0]:=rs[[3]];
affineRootSystem/:rs_affineRootSystem[simpleRoot][n_?NumberQ]/;n<=rs[rank]:=rs[simpleRoots][[n]];

zeroWeight[rs_finiteRootSystem]:=makeFiniteWeight[Table[0,{rs[rank]}]];
zeroWeight[rs_affineRootSystem]:=makeAffineWeight[zeroWeight[rs[finiteRootSystem]],0,0];

toFundamentalChamber[rs_affineRootSystem][vec_affineWeight]:=
    First[NestWhile[Function[v,
		       reflection[Scan[If[#.v<0,Return[#]]&,rs[simpleRoots]]][v]],
	      vec,
	      Head[#]=!=reflection[Null]&]]

marks[rs_affineRootSystem]:=Prepend[Inverse[cartanMatrix[rs[finiteRootSystem]]].(-2*#.rs[simpleRoot][0]/(#.#)&)/@rs[realRoots],1]

comarks[rs_affineRootSystem]:=marks[rs]*Map[#.#/2&,rs[simpleRoots]]

b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]


fundamentalWeights[rs_affineRootSystem]:=Map[makeAffineWeight[#[[1]],#[[2]],0]&,
					     Transpose[{Prepend[fundamentalWeights[rs[finiteRootSystem]],
								0*rs[finiteRootSystem][simpleRoot][1]],
							comarks[rs]}]]

(* fundamentalWeights[b2a] *)

orbit[rs_affineRootSystem][{weights__affineWeight},gradelimit_?NumberQ]:=
			   NestWhileList[
			       Function[x,
					Union[Flatten[
					    Map[Function[y,
							 Map[reflection[#][y]&,
							     Cases[rs[simpleRoots],z_ /; And[z.y>0,Abs[reflection[z][y][grade]]<gradelimit]]
							    ]],
						x]],
					      SameTest->(#1==#2&)]],
			       {weights},
			       #=!={}&]

orbitWithEps[rs_affineRootSystem][weight_affineWeight,gradelimit_?NumberQ]:=Flatten[Most[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[rs][{weight},gradelimit]]],1]

positiveRoots[rs_affineRootSystem,gradelimit_?NumberQ]:=Map[-#&,Flatten[orbit[rs][Map[-#&,rs[simpleRoots]],gradelimit]]]


racahMultiplicities[rs_affineRootSystem][highestWeight_affineWeight,gradelimit_?NumberQ]:=
    Module[{rh=rho[rs],weights,mults,c,insideQ,
	    fan,
	    toFC=toFundamentalChamber[rs]},
	   fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[rs][rh,gradelimit]]];
	   weights=Sort[ Rest[Flatten[weightSystem[rs][highestWeight,gradelimit]]], #1.rh>#2.rh&];
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 Plus@@(fan /. {x_finiteWeight,e_Integer}:> If[insideQ[v+x],-e*mults[toFC[v+x]],0])],
		weights];
	   mults]
				     
