clear[standardBase]

makeSimpleRootSystem[A,rank_Integer]:=standardBase @@@ Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank},{j,1,rank+1}];
makeSimpleRootSystem[B,rank_Integer]:=standardBase @@@ Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],1]];
makeSimpleRootSystem[C,rank_Integer]:=standardBase @@@ Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Table[0,{rank-1}],2]];
makeSimpleRootSystem[D,rank_Integer]:=standardBase @@@ Append[Table[If[i==j,1,If[i==j-1,-1,0]],{i,1,rank-1},{j,1,rank}],Append[Append[Table[0,{rank-2}],1],1]];

Unprotect[Dot];
Dot[x_standardBase,y_standardBase]:=standardBase @@ Dot[List@@x,List@@y]

Unprotect[Plus];
Plus[x_standardBase,y_standardBase]:=standardBase @@ Plus[List@@x,List@@y]

Unprotect[Times];
Times[x_,y_standardBase]:=standardBase @@ Times[x,List@@y]

reflection[x_standardBase]:=Function[y, y-2*(x.y)/(x.x)*x]

coroot[x_standardBase]:=2*x/(x.x)

cartanMatrix[{x__standardBase}]:=Transpose[Outer[Dot,{x},coroot/@{x}]]

weights[
coroot[b2[[2]]]

coroot/@b2



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