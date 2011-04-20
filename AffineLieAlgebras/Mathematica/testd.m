AppendTo[$Path,"/home/anton/study/2011/Articles/AffineLieAlgebras/Mathematica/"];

<<datastructures.m;

orbitWithEps[rs_?rootSystemQ][weight_?weightQ]:=Flatten[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[rs][weight]],1];

anomalousWeights[makeSimpleRootSystem[B,2]][makeFiniteWeight[{2,1}]][weights]

Out[21]= {finiteWeight[2, {2, 1}], finiteWeight[2, {0, 3}], 
 
>    finiteWeight[2, {2, -2}], finiteWeight[2, {-3, 3}], 
 
>    finiteWeight[2, {0, -4}], finiteWeight[2, {-5, 1}], 
 
>    finiteWeight[2, {-3, -4}], finiteWeight[2, {-5, -2}]}

Out[20]= {1, -1, -1, 1, 1, -1, -1, 1}

anomalousWeights[rs_?rootSystemQ][hweight_?weightQ]:=
    makeFormalElement @@ Transpose[{#[[1]]-rho[rs],#[[2]]}& /@ orbitWithEps[rs][hweight+rho[rs]]]


Out[18]= {1, -1, -1, 1, 1, -1, -1, 1}

Out[17]= {finiteWeight[2, {2, 1}], finiteWeight[2, {1, 2}], 
 
>    finiteWeight[2, {2, -1}], finiteWeight[2, {-1, 2}], 
 
>    finiteWeight[2, {1, -2}], finiteWeight[2, {-2, 1}], 
 
>    finiteWeight[2, {-1, -2}], finiteWeight[2, {-2, -1}]}

Out[11]= {finiteWeight[2, {2, 1}], finiteWeight[2, {1, 2}], 
 
>    finiteWeight[2, {2, -1}], finiteWeight[2, {-1, 2}], 
 
>    finiteWeight[2, {1, -2}], finiteWeight[2, {-2, 1}], 
 
>    finiteWeight[2, {-1, -2}]}

Out[10]= {1, -1, -1, 1, 1, -1, -1}

Rest[orbit[makeSimpleRootSystem[B,2]][makeFiniteWeight[{2,1}]]]

Out[15]= {{finiteWeight[2, {1, 2}], finiteWeight[2, {2, -1}]}, 
 
>    {finiteWeight[2, {-1, 2}], finiteWeight[2, {1, -2}]}, 
 
>    {finiteWeight[2, {-2, 1}], finiteWeight[2, {-1, -2}]}, 
 
>    {finiteWeight[2, {-2, -1}]}}

Out[14]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 2}], finiteWeight[2, {2, -1}]}, 
 
>    {finiteWeight[2, {-1, 2}], finiteWeight[2, {1, -2}]}, 
 
>    {finiteWeight[2, {-2, 1}], finiteWeight[2, {-1, -2}]}}

orbit[makeSimpleRootSystem[B,2]][makeFiniteWeight[{2,1}]]

Out[13]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 2}], finiteWeight[2, {2, -1}]}, 
 
>    {finiteWeight[2, {-1, 2}], finiteWeight[2, {1, -2}]}, 
 
>    {finiteWeight[2, {-2, 1}], finiteWeight[2, {-1, -2}]}, 
 
>    {finiteWeight[2, {-2, -1}]}}

orbitWithEps[makeSimpleRootSystem[B,2]][makeFiniteWeight[{2,1}]]

Out[9]= {{finiteWeight[2, {2, 1}], 1}, {finiteWeight[2, {1, 2}], -1}, 
 
>    {finiteWeight[2, {2, -1}], -1}, {finiteWeight[2, {-1, 2}], 1}, 
 
>    {finiteWeight[2, {1, -2}], 1}, {finiteWeight[2, {-2, 1}], -1}, 
 
>    {finiteWeight[2, {-1, -2}], -1}}

Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["Partial orbit test", True, partialOrbit[b2][{rho[b2]}]==
	      {{makeFiniteWeight[{3/2,1/2}]},{makeFiniteWeight[{1/2,3/2}],makeFiniteWeight[{3/2,-1/2}]},
	       {makeFiniteWeight[{-1/2,3/2}],makeFiniteWeight[{1/2,-3/2}]},
	       {makeFiniteWeight[{-3/2,1/2}],makeFiniteWeight[{-1/2,-3/2}]},
	       {makeFiniteWeight[{-3/2,-1/2}]}}]]


Out[7]= {{finiteWeight[2, {1, 0}], 1}, {finiteWeight[2, {0, 1}], -1}, 
 
>    {finiteWeight[2, {0, -1}], 1}}

Out[6]= {{finiteWeight[2, {1, 1}], 1}, {finiteWeight[2, {1, -1}], -1}, 
 
>    {finiteWeight[2, {-1, 1}], 1}}

Out[5]= {1, -1, 1}

Out[4]= {finiteWeight[2, {1, 1}], finiteWeight[2, {1, -1}], 
 
>    finiteWeight[2, {-1, 1}]}

Out[3]= formalElement[table$147]

Mathematica 8.0 for Linux x86 (32-bit)
Copyright 1988-2010 Wolfram Research, Inc.

fn=fan[makeSimpleRootSystem[B,2],makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]]

                  1  1
{finiteWeight[2, {-, -}], finiteWeight[2, {1, 1}]}
                  2  2
{2, 0}

Out[97]= formalElement[table$910]

                  1  1
{finiteWeight[2, {-, -}], finiteWeight[2, {1, 1}]}
                  2  2
{2, 0}

makeFormalElement[{makeFiniteWeight[{1/2,1/2}],makeFiniteWeight[{1/2,1/2}]}]

Expand[makeFormalElement[{zeroWeight[makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]]},{1}]*(1-Exp[makeFiniteWeight[{1/2,1/2}]])^2][multiplicities]

orbitWithEps[

Out[78]= {1, 0, 0}

ze=makeFormalElement[{zeroWeight[makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]]}]


                                                    1  1
Out[77]= {finiteWeight[2, {0, 0}], finiteWeight[2, {-, -}], 
                                                    2  2
 
>    finiteWeight[2, {1, 1}]}

(Expand[(1-Exp[makeFiniteWeight[{1/2,1/2}]])^2][[2]]*ze)[multiplicities]

Out[96]= {-2}

(2*ze*Exp[makeFiniteWeight[{1/2,1/2}]])[multiplicities]

Out[95]= {2}

ze[makeFiniteWeight[{0,0}]]

?MapAt

MapAt[f, expr, n] applies f to the element at position n
      in expr. If n is negative, the position is counted from the end. 

        MapAt[f, expr, {i, j, ...}]
         applies f to the part of expr
           at position {i, j, ...}
           . MapAt[f, expr, {{i , j , ...}, {i , j , ...}, ...}]
                               1   1          2   2
             applies f to parts of expr at several positions. 

Out[92]= 1

Out[91]= 0

Out[90]= {0}

Out[89]= {2}

Out[88]= formalElement[table$851]

Out[87]= {0}

                           1  1
Out[86]= {finiteWeight[2, {-, -}]}
                           2  2

Out[85]= formalElement[table$846]

             finiteWeight[2, {1/2, 1/2}]
Out[84]= -2 E

                finiteWeight[2, {1/2, 1/2}]    finiteWeight[2, {1, 1}]
Out[79]= 1 - 2 E                            + E

Out[76]= formalElement[table$824]

Out[75]= {1}

Out[74]= {finiteWeight[2, {0, 0}]}

Out[73]= formalElement[table$814]

Out[72]= makeFormalElement[{zeroWeight[subs]}, {1}]

Out[71]= formalElement[table$813]

fn[weights]

                                                    1  1
Out[98]= {finiteWeight[2, {0, 0}], finiteWeight[2, {-, -}], 
                                                    2  2
 
>    finiteWeight[2, {1, 1}]}

                                                    1  1
Out[68]= {finiteWeight[2, {0, 0}], finiteWeight[2, {-, -}], 
                                                    2  2
 
>    finiteWeight[2, {1, 1}]}

fn[multiplicities]

b2=makeSimpleRootSystem[B,2]

orbitWithEps[b2][weight[b2][1,1]]

                             3  1                          1  3
Out[104]= {{finiteWeight[2, {-, -}], 1}, {finiteWeight[2, {-, -}], -1}, 
                             2  2                          2  2
 
                       3    1                              1   3
>    {finiteWeight[2, {-, -(-)}], -1}, {finiteWeight[2, {-(-), -}], 1}, 
                       2    2                              2   2
 
                       1    3                             3   1
>    {finiteWeight[2, {-, -(-)}], 1}, {finiteWeight[2, {-(-), -}], -1}, 
                       2    2                             2   2
 
                         1     3
>    {finiteWeight[2, {-(-), -(-)}], -1}}
                         2     2

Out[99]= {1, -2, 1}

Out[69]= {1, 0, 0}

Out[67]= formalElement[table$768]

Out[66]= formalElement[table$752]

formalElement[table$707]

makeFormalElement[positiveRoots[makeSimpleRootSystem[B,2]]]

Scan::nonopt: Options expected (instead of finiteWeight[2, {1, 0}]
    ) beyond position 3 in Scan[(res$669[hashtable][#1] = 
        res$669[#1] + 1) & , finiteWeight[2, {1, -1}], 
     finiteWeight[2, {0, 1}], finiteWeight[2, {1, 1}], 
     finiteWeight[2, {1, 0}]]. An option must be a rule or a list of rules.

Out[62]= formalElement[table$670]


makeFormalElement[positiveRoots[makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]]]

Out[47]= formalElement[table$449]

Out[46]= {finiteWeight[2, {1, 1}]}

Scan::nonopt: Options expected (instead of finiteWeight[2, {1, 1}]
    ) beyond position 3 in Scan[(res$443[hashtable][#1] := 
                                              1  1
        res$443[#1] + 1) & , finiteWeight[2, {-, -}], 
                                              2  2
                      1  1
     finiteWeight[2, {-, -}], finiteWeight[2, {1, 1}]]. An option must be a
                      2  2
     rule or a list of rules.
formalElement[table$444] + makeFormalElement[{2, {1, 1}}, 
 
>    {-1 - formalElement[table$446][2], 
 
>     -1 - formalElement[table$446][{1, 1}]}]

                    weights
Out[45]= Power[1 - E       , (formalElement[table$444] + 
 
>        makeFormalElement[{2, {1, 1}}, 
 
>         {-1 - formalElement[table$446][2], 
 
>          -1 - formalElement[table$446][{1, 1}]}])[weights]] 
 
>    formalElement[table$447]

subElement[makeFormalElement[positiveRoots[makeSimpleRootSystem[B,2]],{1,2,3,4}],{makeFiniteWeight[{1,-1}]}][multiplicities]

fe=makeFormalElement[positiveRoots[makeSimpleRootSystem[B,2]],{1,2,3,4}]

projection[makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]][fe]



Out[41]= formalElement[table$435]

Expand[fe*(1-Exp[makeFiniteWeight[{1,-1}]])^3][weights]

?FoldList

FoldList[f, x, {a, b, ...}] gives {x, f[x, a], f[f[x, a], b], ...}. 

Fold[f, x, list] gives the last element of FoldList[f, x, list]. 

Out[32]= {finiteWeight[2, {0, 1}], finiteWeight[2, {1, -1}], 
 
>    finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>    finiteWeight[2, {2, -2}], finiteWeight[2, {2, -1}], 
 
>    finiteWeight[2, {2, 0}], finiteWeight[2, {3, -3}], 
 
>    finiteWeight[2, {3, -2}], finiteWeight[2, {3, -1}], 
 
>    finiteWeight[2, {4, -4}], finiteWeight[2, {4, -3}], 
 
>    finiteWeight[2, {4, -2}]}

Out[31]= formalElement[table$315]

Out[30]= formalElement[table$303]

                 finiteWeight[2, {1, -1}]      finiteWeight[2, {2, -2}]
Out[29]= (1 - 3 E                         + 3 E                         - 
 
        finiteWeight[2, {3, -3}]
>      E                        ) formalElement[table$291]

                 finiteWeight[2, {1, -1}]      finiteWeight[2, {2, -2}]
Out[28]= (1 - 3 E                         + 3 E                         - 
 
        finiteWeight[2, {3, -3}]
>      E                        ) formalElement[table$291]

               finiteWeight[2, {1, -1}] 3
Out[27]= (1 - E                        )  formalElement[table$291]

Out[26]= formalElement[table$291]

(makeFormalElement[{makeAffineWeight[{1,1},1,10]},{1}]*Exp[makeAffineWeight[{1,2},1,1]]*2)[weights]

Out[25]= {}

Syntax::sntxb: Expression cannot begin with 
    "[makeFormalElement[{makeAffineWeight[{1,1},1,10]},{1}]*Exp[makeAffineWeig
      ht[{1,2},1,1]]*2][weights]".

Out[24]= formalElement[table$285]

Out[23]= formalElement[table$282]

Out[22]= formalElement[table$280]

Out[20]= makeFormalElement[{}, {}]

Out[18]= subElement[formalElement[table$219], {}]

Out[16]= subElement[formalElement[table$189], checkGrade[]]

Out[15]= {affineWeight[2, finiteWeight[2, {1, 1}], 1, 10]}

Out[14]= formalElement[table$187]

Out[12]= {1}

Out[11]= {finiteWeight[2, {1, -1}]}

Out[9]= makeFormalElement[finiteWeight[2, {1, -1}], 
 
>     finiteWeight[formalElement[table$153][2], 
 
>      formalElement[table$153][{1, -1}]]][weights]

Out[8]= makeFormalElement[finiteWeight[2, {1, -1}], 
 
>     finiteWeight[formalElement[table$152][2], 
 
>      formalElement[table$152][{1, -1}]]][weights]

Out[7]= {1, 2, 3, 4}

Out[6]= {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}], 
 
>    finiteWeight[2, {1, 1}], finiteWeight[2, {1, 0}]}

Out[5]= formalElement[table$149]

fan[makeAffineExtension[makeSimpleRootSystem[A,2]],makeAffineExtension[makeSimpleRootSystem[A,1]]]

?Complement

?Select

Select[list, crit] picks out all elements e
                                           i
     of list for which crit[e ] is True. Select[list, crit, n]
                             i
        picks out the first n elements for which crit[e ] is True. 
                                                       i

Complement[e   , e , e , ...] gives the elements in e
            all   1   2                              all
     which are not in any of the e . 
                                  i

[Mathematica segmentation fault.]

$RecursionLimit::reclim: Recursion depth of 256 exceeded.

Out[8]= (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          (((Hold[r$144 - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144] - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                            Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>                            0, 10]] r$144) - 
 
>                          Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 
 
>                            10]] r$144) - 
 
>                        Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 
 
>                           10]] r$144) - 
 
>                      Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 
 
>                         10]] r$144) - 
 
>                    Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 10]] 
 
>                     r$144) - Exp[affineWeight[3, 
 
>                     finiteWeight[3, {0, 0, 0}], 0, 10]] r$144) - 
 
>                Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 10]] r$144
 
>                ) - Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 10]] 
 
>               r$144) - Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 
 
>               10]] r$144) - Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 
 
>             0, 10]] r$144) - Exp[affineWeight[3, 
 
>           finiteWeight[3, {0, 0, 0}], 0, 10]] r$144) - 
 
>      Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 10]] r$144) - 
 
>    Exp[affineWeight[3, finiteWeight[3, {0, 0, 0}], 0, 10]] r$144

[Mathematica segmentation fault.]

[Mathematica segmentation fault.]


Throw::nocatch: 
   Uncaught Throw[checkGrade affine Lie algebra: GOT UNEXPECTED VALUE True
      INSTEAD OF False, assertion exception] returned to top level.

Out[163]= Hold[Throw[checkGrade affine Lie algebra: GOT UNEXPECTED VALUE True\
 
>      INSTEAD OF False, assertion exception]]

Module[{b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]},
       b2a[gradeLimit]=2;
       positiveRoots[b2a]
       Expect["Positive roots of affine B2",True, positiveRoots[b2a]=={makeAffineWeight[{-1, -1}, 0, 1], 
								       makeAffineWeight[{1, -1}, 0, 0], 
								       makeAffineWeight[{0, 1}, 0, 0], 
								       makeAffineWeight[{1, 1}, 0, 0], 
								       makeAffineWeight[{1, 0}, 0, 0], 
								       makeAffineWeight[{0, 0}, 0, 0], 
								       makeAffineWeight[{0, 0}, 0, 1], 
								       makeAffineWeight[{0, 0}, 0, 0], 
								       makeAffineWeight[{0, 0}, 0, 1]}]]


pr=positiveRoots[makeAffineExtension[makeSimpleRootSystem[B,2]]];
fe=makeFormalElement[pr,Table[1,{Length[pr]}]];

Union[pr]



b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]

Out[3]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]

b2a[gradeLimit]=2

Out[4]= 2

fan[b2a,b2a]

$RecursionLimit::reclim: Recursion depth of 256 exceeded.

Out[5]= (((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          ((((((((((((((((((((((((((((((((((((((((((((((
 
>                          (((Hold[r$143 - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143] - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                            2]] r$143) - 
 
>                          Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                           2]] r$143) - 
 
>                        Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] 
 
>                         r$143) - 
 
>                      Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] 
 
>                       r$143) - 
 
>                    Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143
 
>                    ) - Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] 
 
>                   r$143) - Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 
 
>                   2]] r$143) - 
 
>              Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143) - 
 
>            Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143) - 
 
>          Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143) - 
 
>        Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143) - 
 
>      Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143) - 
 
>    Exp[affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]] r$143

[Mathematica segmentation fault.]

Out[10]= 2

(-fe*Exp[makeAffineWeight[{1,1},1,1]])[multiplicities]

Out[8]= {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
 
>    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
 
>    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
 
>    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 
 
>    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1}

Out[7]= {affineWeight[2, finiteWeight[2, {0, 1}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 11], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 3]}

Out[6]= {affineWeight[2, finiteWeight[2, {2, 0}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 10], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {0, 2}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 11], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 2}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 7], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 1, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 6], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 1, 4], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 1, 3], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 1, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 1, 8]}

Out[5]= formalElement[table$117]

makeAffineExtension[makeSimpleRootSystem[B,2]][gradeLimit]

Unassign[ht[2]]



Out[162]= Unassign[2]

Clear::ssym: ht[2] is not a symbol or a string.

Out[160]= 2

Out[159]= 0

Out[158]= 11

Out[157]= affineRootSystem[2, finiteRootSystem[2, 3, 
 
>      {finiteWeight[3, {1, -1, 0}], finiteWeight[3, {0, 1, -1}]}], 
 
>     affineWeight[3, finiteWeight[3, {-1, 0, 1}], 0, 1], 
 
>     {affineWeight[3, finiteWeight[3, {1, -1, 0}], 0, 0], 
 
>      affineWeight[3, finiteWeight[3, {0, 1, -1}], 0, 0]}][gradeLimit]

b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]

b2a[gradeLimit]=11

Out[156]= 11

Out[155]= 10

Out[154]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]

Out[153]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]

Out[152]= 2

Out[147]= {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1]}

Length[pr]

Out[166]= 98

Out[146]= 9

Out[145]= {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1]}

fe=makeFormalElement[pr,Table[1,{Length[pr]}]];

(fe+fe)[multiplicities]

Out[179]= {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
 
>    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
 
>    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
 
>    2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2}

Out[177]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    res$288[#1] + formalElement[table$277][#1] & , 1, 1}

formalElement/:x_formalElement + y_formalElement:=Module[{res},
							 res=makeFormalElement[x[weights],x[multiplicities]];
							 Scan[(res[[1]][#]=res[#]+y[#])&,y[weights]];
							 res];


Out[172]= formalElement[table$285]

Out[171]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}

Out[170]= {affineWeight[2, finiteWeight[2, {-1, 0}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 4], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 7], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {-1, 1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 8], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 10], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 3], 
 
>    affineWeight[2, finiteWeight[2, {-1, 0}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 6], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 9], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 5], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 2], 
 
>    affineWeight[2, finiteWeight[2, {0, -1}], 0, 6]}

Out[169]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}

Out[168]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}

Out[144]= {1, 1, 1, 1, 1, 1, 1, 1, 1}

Out[142]= {1, 1, 1, 1, 1, 1, 1, 1, 1}

Out[136]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}

fe[makeAffineWeight[{1,2},1,1]]

Out[143]= 0

Out[138]= formalElement[table$193][affineWeight[2, finiteWeight[2, {1, 2}], 
 
>     1, 1]]

Out[137]= formalElement[table$221]

Out[133]= Null[multiplicities]

Out[131]= {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
 
>    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}



Message::name: Message name MessageName[fan[(rs_)?rootSystemQ, 
      (subs_)?rootSystemQ], usage] is not of the form symbol::name or
     symbol::name::language.

t=hh/@{a,b,c}

t

Map[hh[#]:=1&,{f,g,p}]

Out[15]= {Null[f], Null[g], Null[p]}

Out[13]= {1, 1, 1}

hh[f]

keys::"usage"="keys[hashtable] gives all the keys in hashtable";
keys = DownValues[#,Sort->False][[All,1,1,1]]&;

DownValues[ff]={HoldPattern[ff[1]] :> 2, HoldPattern[ff[2]] :> 4, HoldPattern[ff[3]] :> 2, HoldPattern[ff[4]] :> 4}

makeHashtable[keys_List,values_List]:=Module[{table},
					     DownValues[table]=((table[#[[1]]] -> #[[2]] &)/@ Transpose[{keys,values}]);
					     table]

ht=makeHashtable[{1,2,3},{4,5,6,7}]

Module[{ht,tt},
       ht=makeHashtable[{"a",2,tt},{1,2,3}];
       Expect["Hashtable test",1,ht["a"]];
       Expect["Hashtable test",3,ht[tt]];
       Expect["Hashtable test",2,ht[2]]]

ht[n_?NumberQ]=0

ht[3]

Out[91]= 0

Out[90]= 2

Out[89]= 0

Most[keys[ht]]

?Sort

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


?ExpandDenominator

values = Function[ht,(ht[#]&)/@keys[ht]]

values[ht]

keys[ht]

?Union

ht[bb]=ht[bb]+1

?NumberQ

keys[ht]

ht[14]=22

ht[12]

Out[122]= 0

Out[121]= 0

Out[120]= 0

Out[119]= 22

Out[118]= 22

Out[117]= {a, 2, tt, bb, (n_)?NumberQ}

NumberQ[expr] gives True if expr is a number, and False otherwise. 

Information::notfound: Symbol Assigned not found.

Out[113]= table$114[bb]

Out[112]= 0

Union[list , list , ...] gives a sorted list of all the distinct elements that
          1      2
     appear in any of the list . Union[list]
                              i
      gives a sorted version of a list, in which all duplicated elements have
      been dropped. 

Out[110]= {a, 2, tt, (n_)?NumberQ}

Out[109]= {1, 2, 3, table$114[(n_)?NumberQ]}

Out[108]= Function[ht, (ht[#1] & ) /@ keys[ht]]

ExpandDenominator[expr] expands out products and powers that appear as
     denominators in expr. 

Expand            ExpandDenominator ExpandFileName    ExpandNumerator
ExpandAll

Expand[expr] expands out products and positive integer powers in expr
    . Expand[expr, patt] leaves unexpanded any parts of expr

       that are free of the pattern patt. 

Out[103]= NestWhileList

Sort[list] sorts the elements of list
     into canonical order. Sort[list, p] sorts using the ordering function p. 


System`MergeDifferences

Attributes[MergeDifferences] = {Protected}

TopologicalSort[g] gives a list of vertices of g
     in topologically sorted order for a directed acyclic graph g.

MenuSortingValue Sort             SortBy           TopologicalSort

Sort   SortBy

Join           Joined         JoinedCurve    JoinedCurveBox JoinForm

Join[list , list , ...]
         1      2
      concatenates lists or other expressions that share the same head.Join[

     list , list , ..., n] joins the objects at level n in each of the list .
         1      2                                                          i

Attributes[Join] = {Flat, OneIdentity, Protected}

Join[list , list , ...]
         1      2
      concatenates lists or other expressions that share the same head.Join[

     list , list , ..., n] joins the objects at level n in each of the list .
         1      2                                                          i

Out[94]= {a, 2, tt}

ht[d]

Out[92]= {a, 2, tt, (n_)?NumberQ}

ht[2]

Out[88]= 2

Out[87]= {a, 2, tt, (n_)?NumberQ}

Out[86]= 0

Out[84]= table$114

                                    
Out[81]= Expect[Hashtable test, 2, 2]

ht=makeHashtable[{"a",2,tt},{1,2,3}]

Expect["Hashtable test",4,ht[2]]

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


Out[80]= Expect[Hashtable test, 4, 2]

ht[tt]

Out[79]= 3

Out[78]= 2

Out[77]= table$114[4]

Out[76]= table$114
                                    
Out[75]= Expect[Hashtable test, 4, 2]

                                    
Out[74]= Expect[Hashtable test, 2, 2]
       


Transpose::nmtx: 
   The first two levels of the one-dimensional list {{1, 2, 3}, {4, 5, 6, 7}}
     cannot be transposed.

DownValues::vlist: 
   Cannot set DownValues[table$97] to 
    Transpose[table$97[{1, 2, 3}] -> {4, 5, 6, 7}]
    ; value must be a list of rules.

Out[69]= table$97

Out[63]= table$95

Out[59]= table$94

Out[52]= table$93

DownValues::vrule: 
   Cannot set DownValues[table$90] to {4, 5, 6}; value contains 4, which is
     not a rule.

Out[50]= table$90

Out[45]= table$89

DownValues[ht]

Out[46]= {}

keys[ht]

Out[64]= {1, 2, 3}

Out[60]= {1, 2, 3}

Out[53]= {1, 2, 3}

ht[1]

ht2=makeHashtable[{a1,a2,a4},{b8,b2,b1}]

keys[ht2]

ht2[a4]

Out[68]= b1

Out[67]= {a1, a2, a4}

Out[66]= table$96

Out[65]= 4

Out[61]= 4

Out[54]= 

DownValues[table$93]

Out[57]= {HoldPattern[table$93[{1, 4}[[1]]]] :> {1, 4}[[2]], 
 
>    HoldPattern[table$93[{2, 5}[[1]]]] :> {2, 5}[[2]], 
 
>    HoldPattern[table$93[{3, 6}[[1]]]] :> {3, 6}[[2]]}

Out[56]= table$93

Out[55]= table$93[1]

Out[48]= table$89[1]

Out[47]= {1, 2, 3}

Out[43]= {1, 2, 3}

Out[42]= {}

Out[41]= table$88[2]

Out[40]= {1, 2, 3}

Out[39]= table$88

(HoldPattern[#[[1]]] :> #[[2]] &)/@ Transpose[{{1,2,3},{4,5,6}}]

?HoldPattern

HoldPattern[expr] is equivalent to expr
     for pattern matching, but maintains expr in an unevaluated form. 

Out[36]= {HoldPattern[{1, 4}[[1]]] :> {1, 4}[[2]], 
 
>    HoldPattern[{2, 5}[[1]]] :> {2, 5}[[2]], 
 
>    HoldPattern[{3, 6}[[1]]] :> {3, 6}[[2]]}

DownValues[ht]

Out[35]= {}

(#[[1]]+#[[2]]&)/@ Transpose[{{1,2,3},{4,5,6}}]

Out[34]= {5, 7, 9}

Out[33]= {{1, 4}, {2, 5}, {3, 6}}

Transpose::tperm: 
   Permutation {4, 5, 6} is longer than the dimensions {3} of the array.

Out[32]= Transpose[{1, 2, 3}, {4, 5, 6}]

keys[ht]

Out[31]= {{1, 4}, {2, 5}, {3, 6}}

Out[30]= table$83[1]

Out[29]= table$83

keys[ff]

HoldPattern[

Out[27]= {1, 2, 3, 4}

Out[25]= keys[ff]

Out[24]= {HoldPattern[ff[1]] :> 2, HoldPattern[ff[2]] :> 4, 
 
>    HoldPattern[ff[3]] :> 2, HoldPattern[ff[4]] :> 4}

Out[23]= {HoldPattern[ff[1]] :> 2, HoldPattern[ff[2]] :> 4}

Out[21]= {HoldPattern[ff[1]] :> 2}

ff[2]=4

Out[22]= 4

Out[20]= 2

Out[19]= {}

Out[17]= {HoldPattern[hh[#1]] :> (1 & )}

Out[16]= hh[f]

Out[14]= hh[f]

Out[12]= hh[1]

Out[11]= {1, 1, 1}

Out[10]= {hh[1], hh[2], hh[3]}

Out[9]= {hh[1], hh[2], hh[3]}

Out[8]= {hh[1], hh[2], hh[3]}

Out[6]= {hh[1], hh[2], hh[3]}

Set::write: Tag Map in hh /@ {1, 2, 3} is Protected.

Out[5]= {1, 2, 3}

Set::write: Tag Map in h /@ {1, 2, 3} is Protected.

Out[4]= {1, 2, 3}

b

Out[3]= 2

Out[2]= 1

Out[1]= {1, 2, 3}

b2=makeSimpleRootSystem[B,2]

makeFiniteRootSystem[{{1,-1,0,0,0,0},{0,1,0,0,0,0},{0,0,1,-1,0,0},{0,0,0,1,-1,0},{0,0,0,0,1,-1}}]

Module[{b2=makeSimpleRootSystem[B,2],a3=makeSimpleRootSystem[A,3]},
       Expect["Direct sum of finite-dimensional Lie algebras",True,
	      b2+a3==makeFiniteRootSystem[{{1,-1,0,0,0,0},
					   {0,1,0,0,0,0},
					   {0,0,1,-1,0,0},
					   {0,0,0,1,-1,0},
					   {0,0,0,0,1,-1}}]]]

Module[{b2=makeSimpleRootSystem[B,2],a3=makeSimpleRootSystem[A,3]},
       Expect["Direct sum of affine Lie algebras",True,
	      makeAffineExtension[b2] + makeAffineExtension[a3]==makeAffineExtension[makeFiniteRootSystem[{{1,-1,0,0,0,0},
					   {0,1,0,0,0,0},
					   {0,0,1,-1,0,0},
					   {0,0,0,1,-1,0},
					   {0,0,0,0,1,-1}}]]]]

Expect["checkGrade finite dimensional Lie algebra", True, checkGrade[makeAffineExtension[makeSimpleRootSystem[B,2]]][makeFiniteWeight[{2,1}]]]

Expect["checkGrade affine Lie algebra", True, checkGrade[makeAffineExtension[makeSimpleRootSystem[B,2]]][makeAffineWeight[{2,1},1,10]]]

b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]

b2a[gradeLimit]=15

b2a

Out[30]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]

Out[29]= 15

Out[28]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]

Module[{b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]},
       b2a[gradeLimit]=15;
       checkGrade[b2a][makeAffineWeight[{2,1},1,10]]]

Module[{b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]},
       b2a[gradeLimit]=15;
       Expect["checkGrade affine Lie algebra", True, checkGrade[b2a][makeAffineWeight[{2,1},1,10]]]]

Expect["weight predicate on finite weights", False, weightQ[{1,2,3}]]

Expect["weight predicate on finite weights", False, weightQ[makeFiniteWeight[1,2,3]]]

Expect["weight predicate on affine weights", False, weightQ[makeAffineWeight[1,2,3]]]

b2a[gradeLimit]=1

Out[44]= 1

positiveRoots[b2a]

Out[45]= {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 0, 1]}

Module[{b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]},
       b2a[gradeLimit]=1;
       Expect["Positive roots of affine B2",True, positiveRoots[b2a]=={makeAffineWeight[{-1, -1}, 0, 1], 
								       makeAffineWeight[{1, -1}, 0, 0], 
								       makeAffineWeight[{0, 1}, 0, 0], 
								       makeAffineWeight[{1, 1}, 0, 0], 
								       makeAffineWeight[{1, 0}, 0, 0], 
								       makeAffineWeight[{0, 0}, 0, 0], 
								       makeAffineWeight[{0, 0}, 0, 1], 
								       makeAffineWeight[{0, 0}, 0, 0], 
								       makeAffineWeight[{0, 0}, 0, 1]}]]

b2a=makeAffineExtension[makeSimpleRootSystem[B,2]]

Out[47]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]


b2a[gradeLimit]=3;

?Product

?Times

x * y * z, x × y × z or x y z represents a product of terms. 

                              Mutliplication by numbers is defined for the
                                elements of weight space of affine and
                                finite-dimensional Lie algebras. 
                                 
                              For example
                                makeFiniteWeight[{1,2,3}]*2==makeFiniteWeight[
                                {2,4,6}]]

                                               i
Product[f, {i, i   }] evaluates the product     max    f
                max                         \[Product]
                                              i = 1
    . Product[f, {i, i   , i   }]
                      min   max
      starts with i = i   . Product[f, {i, i   , i   , di}]
                       min                  min   max
        uses steps di. Product[f, {i, {i , i , ...}}]
                                        1   2
          uses successive values i
                                  1
          , i , ....Product[f, {i, i   , i   }, {j, j   , j   }, ...]
             2                      min   max        min   max
                                                i          j
                                                 max        max
              evaluates the multiple product \[Product] \[Product] ... f
                                              i = i      j = j
                                                   min        min
              . Product[f, i] gives the indefinite product \[Product] f.
                                                               i

Product             ProductDistribution ProductLog

1-Exp[{1,2,3}]

                      2       3
Out[56]= {1 - E, 1 - E , 1 - E }

Length[positiveRoots[b2a]]
Out[10]= pr

pr

[Calculating...]

Expand[Times @@ (1-Exp[positiveRoots[b2a]])]

Out[6]= $Aborted

?Expand

Expand[(1-a)(1-b)(1-c)]

Out[3]= 1 - a - b + a b - c + a c + b c - a b c

Expand[expr] expands out products and positive integer powers in expr
    . Expand[expr, patt] leaves unexpanded any parts of expr

       that are free of the pattern patt. 

               affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1]
Out[57]= {1 - E                                                , 
 
          affineWeight[2, finiteWeight[2, {1, -1}], 0, 0]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {1, 1}], 0, 0]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {1, 0}], 0, 0]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {-1, 1}], 0, 1]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {-1, 0}], 0, 1]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {1, -1}], 0, 1]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {0, -1}], 0, 1]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {-1, -1}], 0, 2]
>    1 - E                                                , 
 
          affineWeight[2, finiteWeight[2, {1, 1}], 0, 1]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 1}], 0, 1]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {-1, 1}], 0, 2]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {1, 0}], 0, 1]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {1, -1}], 0, 2]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {-1, 0}], 0, 2]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {1, 1}], 0, 2]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, -1}], 0, 2]
>    1 - E                                               , 
 
          affineWeight[2, finiteWeight[2, {0, 1}], 0, 2]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {1, 0}], 0, 2]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 0]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 1]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 3]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 0]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 1]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 2]
>    1 - E                                              , 
 
          affineWeight[2, finiteWeight[2, {0, 0}], 0, 3]
>    1 - E                                              }

checkGrade[makeAffineExtension[makeSimpleRootSystem[B,2]]][makeAffineWeight[{2,1},1,20]]

Out[37]= False

Out[36]= False

Out[35]= True

Out[34]= True

checkGrade[rs_?rootSystemQ][w_affineWeight]:=(Abs[w[grade]]<rs[gradeLimit]) /. rs[gradeLimit]->10;

Out[32]= True

Out[31]= False
                  
Out[27]= False

                  
Out[26]= False

                  
Out[25]= False

       Expect["checkGrade affine Lie algebra", True, checkGrade[b2a][makeAffineWeight[{2,1},1,10]]]]

                  
Throw::nocatch: 
   Uncaught Throw[checkGrade affine Lie algebra: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[24]= Hold[Throw[checkGrade affine Lie algebra: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]

makeAffineExtension[b2]

Out[11]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]

a3=makeSimpleRootSystem[A,3]

Out[13]= finiteRootSystem[3, 4, 
 
>    {finiteWeight[4, {1, -1, 0, 0}], finiteWeight[4, {0, 1, -1, 0}], 
 
>     finiteWeight[4, {0, 0, 1, -1}]}]

makeAffineExtension[a3]

Out[14]= affineRootSystem[3, finiteRootSystem[3, 4, 
 
>     {finiteWeight[4, {1, -1, 0, 0}], finiteWeight[4, {0, 1, -1, 0}], 
 
>      finiteWeight[4, {0, 0, 1, -1}]}], 
 
>    affineWeight[4, finiteWeight[4, {-1, 0, 0, 1}], 0, 1], 
 
>    {affineWeight[4, finiteWeight[4, {1, -1, 0, 0}], 0, 0], 
 
>     affineWeight[4, finiteWeight[4, {0, 1, -1, 0}], 0, 0], 
 
>     affineWeight[4, finiteWeight[4, {0, 0, 1, -1}], 0, 0]}]


affineRootSystem/:x_affineRootSystem+y_affineRootSystem:=makeAffineExtension[x[finiteRootSystem]+y[finiteRootSystem]];

Out[11]+Out[14]

Out[18]= affineRootSystem[5, finiteRootSystem[5, 6, 
 
>     {finiteWeight[6, {1, -1, 0, 0, 0, 0}], 
 
>      finiteWeight[6, {0, 1, 0, 0, 0, 0}], 
 
>      finiteWeight[6, {0, 0, 1, -1, 0, 0}], 
 
>      finiteWeight[6, {0, 0, 0, 1, -1, 0}], 
 
>      finiteWeight[6, {0, 0, 0, 0, 1, -1}]}], 
 
>    affineWeight[6, finiteWeight[6, {0, 0, -1, 0, 0, 1}], 0, 1], 
 
>    {affineWeight[6, finiteWeight[6, {1, -1, 0, 0, 0, 0}], 0, 0], 
 
>     affineWeight[6, finiteWeight[6, {0, 1, 0, 0, 0, 0}], 0, 0], 
 
>     affineWeight[6, finiteWeight[6, {0, 0, 1, -1, 0, 0}], 0, 0], 
 
>     affineWeight[6, finiteWeight[6, {0, 0, 0, 1, -1, 0}], 0, 0], 
 
>     affineWeight[6, finiteWeight[6, {0, 0, 0, 0, 1, -1}], 0, 0]}]

Out[16]= affineRootSystem[2, finiteRootSystem[2, 2, 
 
>      {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>     affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>     {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>      affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}] + 
 
>    affineRootSystem[3, finiteRootSystem[3, 4, 
 
>      {finiteWeight[4, {1, -1, 0, 0}], finiteWeight[4, {0, 1, -1, 0}], 
 
>       finiteWeight[4, {0, 0, 1, -1}]}], 
 
>     affineWeight[4, finiteWeight[4, {-1, 0, 0, 1}], 0, 1], 
 
>     {affineWeight[4, finiteWeight[4, {1, -1, 0, 0}], 0, 0], 
 
>      affineWeight[4, finiteWeight[4, {0, 1, -1, 0}], 0, 0], 
 
>      affineWeight[4, finiteWeight[4, {0, 0, 1, -1}], 0, 0]}]

Out[12]= makeAffineExtension[a3]
                                                      
Out[10]= If[True != (affineRootSystem[2, 
 
>         finiteRootSystem[2, 2, 
 
>          {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>         affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>         {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>          affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}] + 
 
>        affineRootSystem[3, finiteRootSystem[3, 4, 
 
>          {finiteWeight[4, {1, -1, 0, 0}], finiteWeight[4, {0, 1, -1, 0}], 
 
>           finiteWeight[4, {0, 0, 1, -1}]}], 
 
>         affineWeight[4, finiteWeight[4, {-1, 0, 0, 1}], 0, 1], 
 
>         {affineWeight[4, finiteWeight[4, {1, -1, 0, 0}], 0, 0], 
 
>          affineWeight[4, finiteWeight[4, {0, 1, -1, 0}], 0, 0], 
 
>          affineWeight[4, finiteWeight[4, {0, 0, 1, -1}], 0, 0]}] == 
 
>       affineRootSystem[5, finiteRootSystem[5, 6, 
 
>         {finiteWeight[6, {1, -1, 0, 0, 0, 0}], 
 
>          finiteWeight[6, {0, 1, 0, 0, 0, 0}], 
 
>          finiteWeight[6, {0, 0, 1, -1, 0, 0}], 
 
>          finiteWeight[6, {0, 0, 0, 1, -1, 0}], 
 
>          finiteWeight[6, {0, 0, 0, 0, 1, -1}]}], 
 
>        affineWeight[6, finiteWeight[6, {0, 0, -1, 0, 0, 1}], 0, 1], 
 
>        {affineWeight[6, finiteWeight[6, {1, -1, 0, 0, 0, 0}], 0, 0], 
 
>         affineWeight[6, finiteWeight[6, {0, 1, 0, 0, 0, 0}], 0, 0], 
 
>         affineWeight[6, finiteWeight[6, {0, 0, 1, -1, 0, 0}], 0, 0], 
 
>         affineWeight[6, finiteWeight[6, {0, 0, 0, 1, -1, 0}], 0, 0], 
 
>         affineWeight[6, finiteWeight[6, {0, 0, 0, 0, 1, -1}], 0, 0]}]), 
 
>    Throw[ToString[Direct sum of affine Lie algebras]<>: GOT UNEXPECTED\
 
>       VALUE <>ToString[affineRootSystem[2, 
 
>          finiteRootSystem[2, 2, 
 
>           {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>          affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>          {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>           affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}] + 
 
>         affineRootSystem[3, finiteRootSystem[3, 4, 
 
>           {finiteWeight[4, {1, -1, 0, 0}], finiteWeight[4, {0, 1, -1, 0}], 
 
>            finiteWeight[4, {0, 0, 1, -1}]}], 
 
>          affineWeight[4, finiteWeight[4, {-1, 0, 0, 1}], 0, 1], 
 
>          {affineWeight[4, finiteWeight[4, {1, -1, 0, 0}], 0, 0], 
 
>           affineWeight[4, finiteWeight[4, {0, 1, -1, 0}], 0, 0], 
 
>           affineWeight[4, finiteWeight[4, {0, 0, 1, -1}], 0, 0]}] == 
 
>        affineRootSystem[5, finiteRootSystem[5, 6, 
 
>          {finiteWeight[6, {1, -1, 0, 0, 0, 0}], 
 
>           finiteWeight[6, {0, 1, 0, 0, 0, 0}], 
 
>           finiteWeight[6, {0, 0, 1, -1, 0, 0}], 
 
>           finiteWeight[6, {0, 0, 0, 1, -1, 0}], 
 
>           finiteWeight[6, {0, 0, 0, 0, 1, -1}]}], 
 
>         affineWeight[6, finiteWeight[6, {0, 0, -1, 0, 0, 1}], 0, 1], 
 
>         {affineWeight[6, finiteWeight[6, {1, -1, 0, 0, 0, 0}], 0, 0], 
 
>          affineWeight[6, finiteWeight[6, {0, 1, 0, 0, 0, 0}], 0, 0], 
 
>          affineWeight[6, finiteWeight[6, {0, 0, 1, -1, 0, 0}], 0, 0], 
 
>          affineWeight[6, finiteWeight[6, {0, 0, 0, 1, -1, 0}], 0, 0], 
 
>          affineWeight[6, finiteWeight[6, {0, 0, 0, 0, 1, -1}], 0, 0]}]]<>
 
>       INSTEAD OF <>ToString[True], assertion exception]]

makeSimpleRootSystem[A,3]

Out[7]= finiteRootSystem[3, 4, {finiteWeight[4, {1, -1, 0, 0}], 
 
>     finiteWeight[4, {0, 1, -1, 0}], finiteWeight[4, {0, 0, 1, -1}]}]

Out[6]= makeFiniteRootSystem[A, 3]

Out[5]= finiteRootSystem[2, 2, {finiteWeight[2, {1, -1}], 
 
>     finiteWeight[2, {0, 1}]}]

Out[4]= makeFiniteRootSystem[{1, -1}, {0, 1}]

Out[3]= finiteRootSystem[2, 2, {finiteWeight[2, {1, -1}], 
 
>     finiteWeight[2, {0, 1}]}]

Message::name: Message name MessageName[fan[(rs_)?rootSystemQ, 
      (subs_)?rootSystemQ], usage] is not of the form symbol::name or
     symbol::name::language.

Syntax::sntue: Unexpected end of file (probably unfinished expression) 
     (line 619 of "datastructures.m").

b2=makeSimpleRootSystem[B,2]

Exp[b2[simpleRoot][1]]*Exp[b2[simpleRoot][1]]

         finiteWeight[2, {2, -2}]
Out[8]= E

           finiteWeight[2, {1, -1}]
Out[7]= 2 E

         finiteWeight[2, {0, 1}]    finiteWeight[2, {1, -1}]
Out[6]= E                        + E

         finiteWeight[2, {1, 0}]
Out[5]= E

         finiteWeight[2, {1, -1}]
Out[4]= E

Out[3]= finiteWeight[2, {1, -1}]

Out[2]= finiteRootSystem[2, 2, {finiteWeight[2, {1, -1}], 
 
>     finiteWeight[2, {0, 1}]}]

b2Out[45]= finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>     finiteWeight[2, {0, 1}]}]

Out[46]= {}

orthogonalSubsystem[b2,makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]]

marks::"usage"="comarks[rs_affineRootSystem] returns comarks of affine Lie algebra";
comarks[rs_affineRootSystem]:=marks[rs]*Map[#.#/2&,rs[simpleRoots]]

fundamentalWeights[b2].{3,2}                                                    1  1
Out[74]= {finiteWeight[2, {1, 0}], finiteWeight[2, {-, -}]}
                                                    2  2

fundamentalWeights[b2]

makeWeight[rs_?rootSystemQ][labels__Integer]:=fundamentalWeights[rs].{labels}

makeWeight[makeAffineExtension[b2]][0,0,1]

weight[rs_?rootSystemQ][labels__Integer]:=fundamentalWeights[rs].{labels}

ww=weight[A2][1,0]




Clear[A2]

A2=makeSimpleRootSystem[A,2]

[Calculating...]

[Calculating...]

Out[87]= finiteRootSystem[2, {finiteWeight[3, {1, -1, 0}], 
 
>     finiteWeight[3, {0, 1, -1}]}]

Out[86]= weight[A2][1, 0]

[Calculating...]

                                          1  1
Out[83]= affineWeight[2, finiteWeight[2, {-, -}], 1, 0]
                                          2  2

Out[82]= affineWeight[2, finiteWeight[2, {0, 0}], 1, 0]

[Calculating...]

                          7  5
Out[78]= finiteWeight[2, {-, -}]
                          2  2

                          3  1
Out[77]= finiteWeight[2, {-, -}]
                          2  2

Out[76]= finiteWeight[2, {1, 0}]

[Calculating...]

[Calculating...]

orthogonalSubsystem[rs_?rootSystemQ,subs_?rootSystemQ]:=Cases[positiveRoots[rs], z_ /; Or @@ (z.#==0 &  /@ subs[simpleRoots])]

Out[67]= {finiteWeight[2, {1, -1}]}

Out[63]= {finiteWeight[2, {1, -1}]}

Out[61]= {finiteWeight[2, {1, -1}]}

Out[59]= {True, False, False, False}

Out[57]= {True, False, False, False}

Out[55]= {{True}, {False}, {False}, {False}}

Out[53]= {{True}, {False}, {False}, {False}}

Out[51]= {}

Out[49]= {}

Out[47]= {}

111

Outer::heads: Heads List and Times at positions 3 and 2
     are expected to be the same.

Throw::nocatch: 
   Uncaught Throw[Weights of [2,1] module of B2: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[1]= Hold[Throw[Weights of [2,1] module of B2: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]



Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["Partial orbit test", True, partialOrbit[b2][{rho[b2]}]==
	      {{makeFiniteWeight[{3/2,1/2}]},{makeFiniteWeight[{1/2,3/2}],makeFiniteWeight[{3/2,-1/2}]},
	       {makeFiniteWeight[{-1/2,3/2}],makeFiniteWeight[{1/2,-3/2}]},
	       {makeFiniteWeight[{-3/2,1/2}],makeFiniteWeight[{-1/2,-3/2}]},
	       {makeFiniteWeight[{-3/2,-1/2}]}}]]

Throw::nocatch: 
   Uncaught Throw[Partial orbit test: GOT UNEXPECTED VALUE False INSTEAD OF
      True, assertion exception] returned to top level.

Out[1]= Hold[Throw[Partial orbit test: GOT UNEXPECTED VALUE False INSTEAD OF\
 
>      True, assertion exception]]

Throw::nocatch: 
   Uncaught Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[16]= Hold[Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]

b2=makeSimpleRootSystem[B,2]

Out[6]= finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>     finiteWeight[2, {0, 1}]}]

b2+b2

[Calculating...]

Out[42]= 2 finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>      finiteWeight[2, {0, 1}]}]

partialOrbit[b2][{rho[b2]}]

                            3  1
Out[16]= {{finiteWeight[2, {-, -}]}, 
                            2  2
 
                       1  3                     3    1
>    {finiteWeight[2, {-, -}], finiteWeight[2, {-, -(-)}]}, 
                       2  2                     2    2
 
                         1   3                     1    3
>    {finiteWeight[2, {-(-), -}], finiteWeight[2, {-, -(-)}]}, 
                         2   2                     2    2
 
                         3   1                       1     3
>    {finiteWeight[2, {-(-), -}], finiteWeight[2, {-(-), -(-)}]}, 
                         2   2                       2     2
 
                         3     1
>    {finiteWeight[2, {-(-), -(-)}]}}
                         2     2

checkGrade[b2][rho[b2]]

Out[12]= checkGrade[finiteRootSystem[2, 
 
>      {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}]][finiteWeight[2, 
 
       3  1
>     {-, -}]]
       2  2

                            3  1
Out[11]= {{finiteWeight[2, {-, -}]}}
                            2  2

                           3  1
Out[7]= {{finiteWeight[2, {-, -}]}}
                           2  2

weightSystem[b2][makeFiniteWeight[{2,1}]]

Out[19]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}}

Out[6]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}, {}}


Module[{b2=makeSimpleRootSystem[B,2]},
       Expect["Weights of [1,1] module of B2",True, weightSystem[b2][makeFiniteWeight[{2,1}]]==
	      {{makeFiniteWeight[{2,1}]},
	      {makeFiniteWeight[{1,0}],makeFiniteWeight[{1,1}],makeFiniteWeight[{2,0}]},
	      {makeFiniteWeight[{0,0}]}}]]


weightSystem[b2][makeFiniteWeight[{2,1}]]

Out[15]= {{finiteWeight[2, {2, 1}]}}
                                    
Throw::nocatch: 
   Uncaught Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[14]= Hold[Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]

Module[{b2=makeSimpleRootSystem[B,2],fm},
       fm=freudenthalMultiplicities[b2][makeFiniteWeight[{2,1}]];
       Expect["Mutliplicites of [2,1] B2 representation",{1, 1, 2, 3, 3},
	      fm[#]&/@keys[fm]]]

freudenthalMultiplicities[makeAffineExtension[b2]][makeAffineWeight[{2,1},1,0]]


orbitWithEps[rs_?rootSystemQ][weight_?weightQ]:=Flatten[Most[MapIndexed[Function[{x,i},Map[{#,(-1)^(i[[1]]+1)}&,x]],orbit[rs][weight]]],1];

orbitWithEps[makeAffineExtension[b2]][makeAffineWeight[{2,1},1,0]]

toFundamentalChamber[makeAffineExtension[b2]][makeAffineWeight[{2,1},1,0]]

Out[4]= affineWeight[2, finiteWeight[2, {1, 0}], 1, 2]

Out[3]= {{affineWeight[2, finiteWeight[2, {1, 0}], 1, 2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, 1}], 1, 2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, -1}], 1, 2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, 0}], 1, 2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, 1}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, 2}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, -1}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, 2}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, -2}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, 1}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, -2}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, 0}], 1, -2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, -1}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, 3}], 1, -2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, 2}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, -3}], 1, -2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, 3}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, -2}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-3, 0}], 1, -2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, 3}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, -3}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {4, 1}], 1, -6], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-3, 2}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, -3}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, 4}], 1, -6], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {4, -1}], 1, -6], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-3, -2}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, 4}], 1, -6], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, -4}], 1, -6], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-4, 1}], 1, -6], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, -4}], 1, -6], 1}}


racahMultiplicities[rs_?rootSystemQ][highestWeight_?weightQ]:=
    Module[{rh=rho[rs],weights,mults,c,insideQ,
	    fan,
	    toFC=toFundamentalChamber[rs]},
	   fan=Map[{rh-#[[1]],#[[2]]}&,Rest[orbitWithEps[rs][rh]]];
	   Print[fan];
	   weights=Sort[ Rest[Flatten[weightSystem[rs][highestWeight]]], #1.rh>#2.rh&];
	   mults[highestWeight]=1;
	   insideQ:=IntegerQ[mults[toFC[#]]]&;
	   Scan[Function[v,
			 mults[v]=
			 Plus@@(fan /. {x_?weightQ,e_Integer}:> If[insideQ[v+x],-e*mults[toFC[v+x]],0])],
		weights];
	   mults]

racahMultiplicities[makeAffineExtension[b2]][makeAffineWeight[{2,1},1,0]]

[Calculating...]

[Calculating...]

Out[7]= {{affineWeight[2, finiteWeight[2, {1, 0}], 1, 2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, 1}], 1, 2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, -1}], 1, 2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, 0}], 1, 2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, 1}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, 2}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, -1}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, 2}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, -2}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, 1}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, -2}], 1, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, 0}], 1, -2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, -1}], 1, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, 3}], 1, -2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, 2}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, -3}], 1, -2], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, 3}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, -2}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-3, 0}], 1, -2], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, 3}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, -3}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {4, 1}], 1, -6], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-3, 2}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-2, -3}], 1, -4], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, 4}], 1, -6], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {4, -1}], 1, -6], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-3, -2}], 1, -4], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, 4}], 1, -6], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, -4}], 1, -6], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {-4, 1}], 1, -6], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {-1, -4}], 1, -6], 1}}

Out[3]= orbitWithEps[makeAffineExtension[b2]][affineWeight[2, 
 
>     finiteWeight[2, {2, 1}], 1, 0]]

Out[24]= $Aborted

Out[22]= freudenthalMultiplicities[affineRootSystem[2, 
 
>      finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>        finiteWeight[2, {0, 1}]}], 
 
>      affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>      {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>       affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]][affineWeight[2, 
 
>     finiteWeight[2, {2, 1}], 1, 0]]

Out[20]= weightSystem[affineRootSystem[2, 
 
>      finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>        finiteWeight[2, {0, 1}]}], 
 
>      affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>      {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>       affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]][{affineWeight[2, 
 
>      finiteWeight[2, {2, 1}], 1, 0]}]



Out[24]= $Aborted

Out[23]= freudenthalMultiplicities[affineRootSystem[2, 
 
>      finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>        finiteWeight[2, {0, 1}]}], 
 
>      affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>      {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>       affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]][affineWeight[2, 
 
>     finiteWeight[2, {2, 1}], 1, 0]]

Out[22]= affineRootSystem[2, finiteRootSystem[2, 
 
>     {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    {affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>     affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}]
       
fm=freudenthalMultiplicities[b2][makeFiniteWeight[{2,1}]]

fm[#]&/@keys[fm]

Out[20]= {1, 1, 2, 3, 3}

Out[19]= {finiteWeight[2, {2, 1}], finiteWeight[2, {2, 0}], 
 
>    finiteWeight[2, {1, 1}], finiteWeight[2, {1, 0}], 
 
>    finiteWeight[2, {0, 0}]}

Out[18]= mults$152

{{makeFiniteWeight[{2,1}]},
	      {makeFiniteWeight[{1,0}],makeFiniteWeight[{1,1}],makeFiniteWeight[{2,0}]},
	      {makeFiniteWeight[{0,0}]}}

                  
Out[15]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}}

weightSystem[b2][makeFiniteWeight[{2,1}]]

Out[14]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}, {}}

                  
Out[13]= False

                  
Out[12]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}}

weightSystem[b2][makeFiniteWeight[{2,1}]]

Out[11]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}, {}}
                                    
Throw::nocatch: 
   Uncaught Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[10]= Hold[Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]

                                        
Throw::nocatch: 
   Uncaught Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[7]= Hold[Throw[Weights of [1,1] module of B2: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]

Out[5]= {{finiteWeight[2, {2, 2}]}, 
 
>    {finiteWeight[2, {1, 1}], finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {0, 0}], finiteWeight[2, {1, 0}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {1, 1}]}, 
 
>    {finiteWeight[2, {0, 0}], finiteWeight[2, {1, 0}]}, {}}

                           3  1                       1  1
Out[4]= {{finiteWeight[2, {-, -}]}, {finiteWeight[2, {-, -}]}, {}}
                           2  2                       2  2

?affineWeight

weightQ[x_finiteWeight]=True;
weightQ[x_affineWeight]=True;
weightQ[_]=False;

weightQ[1]

coroot[x_?weightQ]:=2*x/(x.x);

coroot[makeAffineWeight[{1,0},1,0]]

Expect["Co root of affine [1,0]", True, coroot[makeAffineWeight[{1,0},1,0]]==makeAffineWeight[{2,0},2,0]]

cartanMatrix[rs_?rootSystemQ]:=Transpose[Outer[Dot,rs[simpleRoots],coroot/@rs[simpleRoots]]];

cartanMatrix[makeAffineExtension[makeSimpleRootSystem[A,4]]]//MatrixForm

rs=makeAffineExtension[makeSimpleRootSystem[A,4]]

sr=makeAffineExtension[makeSimpleRootSystem[A,4]][simpleRoots]

rh=rho[makeAffineExtension[makeSimpleRootSystem[A,4]]]

Out[98]= affineWeight[5, finiteWeight[5, {2, 1, 0, -1, -2}], 5, 0]

FullForm[weylGroupElement[0,1,2,1,0,1][rs]]

y=aaa[1,2,3]

Append[y,4]

Out[118]= aaa[1, 2, 3, 4]

Out[117]= append[aaa[1, 2, 3], 4]

Out[116]= aaa[1, 2, 3]


Out[113]//FullForm= 
 
>   affineWeight[5, finiteWeight[5, List[3, 2, 0, -1, -4]], 5, -2]

Range[0,rs[rank]]

Out[112]= {0, 1, 2, 3, 4}

Out[111]= Rnage[0, 4]

Out[110]= 4

modorbit[rs_affineRootSystem][weight__affineWeight,gradelimit_?NumberQ]:=
			   NestWhileList[
			       Function[x,
					Union[Flatten[
					    Map[Function[y,
							 Map[Append[y,#]&,
							     Cases[Range[0,rs[rank]],z_ /; And[rs[simpleRoot][z].y[weight]>0,Abs[reflection[rs[simpleRoots][z]][y[weight]][grade]]<gradelimit]]
							    ]],
						x]],
					      SameTest->(#1==#2&)]],
			       weylGroupElement[],
			       #=!={}&]

modorbit[rs][rh,1]

orbit[rs][{rh},2]

rh

Out[127]= affineWeight[5, finiteWeight[5, {2, 1, 0, -1, -2}], 5, 0]


Out[125]= {{affineWeight[5, finiteWeight[5, {2, 1, 0, -1, -2}], 5, 0]}, 
 
>    {affineWeight[5, finiteWeight[5, {1, 2, 0, -1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 0, 1, -1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 1, -1, 0, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 1, 0, -2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {3, 1, 0, -1, -3}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {0, 2, 1, -1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, 2, -1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 2, -1, 0, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 2, 0, -2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 3, 0, -1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {2, -1, 1, 0, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 0, -1, 1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 0, 1, -2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 1, -2, 0, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 1, -1, -2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {3, 0, 1, -1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 1, -1, 0, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 1, 0, -3, -1}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-1, 2, 1, 0, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, 2, -1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 2, -1, 1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 2, 1, -2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 3, 1, -1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 2, 0, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -1, 2, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, 2, -2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, 3, -1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 2, -2, 0, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 2, -1, -2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 3, -1, 0, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 3, 0, -3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {2, -2, 1, 0, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -1, 0, 1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -1, 1, -2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 0, -2, 1, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 0, -1, -2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 1, -2, -1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {3, -1, 1, 0, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 0, -1, 1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 0, 1, -3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 1, -3, 0, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 1, -1, -3, 0}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-2, 2, 1, 0, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 2, 0, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 2, 0, 1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 2, 1, -2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 3, 1, 0, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 2, 1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -1, 2, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, 2, -2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, 3, -1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 2, -2, 1, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 2, -1, -2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 3, -1, 1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 3, 1, -3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -2, 2, 0, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 0, 2, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 2, -2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 3, 0, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -2, 2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -1, -2, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -1, 3, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, 3, -3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 2, -2, -1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 3, -3, 0, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 3, -1, -3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {2, -2, 0, 1, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -2, 1, -1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -1, -2, 1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -1, 0, -2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, 0, -2, -1, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {3, -3, 1, 0, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -1, 0, 1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -1, 1, -3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 0, -3, 1, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 0, -1, -3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 1, -3, -1, 0}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, 3, 1, 0, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-2, 1, 2, 0, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 2, 0, 1, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 2, 1, -1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 2, 1, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 0, 2, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 2, -2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 3, 0, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 2, -2, 1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 2, 0, -2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 3, 0, 1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 3, 1, -3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -2, 2, 1, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 1, 2, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 2, -2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 3, 1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -2, 2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -1, -2, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -1, 3, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, 3, -3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 2, -2, -1, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 3, -3, 1, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 3, -1, -3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -3, 3, 0, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -2, 0, 2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -2, 2, -1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, -2, 2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 0, -2, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 0, 3, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 3, -3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -3, 3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -2, -1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -1, -3, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 3, -3, -1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {2, -2, -1, 1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -2, 0, -1, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {2, -1, -2, 0, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {3, -3, 0, 1, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -3, 1, -1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -1, -3, 1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -1, 0, -3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, 0, -3, -1, 1}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, 1, 3, 0, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 3, 0, 1, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 3, 1, -1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-2, 0, 2, 1, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 1, 0, 2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 1, 2, -1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 2, -1, 1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 2, 0, -1, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -2, 2, 1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 1, 2, -2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 2, -2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 3, 1, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, -2, 2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 0, -2, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 0, 3, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 3, -3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 2, -2, 0, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 3, -3, 1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 3, 0, -3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -3, 3, 1, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -2, 1, 2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -2, 2, -1, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, -2, 2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 1, -2, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 1, 3, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 3, -3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -3, 3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -2, -1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -1, -3, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 3, -3, -1, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -3, 0, 3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -3, 3, -1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -2, -1, 2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -2, 0, -1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, -3, 3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, -2, 0, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, 0, -3, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, 0, -3, -1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {2, -2, -1, 0, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {3, -3, -1, 1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -3, 0, -1, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -1, -3, 0, 1}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, 0, 3, 1, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 1, 0, 3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 1, 3, -1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 3, -1, 1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 3, 0, -1, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-2, -1, 2, 1, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 0, 1, 2, -1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 0, 2, -1, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 1, -1, 2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 1, 0, -1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 2, -1, 0, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -3, 3, 1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, -2, 1, 2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -2, 2, 0, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, -2, 2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 1, -2, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 1, 3, -3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 3, -3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, -3, 3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, -2, 0, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, 0, -3, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 3, -3, 0, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -3, 1, 3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -3, 3, -1, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -2, -1, 2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -2, 1, -1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, -3, 3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, -2, 1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, 1, -3, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -3, -1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -3, -1, 3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -3, 0, -1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -2, -1, 0, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {1, -1, -3, 0, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {3, -3, -1, 0, 1}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, -1, 3, 1, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 0, 1, 3, -1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 0, 3, -1, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 1, -1, 3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 1, 0, -1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 3, -1, 0, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-2, -1, 1, 2, 0}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, -1, 2, 0, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 0, -1, 2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 0, 1, -1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 1, -1, 0, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -3, 1, 3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, -3, 3, 0, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, -2, 0, 2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -2, 1, 0, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, -3, 3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, -2, 1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, 1, -3, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, 1, -3, 0, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -3, -1, 3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -3, 1, -1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -2, -1, 1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, -1, -3, 1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {1, -3, -1, 0, 3}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, -1, 1, 3, 0}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, -1, 3, 0, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 0, -1, 3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 0, 1, -1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 1, -1, 0, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-2, -1, 0, 2, 1}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, -1, 1, 0, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-2, 0, -1, 1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -3, 0, 3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, -3, 1, 0, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-1, -2, 0, 1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, 0, -3, 1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {0, -3, -1, 1, 3}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, -1, 0, 3, 1}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, -1, 1, 0, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-3, 0, -1, 1, 3}], 5, -1], 
 
>     affineWeight[5, finiteWeight[5, {-2, -1, 0, 1, 2}], 5, 0], 
 
>     affineWeight[5, finiteWeight[5, {-1, -3, 0, 1, 3}], 5, -1]}, 
 
>    {affineWeight[5, finiteWeight[5, {-3, -1, 0, 1, 3}], 5, -1]}, {}}

[Calculating...]

Out[109]= affineWeight[5, finiteWeight[5, {3, 2, 0, -1, -4}], 5, -2]

Out[108]= affineWeight[5, finiteWeight[5, {3, 2, 0, -1, -4}], 5, -2]

Out[107]= affineWeight[5, finiteWeight[5, {1, 2, 0, -1, -2}], 5, 0]

Out[106]= affineWeight[5, finiteWeight[5, {2, 1, 0, -1, -2}], 5, 0]

Out[105]= affineWeight[5, finiteWeight[5, {3, 1, 0, -1, -3}], 5, -1]

Out[104]= affineWeight[5, finiteWeight[5, {2, 1, 0, -1, -2}], 5, 0]

Out[103]= 0

Out[102]= -1

Out[101]= affineWeight[5, finiteWeight[5, {3, 1, 0, -1, -3}], 5, -1]

Out[97]= {affineWeight[5, finiteWeight[5, {-1, 0, 0, 0, 1}], 0, 1], 
 
>    affineWeight[5, finiteWeight[5, {1, -1, 0, 0, 0}], 0, 0], 
 
>    affineWeight[5, finiteWeight[5, {0, 1, -1, 0, 0}], 0, 0], 
 
>    affineWeight[5, finiteWeight[5, {0, 0, 1, -1, 0}], 0, 0], 
 
>    affineWeight[5, finiteWeight[5, {0, 0, 0, 1, -1}], 0, 0]}


Out[96]= affineRootSystem[4, finiteRootSystem[4, 
 
>     {finiteWeight[5, {1, -1, 0, 0, 0}], finiteWeight[5, {0, 1, -1, 0, 0}], 
 
>      finiteWeight[5, {0, 0, 1, -1, 0}], finiteWeight[5, {0, 0, 0, 1, -1}]}]\
 
>     , affineWeight[5, finiteWeight[5, {-1, 0, 0, 0, 1}], 0, 1], 
 
>    {affineWeight[5, finiteWeight[5, {1, -1, 0, 0, 0}], 0, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 1, -1, 0, 0}], 0, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 0, 1, -1, 0}], 0, 0], 
 
>     affineWeight[5, finiteWeight[5, {0, 0, 0, 1, -1}], 0, 0]}]

Out[95]//MatrixForm= 2    -1   0    0    -1

                     -1   2    -1   0    0

                     0    -1   2    -1   0

                     0    0    -1   2    -1

                     -1   0    0    -1   2

Out[94]//MatrixForm= 2    -1   -1

                     -1   2    -1

                     -1   -1   2

Out[93]= {{2, -1, -1}, {-1, 2, -1}, {-1, -1, 2}}

Out[92]= {{2, -1}, {-1, 2}}

Out[91]= cartanMatrix[makeFiniteRootSystem[A, 2]]

Out[88]= affineWeight[2, finiteWeight[2, {2, 0}], 2, 0]

                                          2       2  2
Out[87]= affineWeight[2, finiteWeight[2, {-, 0}], -, -]
                                          3       3  3

Out[85]= False

Out[84]= True

Out[83]= False

Out[82]= True


?Times

x * y * z, x z or x y z represents a product of terms. Mutliplication by
                                numbers is defined for the elements of weight
                                space of affine and finite-dimensional Lie
                                algebras. 
                                 
                              For example
                                makeFiniteWeight[{1,2,3}]*2==makeFiniteWeight[
                                {2,4,6}]]

?finiteWeight

v=3

v+=4

v="aa"

v<>="bb"

Syntax::sntxf: "v<>" cannot be followed by "="bb"".

Out[76]= aa

Out[75]= 7

Out[74]= 7

Out[73]= 3

finiteWeight[dimension_?NomberQ,coordinates_standardBase] represents 
    vector in weight space of finite-dimensional Lie algebra.
 
    finiteWeight[dimension] returns dimension of the space, where weight
   vector is embedded 
    (i.e. for sl_n it is n+1.
 
     finiteWeight[standardBase] returns standard base coordinates of weight of
   finite-dimensional Lie algebra

?==

orbit[b2][{makeFiniteWeight[{2,-1}]}]

ChebyshevU[5,x]

5!

M[p_,k_]/;Mod[p,2]==Mod[k,2]:=p!(k+1)/((p-k)/2)!/((p+k+2)/2)!;
M[p_,k_]/;Mod[p,2]!=Mod[k,2]:=0

Clear[M]

MM:=Table[M[i,j],{i,10},{j,10}]

MM//MatrixForm

Out[53]//MatrixForm= 1   0   0   0   0

                     0   1   0   0   0

                     2   0   1   0   0

                     0   3   0   1   0

                     5   0   4   0   1

Inverse[MM].((2*x)^Range[10])

                  2            3       2       4            3       5
Out[69]= {2 x, 4 x , -4 x + 8 x , -12 x  + 16 x , 6 x - 32 x  + 32 x , 
 
         2       4       6             3        5        7
>    24 x  - 80 x  + 64 x , -8 x + 80 x  - 192 x  + 128 x , 
 
          2        4        6        8
>    -40 x  + 240 x  - 448 x  + 256 x , 
 
                 3        5         7        9
>    10 x - 160 x  + 672 x  - 1024 x  + 512 x , 
 
         2        4         6         8         10
>    60 x  - 560 x  + 1792 x  - 2304 x  + 1024 x  }

                          2            3          2       4
Out[66]= {1, 2 x, -2 + 4 x , -6 x + 8 x , 3 - 16 x  + 16 x }

Range[0,5]

Out[65]= {0, 1, 2, 3, 4, 5}

                  2            3       2       4            3       5
Out[57]= {2 x, 4 x , -4 x + 8 x , -12 x  + 16 x , 6 x - 32 x  + 32 x }

ChebyshevU[#,x]&/@ Range[10]

                       2            3          2       4
Out[70]= {2 x, -1 + 4 x , -4 x + 8 x , 1 - 12 x  + 16 x , 
 
               3       5           2       4       6
>    6 x - 32 x  + 32 x , -1 + 24 x  - 80 x  + 64 x , 
 
                3        5        7          2        4        6        8
>    -8 x + 80 x  - 192 x  + 128 x , 1 - 40 x  + 240 x  - 448 x  + 256 x , 
 
                 3        5         7        9
>    10 x - 160 x  + 672 x  - 1024 x  + 512 x , 
 
              2        4         6         8         10
>    -1 + 60 x  - 560 x  + 1792 x  - 2304 x  + 1024 x  }

                          2            3          2       4
Out[67]= {1, 2 x, -1 + 4 x , -4 x + 8 x , 1 - 12 x  + 16 x }

                       2            3          2       4            3       5
Out[59]= {2 x, -1 + 4 x , -4 x + 8 x , 1 - 12 x  + 16 x , 6 x - 32 x  + 32 x }

Out[58]= {1, 2, 3, 4, 5}

x^Range[10]

              2   3   4   5   6   7   8   9   10
Out[56]= {x, x , x , x , x , x , x , x , x , x  }

Out[55]= {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}

Out[54]//MatrixForm= 1    0    0    0    0

                     0    1    0    0    0

                     -2   0    1    0    0

                     0    -3   0    1    0

                     3    0    -4   0    1

M[3,3]

Simplify[Plus@@Table[M[5,i]*ChebyshevU[i,x],{i,0,5}]]

             5
Out[51]= 32 x

                    3       5                3
Out[50]= 16 x - 32 x  + 32 x  + 4 (-4 x + 8 x )

                                   3                3       5
Out[49]= {0, 10 x, 0, 4 (-4 x + 8 x ), 0, 6 x - 32 x  + 32 x }

                                   3                3       5
Out[48]= {0, 10 x, 0, 4 (-4 x + 8 x ), 0, 6 x - 32 x  + 32 x }

                                   3                3       5
Out[47]= {0, 10 x, 0, 4 (-4 x + 8 x ), 0, 6 x - 32 x  + 32 x }

                                      2
           1024         1024 (-1 + 4 x )               3
Out[44]= {------, 10 x, ----------------, 4 (-4 x + 8 x ), 
          105 Pi             63 Pi
 
                   2       4
     5120 (1 - 12 x  + 16 x )            3       5
>    ------------------------, 6 x - 32 x  + 32 x }
              693 Pi

Out[43]= 1

Out[42]= 2

Out[41]//MatrixForm= 1    0    0    0    0

                     0    1    0    0    0

                     -2   0    1    0    0

                     0    -3   0    1    0

                     3    0    -4   0    1

Out[40]//MatrixForm= 1   0   0   0   0

                     0   1   0   0   0

                     2   0   1   0   0

                     0   3   0   1   0

                     5   0   4   0   1

Out[38]//MatrixForm= 0   0   0   0   0

                     0   1   0   0   0

                     0   0   0   0   0

                     0   3   0   1   0

                     0   0   0   0   0

                                   1
Out[37]= {{0, 0, 0, 0, 0}, {0, ----------, 0, 0, 0}, {0, 0, 0, 0, 0}, 
                               MatrixForm
 
             3              1
>    {0, ----------, 0, ----------, 0}, {0, 0, 0, 0, 0}}
         MatrixForm     MatrixForm

?Mod

?&&

e  && e  && ... is the logical AND function. It evaluates its arguments in
 1     2
    order, giving False immediately if any of them are False, and True if they
    are all True. 

Mod[m, n] gives the remainder on division of m
     by n. Mod[m, n, d] uses an offset d. 


Out[33]= {{15, 0, 0, 0, 0}, {-45, 15, 0, 0, 0}, {90, -75, 15, 0, 0}, 
 
>    {-150, 225, -105, 15, 0}, {225, -525, 420, -135, 15}}

Out[32]= 5 Null

Inverse[MM]//MatrixForm

Out[31]//MatrixForm= 1     0     0     0     0

                     -3    1     0     0     0

                     6     -5    1     0     0

                     -10   15    -7    1     0

                     15    -35   28    -9    1

Out[30]//MatrixForm= 
 
>   MatrixInverse[{{1, 0, 0, 0, 0}, {3, 1, 0, 0, 0}, {9, 5, 1, 0, 0}, 
     
    >    {28, 20, 7, 1, 0}, {90, 75, 35, 9, 1}}]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Out[29]//MatrixForm= 
 


>   1                 ComplexInfinity   ComplexInfinity   ComplexInfinity
     


    >       ComplexInfinity

    1
    -
    3                 1                 ComplexInfinity   ComplexInfinity
     


    >       ComplexInfinity

    1                 1
    -                 -
    9                 5                 1                 ComplexInfinity
     


    >       ComplexInfinity

    1                 1                 1
    --                --                -
    28                20                7                 1
     


    >       ComplexInfinity

    1                 1                 1                 1
    --                --                --                -
    90                75                35                9
     


    >       1

Out[28]//MatrixForm= 1    0    0    0    0

                     3    1    0    0    0

                     9    5    1    0    0

                     28   20   7    1    0

                     90   75   35   9    1

Out[27]= {{1, 0, 0, 0, 0}, {3, 1, 0, 0, 0}, {9, 5, 1, 0, 0}, 
 
>    {28, 20, 7, 1, 0}, {90, 75, 35, 9, 1}}

               8        -8         8        -8         8
Out[25]= {{1, ----, 0, -----, 0, -----, 0, -----, 0, ------}, 
              5 Pi     21 Pi     45 Pi     77 Pi     117 Pi
 
       64        128        -64         256        -320
>    {-----, 1, ------, 0, ------, 0, -------, 0, -------, 0}, 
      15 Pi     105 Pi     315 Pi     3465 Pi     9009 Pi
 
          192       64        -64         192         -64
>    {2, -----, 1, -----, 0, ------, 0, -------, 0, -------}, 
         35 Pi     63 Pi     495 Pi     5005 Pi     4095 Pi
 
       1024       2048       1024        -4096         1024
>    {------, 3, ------, 1, -------, 0, --------, 0, --------, 0}, 
      105 Pi     315 Pi     1155 Pi     45045 Pi     45045 Pi
 
         1024       5120       1024        -1024         1024
>    {5, -----, 4, ------, 1, -------, 0, --------, 0, --------}, 
         63 Pi     693 Pi     1287 Pi     15015 Pi     69615 Pi
 
       8192      16384       8192        32768         -8192
>    {------, 9, ------, 5, -------, 1, --------, 0, ---------, 0}, 
      315 Pi     693 Pi     1001 Pi     45045 Pi     153153 Pi
 
           8192        40960       57344        8192         -8192
>    {14, ------, 14, -------, 6, -------, 1, --------, 0, ---------}, 
          165 Pi      1287 Pi     6435 Pi     12155 Pi     188955 Pi
 
      262144       524288       262144       1048576       262144
>    {-------, 28, -------, 20, -------, 7, ---------, 1, ---------, 0}, 
      3465 Pi      6435 Pi      6435 Pi     109395 Pi     415701 Pi
 
          786432       262144       1835008       2359296       262144
>    {42, -------, 48, -------, 27, --------, 8, ---------, 1, ---------}, 
          5005 Pi      2145 Pi      36465 Pi     230945 Pi     440895 Pi
 
      2097152      4194304       2097152        8388608      10485760
>    {-------, 90, --------, 75, --------, 35, ---------, 9, ---------, 1}}
      9009 Pi      15015 Pi      12155 Pi      138567 Pi     969969 Pi

Out[23]= 2

Out[21]= 120

                   3       5
Out[20]= 6 x - 32 x  + 32 x

ChebyshevU[n, x] gives the Chebyshev polynomial of the second kind U (x). 
                                                                    n

ChebyshevDistance ChebyshevT        ChebyshevU

Out[17]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 2}], finiteWeight[2, {2, -1}]}, 
 
>    {finiteWeight[2, {-1, 2}], finiteWeight[2, {1, -2}]}, 
 
>    {finiteWeight[2, {-2, 1}], finiteWeight[2, {-1, -2}]}, 
 
>    {finiteWeight[2, {-2, -1}]}}

partialOrbit[b2][{makeFiniteWeight[{2,-1}]}]

Out[16]= {{finiteWeight[2, {2, -1}]}, {finiteWeight[2, {-1, 2}]}, 
 
>    {finiteWeight[2, {-1, -2}]}, {finiteWeight[2, {-2, -1}]}}

Out[15]= {{finiteWeight[2, {1, 0}]}, {finiteWeight[2, {0, 1}]}, 
 
>    {finiteWeight[2, {0, -1}]}, {finiteWeight[2, {-1, 0}]}}

Out[14]= {{finiteWeight[2, {-1, 0}]}}



Out[11]= partialOrbit[b2][{finiteWeight[2, {-1, 0}]}]

Out[10]= partialOrbit[b2][finiteWeight[2, {-1, 0}]]

lhs == rhs returns True if lhs and rhs
      are identical. 
      It is defined for weights of finite and affine Lie algebras

?Dot

?Plus

x + y + z represents a sum of terms. 
          It is defined for weights of finite and affine Lie algebras

x + y + z represents a sum of terms. 

a . b . c or Dot[a, b, c] gives products of vectors, matrices and tensors. Dot
     product for finite weights

Dot product for finite weightsDot product for finite weights

Dot product for finite weights

a . b . c or Dot[a, b, c] gives products of vectors, matrices and tensors. 

finiteWeight[dimension_?NomberQ,coordinates_standardBase] represents vector in
   weight space of finite-dimensional Lie algebra.
 finiteWeight[dimension] returns dimension of the space, where weight vector
   is embedded (i.e. for sl_n it is n+1.
 finiteWeight[standardBase] returns standard base coordinates of weight of
   finite-dimensional Lie algebra

finiteWeight[dimension_?NomberQ,coordinates_standardBase] represents vector in
   weight space of finite-dimensional Lie algebra. Dimension is saved
   separately for simplification of testing

finiteWeight[2, {-1, -1}]
finiteWeight[2, {-1, 0}]
finiteWeight[2, {0, -1}]
finiteWeight[2, {-1, 1}]
finiteWeight[2, {1, -1}]
finiteWeight[2, {0, 1}]
finiteWeight[2, {1, 1}]
finiteWeight[2, {1, 0}]
                 3  1
finiteWeight[2, {-, -}]
                 2  2
{Null}
                 1  3
finiteWeight[2, {-, -}]
                 2  2
                 3    1
finiteWeight[2, {-, -(-)}]
                 2    2
{Null, Null}
                   1   3
finiteWeight[2, {-(-), -}]
                   2   2
                 1    3
finiteWeight[2, {-, -(-)}]
                 2    2
{Null, Null}
                   3   1
finiteWeight[2, {-(-), -}]
                   2   2
                   1     3
finiteWeight[2, {-(-), -(-)}]
                   2     2
{Null, Null}
                   3     1
finiteWeight[2, {-(-), -(-)}]
                   2     2
{Null}

Map[toFundamentalChamber[makeSimpleRootSystem[B,2]],{makeFiniteWeight[{-1,-1}],makeFiniteWeight[{-2,-1}]}]

Out[4]= {finiteWeight[2, {1, 1}], finiteWeight[2, {2, 1}]}

Expect["We can use this and other functions for mapping",True,
       Map[toFundamentalChamber[makeSimpleRootSystem[B,2]],{makeFiniteWeight[{-1,-1}],makeFiniteWeight[{-2,-1}]}]==
       {makeFiniteWeight[{1, 1}], makeFiniteWeight[{2, 1}]}]

b2=makeSimpleRootSystem[B,2];

orbit[b2][{makeFiniteWeight[{1,0}],makeFiniteWeight[{1,1}]}]

Out[17]= {{finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}]}, 
 
>    {finiteWeight[2, {0, 1}], finiteWeight[2, {1, -1}]}, 
 
>    {finiteWeight[2, {-1, 1}], finiteWeight[2, {0, -1}]}, 
 
>    {finiteWeight[2, {-1, -1}], finiteWeight[2, {-1, 0}]}}

Out[16]= {{finiteWeight[2, {1, 0}]}, {finiteWeight[2, {0, 1}]}, 
 
>    {finiteWeight[2, {0, -1}]}, {finiteWeight[2, {-1, 0}]}}

Out[15]= {{finiteWeight[2, {1, 0}]}, {finiteWeight[2, {0, 1}]}, 
 
>    {finiteWeight[2, {0, -1}]}, {finiteWeight[2, {-1, 0}]}}

Out[12]= {{finiteWeight[2, {1, 0}]}, {finiteWeight[2, {0, 1}]}, 
 
>    {finiteWeight[2, {0, -1}]}, {finiteWeight[2, {-1, 0}]}}

Out[11]= {{finiteWeight[2, {0, 0}]}}

Out[9]= {{finiteWeight[2, {0, 0}]}, {}}

Out[8]= {{finiteWeight[2, {1, 0}]}, {finiteWeight[2, {0, 1}]}, 
 
>    {finiteWeight[2, {0, -1}]}, {finiteWeight[2, {-1, 0}]}, {}}

Out[7]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 2}], finiteWeight[2, {2, -1}]}, 
 
>    {finiteWeight[2, {-1, 2}], finiteWeight[2, {1, -2}]}, 
 
>    {finiteWeight[2, {-2, 1}], finiteWeight[2, {-1, -2}]}, 
 
>    {finiteWeight[2, {-2, -1}]}, {}}

Out[3]= {toFundamentalChamber[makeFiniteRootSystem[B, 2]][finiteWeight[2, 
 
>      {-1, -1}]], toFundamentalChamber[makeFiniteRootSystem[B, 2]][
 
>     finiteWeight[2, {-2, -1}]]}

finiteWeight[2, {1, -1}]
finiteWeight[2, {0, 1}]
finiteWeight[2, {1, 1}]
finiteWeight[2, {1, 0}]
                 3  1
finiteWeight[2, {-, -}]
                 2  2
{Null}
                 1  3
finiteWeight[2, {-, -}]
                 2  2
                 3    1
finiteWeight[2, {-, -(-)}]
                 2    2
{Null, Null}
                   1   3
finiteWeight[2, {-(-), -}]
                   2   2
                 1    3
finiteWeight[2, {-, -(-)}]
                 2    2
{Null, Null}
                   3   1
finiteWeight[2, {-(-), -}]
                   2   2
                   1     3
finiteWeight[2, {-(-), -(-)}]
                   2     2
{Null, Null}
                   3     1
finiteWeight[2, {-(-), -(-)}]
                   2     2
{Null}
{}

fundamentalWeights[rs_affineRootSystem]:=Map[makeAffineWeight[#[[1]],#[[2]],0]&,
					     Transpose[{Prepend[fundamentalWeights[rs[finiteRootSystem]],
								0*rs[finiteRootSystem][simpleRoot][1]],
							comarks[rs]}]]

rs=b2a


Transpose[{Prepend[fundamentalWeights[rs[finiteRootSystem]],
								0*rs[finiteRootSystem][simpleRoot][1]],
							comarks[rs]}]

                  
Out[41]= {{finiteWeight[2, {0, 0}], finiteWeight[2, {1, 0}], 
 
                       1  1
>     finiteWeight[2, {-, -}]}, {1, 1, 1}}
                       2  2

Prepend[fundamentalWeights[rs[finiteRootSystem]],
								0*rs[finiteRootSystem][simpleRoot][1]]

         
Out[40]= {finiteWeight[2, {0, 0}], finiteWeight[2, {1, 0}], 
 
                      1  1
>    finiteWeight[2, {-, -}]}
                      2  2


(* finiteWeight/:0*y_finiteWeight:=makeFiniteWeight[0*y[standardBase]]; *)

fundamentalWeights[b2a]

Out[49]= {affineWeight[2, finiteWeight[2, {0, 0}], 1, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 1, 0], 
 
                                      1  1
>    affineWeight[2, finiteWeight[2, {-, -}], 1, 0]}
                                      2  2

comarks[b2a]

Out[39]= {1, 1, 1}

marks[rs]*Map[#.#/2&,rs[simpleRoots]]

b2a[simpleRoots]

Set::write: Tag Times in {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
      affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
      affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]} Out[302] is Protected.

          1 . 1  1 . 1
Out[31]= {-----, -----, 1 . 1}
            2      2

Set::write: Tag Times in {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
      affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
      affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]} Out[302] is Protected.

          1 . 1  1 . 1
Out[30]= {-----, -----, 1 . 1}
            2      2

Out[29]= {1, 1, 2}

rho[rs_affineRootSystem]:=Plus@@fundamentalWeights[rs]



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

weightSystem[b2a][rho[b2a],1]

rh=rho[b2a];
rs=b2a;
highestWeight=2*rh;gradelimit=1;
toFundamentalChamber[rs]/@Sort[ Rest[Flatten[weightSystem[rs][highestWeight,gradelimit]]], #1.rh>#2.rh&]

[Calculating...]

                           
Out[78]= {affineWeight[2, finiteWeight[2, {3, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 6, 0]}


mults[highestWeight]=1;
insideQ:=IntegerQ[mults[toFC[#]]]&;

                           
Out[77]= {affineWeight[2, finiteWeight[2, {3, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 6, 0]}

                           
Out[76]= {affineWeight[2, finiteWeight[2, {3, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {3, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 2}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {2, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 1}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {1, 0}], 6, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 0}], 6, 0]}

                           
                                           3  1
Out[75]= {affineWeight[2, finiteWeight[2, {-, -}], 3, 0], 
                                           2  2
 
                                      1  1
>    affineWeight[2, finiteWeight[2, {-, -}], 3, 0]}
                                      2  2

                           
                                           1  1
Out[74]= {affineWeight[2, finiteWeight[2, {-, -}], 3, 0]}
                                           2  2

Out[66]= mults$123

mts:=racahMultiplicities[b2a][rho[b2a],1];

mts[makeAffineWeight[{1/2,1/2},3,0]]


Out[71]= {affineWeight[2, finiteWeight[2, {12, 4}], 0, 0], -1}

                                            3  1
Out[72]= {{affineWeight[2, finiteWeight[2, {-, -}], 3, 0]}, 
                                            2  2
 
                                       1  1
>    {affineWeight[2, finiteWeight[2, {-, -}], 3, 0]}, {}}
                                       2  2

                                           3  1
Out[69]= {affineWeight[2, finiteWeight[2, {-, -}], 3, 0], 
                                           2  2
 
                                      1  1
>    affineWeight[2, finiteWeight[2, {-, -}], 3, 0]}
                                      2  2

[Calculating...]

[Calculating...]

Out[62]= {{affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {0, 1}], 0, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, -1}], 0, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {1, 2}], 0, 0], 1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, 0}], 0, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {2, 2}], 0, 0], -1}, 
 
>    {affineWeight[2, finiteWeight[2, {3, 1}], 0, 0], 1}}



Module[{b2=makeSimpleRootSystem[B,2],fm,rm},
       fm=freudenthalMultiplicities[b2][makeFiniteWeight[{2,1}]];
       rm=freudenthalMultiplicities[b2][makeFiniteWeight[{2,1}]];
       Expect["Racah and Freudenthal formulae should give the same result",rm[#]&/@keys[rm],
	      fm[#]&/@keys[fm]]]

Module[{b2=makeSimpleRootSystem[B,2]},
       Print[positiveRoots[b2]];
       weightSystem[b2][makeFiniteWeight[{2,1}]]]

                  {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}], finiteWeight[2, {1, 1}], 
 
>   finiteWeight[2, {1, 0}]}

Out[11]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}}

                {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}], finiteWeight[2, {1, 1}], 
 
>   finiteWeight[2, {1, 0}]}

Out[9]= {{finiteWeight[2, {2, 1}]}, 
 
>    {finiteWeight[2, {1, 0}], finiteWeight[2, {1, 1}], 
 
>     finiteWeight[2, {2, 0}]}, {finiteWeight[2, {0, 0}]}, {}}

                positiveRoots[finiteRootSystem[2, 
 
>    {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}]}]]

Outer::heads: Heads List and Times at positions 3 and 2
     are expected to be the same.

Out[6]= {{finiteWeight[2, {2, 1}]}, {}}

        
Outer::heads: Heads List and Times at positions 3 and 2
     are expected to be the same.

Out[5]= {{finiteWeight[2, {2, 1}]}, {}}

                                
Outer::heads: Heads List and Times at positions 3 and 2
     are expected to be the same.

Throw::nocatch: 
   Uncaught Throw[Weights of [2,1] module of B2: GOT UNEXPECTED VALUE False
      INSTEAD OF True, assertion exception] returned to top level.

Out[2]= Hold[Throw[Weights of [2,1] module of B2: GOT UNEXPECTED VALUE False\
 
>      INSTEAD OF True, assertion exception]]

positiveRoots[rs]

weightSystem[rs_?rootSystemQ][higestWeight_?weightQ]:=Module[{minusPosRoots=-positiveRoots[rs]},
							     Most[NestWhileList[Function[x,Complement[
										 Cases[Flatten[Outer[Plus,minusPosRoots,x]],y_/;
										       And[checkGrade[rs][y],And@@(#.y>=0&/@rs[simpleRoots])]]
										 ,x]],{higestWeight},#=!={}&]]];



higestRoot::"usage"="returns highest root of root system";
higestRoot[rs_finiteRootSystem]:=toFundamentalChamber[rs][rs[simpleRoot][Ordering[(#.#&)/@rs[simpleRoots],-1][[1]]]]

?higestRoot

makeAffineExtension[makeSimpleRootSystem[B,2]][simpleRoots]

Out[7]= {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 0]}

Out[6]= {affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1], 
 
>    affineWeight[2, finiteWeight[2, {1, -1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {0, 1}], 0, 0], 
 
>    affineWeight[2, finiteWeight[2, {-1, -1}], 0, 1]}

returns highest root of root system


orthogonalSubsystem[rs_?rootSystemQ,subs_?rootSystemQ]:=Cases[Flatten[positiveRoots[rs]], z_?weightQ /; Or[z.#==0& /@ subs[simpleRoots]]]

subs=makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]

Out[28]= finiteRootSystem[1, {finiteWeight[2, {1, 1}]}]

Flatten[positiveRoots[b2]] /. z_?weightQ -> z.#==0& /@ subs[simpleRoots]

Select[Flatten[positiveRoots[b2]], Function[z,Or[(z.#==0&) /@ subs[simpleRoots]]]]

Select[Flatten[positiveRoots[b2]], #.makeFiniteWeight[{1,1}]==0&]

Out[41]= 2 finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>      finiteWeight[2, {0, 1}]}]

Out[40]= {finiteWeight[2, {1, -1}]}

z=finiteWeight[2, {1, 1}]

(z.#==0&) /@ subs[simpleRoots]

[Calculating...]

Out[33]= {}

Out[32]= {}

Out[31]= {}

Out[30]= {{True, False, False, False}}

Out[29]= {{finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}], 
 
>       finiteWeight[2, {1, 1}], finiteWeight[2, {1, 0}]} . 
 
>      finiteWeight[2, {1, 1}] == 0}

Out[27]= {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}], 
 
>    finiteWeight[2, {1, 1}], finiteWeight[2, {1, 0}]}

b2=makeSimpleRootSystem[B,2]

Out[10]= finiteRootSystem[2, {finiteWeight[2, {1, -1}], 
 
>     finiteWeight[2, {0, 1}]}]

orthogonalSubsystem[b2,makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}]]

Out[26]= {}

Out[22]= {}

Out[20]= {{True}, {False}, {False}, {False}}

Out[18]= {}

positiveRoots[b2]

Out[13]= {finiteWeight[2, {1, -1}], finiteWeight[2, {0, 1}], 
 
>    finiteWeight[2, {1, 1}], finiteWeight[2, {1, 0}]}

Out[12]= {}

Out[11]= finiteRootSystem[1, {finiteWeight[2, {1, 1}]}]


makeFiniteRootSystem[{makeFiniteWeight[{1,1}]}][simpleRoots]

Out[14]= {finiteWeight[2, {1, 1}]}