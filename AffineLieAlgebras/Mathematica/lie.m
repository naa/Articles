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
fundamentalWeights[{simpleRoots__standardBase}]:=Plus@@(Inverse[cartanMatrix[{simpleRoots}]]*{simpleRoots})

rho[{simpleRoots__standardBase}]:=Plus@@fundamentalWeights[simpleRoots]

rho[b2]

                      3  1
Out[98]= standardBase[-, -]
                      2  2


orbit[{simpleRoots__standardBase}][weights__standardBase]:=Module[{total={}},
       NestWhileList[Function[x,Complement[Union[Flatten[Outer[#1[#2]&,reflection/@{simpleRoots},x]]],total]],
		     {weights},
		     Function[y,Module[{t=y=!={}},total=Union[total,y];t]]]]

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