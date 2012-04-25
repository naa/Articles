ExpandNCM[(h : NonCommutativeMultiply)[a___, b_Plus, c___]] :=  Distribute[h[a, b, c], Plus, h, Plus, ExpandNCM[h[##]] &];
ExpandNCM[a_] := ExpandAll[a];

L/:NonCommutativeMultiply[L[0],phi] := h*phi;
L/:NonCommutativeMultiply[L[n_Integer],phi]/;n>0 := 0;
L/:NonCommutativeMultiply[L[n_Integer],phi]/;n>0 := 0;
J/:NonCommutativeMultiply[J[n_Integer][a_Integer],phi]/;n>0 := 0;
J/:NonCommutativeMultiply[J[0][3],phi] := j*phi;
L/:NonCommutativeMultiply[L[n_Integer],L[m_Integer]]/;And[n>=0, m<=0, n==-m] := c/12*(n^3-n);
L/:NonCommutativeMultiply[L[n_Integer],L[m_Integer]]/;And[n>=0, m<=0] := (n-m)*L[n+m];
L/:NonCommutativeMultiply[L[n_Integer],J[m_Integer][a_Integer]]:=-m*J[n+m][a];
J/:NonCommutativeMultiply[J[m_Integer][a_Integer],J[n_Integer][b_Integer]]/;And[m>=0,n==-m,a==b] := J[n][a]**J[m][a]+k*m;
J/:NonCommutativeMultiply[J[m_Integer][a_Integer],J[n_Integer][b_Integer]]/;And[m>=0,n<=0] := Sum[I*Signature[a,b,c]*J[m+n][c],{c,3}];
L/:NonCommutativeMultiply[q_?NumericQ,L[n_Integer]]:=q*L[n];
J/:NonCommutativeMultiply[q_?NumericQ,J[n_Integer][a_Integer]]:=q*J[n][a];
phi/:NonCommutativeMultiply[q_?NumericQ,phi]:=q*psi;
psi=(-2*L[-2]+1/2*kappa*L[-1]**L[-1]+1/2*tau*Sum[J[-1][a]**J[-1][a],{a,3}])**phi;

ExpandNCM[L[2]**psi]

                                              kappa L[-1] ** L[-1]
Out[11]= L[2] ** (-2 L[-2]) ** phi + L[2] ** (--------------------) ** phi + 
                                                       2
 
              tau J[-1][1] ** J[-1][1]   tau J[-1][2] ** J[-1][2]
>    L[2] ** (------------------------ + ------------------------ + 
                         2                          2
 
        tau J[-1][3] ** J[-1][3]
>       ------------------------) ** phi
                   2

                             kappa L[-1] ** L[-1]
Out[10]= L[2] ** (-2 L[-2] + -------------------- + 
                                      2
 
       tau (J[-1][1] ** J[-1][1] + J[-1][2] ** J[-1][2] + 
        
       >     J[-1][3] ** J[-1][3])
>      ------------------------------------------------------------------------
                                          2
 
>      ) ** phi

                                                                       
SetDelayed::write: Tag NonCommutativeMultiply in L[0] ** phi is Protected.

SetDelayed::write: 
   Tag NonCommutativeMultiply in L[n_Integer] ** phi /; n > 0 is Protected.

SetDelayed::write: 
   Tag NonCommutativeMultiply in L[n_Integer] ** phi /; n > 0 is Protected.

General::stop: Further output of SetDelayed::write
     will be suppressed during this calculation.

J[n][a]**J[m][a]+k*m;


?I

I represents the imaginary unit Sqrt[-1]. 

?Sum

                                     i
Sum[f, {i, i   }] evaluates the sum   max  f
            max                     \[Sum]
                                    i = 1
    . Sum[f, {i, i   , i   }] starts with i = i
                  min   max                    min
      . Sum[f, {i, i   , i   , di}]
                    min   max
        uses steps di. Sum[f, {i, {i , i , ...}}]
                                    1   2
          uses successive values i , i , ...
                                  1   2
          .Sum[f, {i, i   , i   }, {j, j   , j   }, ...]
                       min   max        min   max
                                         i        j
                                          max      max
            evaluates the multiple sum  \[Sum]   \[Sum]  ...f
                                       i = i    j = j
                                            min      min
            . Sum[f, i] gives the indefinite sum \[Sum] f.
                                                   i


kappa:=((8*h+c)*(cm-cn)-6*h*(k*dimg-k*xe*dima))/(3*h*(cm-cn)-h*(2h+1)*(k*dimg-k*xe*dima))

tau:=12*h/(cm-cn)-2*h*(2*h+1)/(cm-cn)*kappa

Simplify[{kappa,tau} /. {h->cm/(2*(k+hg))-cn/(2*(k*xe+ha)), c->k*dimg/(2*(k+hg))-k*xe*dima/(2*(k*xe+ha))}]/.{cn->1,k->2,dimg->3,dima->1,hg->2,ha->0,xe->1}/.{cm->4}



Solve[{(3*kp-8)*h-c+tu*(k*dimg-xe*k*dima)==0,-12*h+2*kp*h*(2*h+1)+tu*(cm-cn)==0},{kp,tu}]/.{h->cm/(2*(k+hg))-cn/(2*(k*xe+ha)), c->k*dimg/(2*(k+hg))-k*xe*dima/(2*(k*xe+ha))}/.{cn->1,k->2,dimg->3,dima->1,hg->2,ha->0,xe->1}/.{cm->1}

                          9
Out[34]= {{kp -> 8, tu -> --}}
                          16

Syntax::sntxf: "(k*xe+ha))}/.{cn->1,k->2,dimg->3,dima->1,hg->2,ha->0,xe->1}/.{
      cm->" cannot be followed by "}".

                 92        29
Out[31]= {{kp -> --, tu -> --}}
                 11        88

                            3
Out[30]= {{kp -> 7, tu -> -(-)}}
                            4

                                                                          2
                   -3 (-c - 8 h) - 48 h         2 (c - 10 h + 2 c h + 16 h )
Out[29]= {{kp -> -(--------------------), tu -> ----------------------------}}
                   -9 h + 8 h (1 + 2 h)                  -1 + 16 h

                   -((cm - cn) (-c - 8 h)) - 12 h (dimg k - dima k xe)
Out[28]= {{kp -> -(---------------------------------------------------), 
                   -3 (cm - cn) h + 2 h (1 + 2 h) (dimg k - dima k xe)
 
>     tu -> 
 
                                                     2
                           2 (c - 10 h + 2 c h + 16 h )
>       ------------------------------------------------------------------}}
        -3 cm + 3 cn + 2 dimg k + 4 dimg h k - 2 dima k xe - 4 dima h k xe

??Solve

Solve[expr, vars] attempts to solve the system expr
     of equations or inequalities for the variables vars
     . Solve[expr, vars, dom] solves over the domain dom

       . Common choices of dom are Reals, Integers, and Complexes.

Attributes[Solve] = {Protected}
 
Options[Solve] = {Cubics -> True, GeneratedParameters -> C, 
    InverseFunctions -> Automatic, MaxExtraConditions -> 0, 
    Method -> Automatic, Modulus -> 0, Quartics -> True, 
    VerifySolutions -> Automatic, WorkingPrecision -> Infinity}

Solve[expr, vars] attempts to solve the system expr
     of equations or inequalities for the variables vars
     . Solve[expr, vars, dom] solves over the domain dom

       . Common choices of dom are Reals, Integers, and Complexes.

Set::write: Tag Plus in -c + h (-8 + 3 kp) + tu (dimg k - dima k xe)
     is Protected.

Set::write: Tag Plus in -12 h + 2 h (1 + 2 h) kp + (cm - cn) tu is Protected.

Solve::naqs: 0 && 0 is not a quantified system of equations and inequalities.

Out[25]= Solve[{0, 0}, {kp, tu}]

Solve::naqs: -c + h (-8 + 3 kp) + tu (dimg k - dima k xe) && 
     -12 h + 2 h (1 + 2 h) kp + (cm - cn) tu is not a quantified system of
     equations and inequalities.

Out[24]= Solve[{-c + h (-8 + 3 kp) + tu (dimg k - dima k xe), 
 
>     -12 h + 2 h (1 + 2 h) kp + (cm - cn) tu}, {kp, tu}]

             3
Out[23]= {1, -}
             4

          38  17
Out[22]= {--, --}
          21  28

               3
Out[21]= {8, -(-)}
               8

              29
Out[20]= {11, --}
              32

5/2-1/4

         9
Out[19]= -
         4

          5  1
Out[18]= {-, -}
          2  4

                31
Out[17]= {10, -(--)}
                16

                31
Out[15]= {10, -(--)}
                16

                           2                       2
          -4 (44 cm - 16 cm )  -32 + 152 cm - 32 cm
Out[14]= {-------------------, ---------------------}
           cm (-32 + 16 cm)       8 (-32 + 16 cm)

                          2
         -4 (44 cm - 16 cm )
Out[12]= -------------------
          cm (-32 + 16 cm)

Out[11]= 10

Out[10]= 7

Out[8]= -15

                        2
        4 (44 cm + 16 cm )
Out[7]= ------------------
         cm (-32 + 16 cm)

Out[6]= -34

                                   2
Out[5]= ((hg + k) (ha + k xe) (8 cm  (ha + k xe) + 
 
>        cn (8 cn (hg + k) + k (dimg (ha - 6 hg + k (-6 + xe)) + 
 
>              5 dima (hg + k) xe)) + 
 
>        cm (-8 cn (ha + hg + k + k xe) + 
 
>           k (dima xe (-6 ha + hg + k - 6 k xe) + 5 dimg (ha + k xe))))) / 
 
>    ((-(cn (hg + k)) + cm (ha + k xe)) 
 
>      (cm (ha + k xe) (3 hg + k (3 - dimg + dima xe)) - 
 
>        (hg + k) (k (dimg - dima xe) (ha + k xe) + 
 
>           cn (3 ha + k (-dimg + (3 + dima) xe)))))

                                     cm            cn
Out[4]= (6 (dimg k - dima k xe) (---------- - -------------) + 
                                 2 (hg + k)   2 (ha + k xe)
 
                  -(dimg k)      dima k xe            cm            cn
>      (cm - cn) (---------- + ------------- + 8 (---------- - -------------))
                  2 (hg + k)   2 (ha + k xe)      2 (hg + k)   2 (ha + k xe)
 
                             cm            cn
>      ) / (3 (cm - cn) (---------- - -------------) - 
                         2 (hg + k)   2 (ha + k xe)
 
                                 cm            cn
>      (dimg k - dima k xe) (---------- - -------------) 
                             2 (hg + k)   2 (ha + k xe)
 
                    cm            cn
>       (1 + 2 (---------- - -------------)))
                2 (hg + k)   2 (ha + k xe)



Simplify[-12*h+2*kp*h*(2*h+1)/.h->(6-kp)/(2*kp)]

Out[33]= 0

              6 - kp             6 (6 - kp)
Out[32]= (1 + ------) (6 - kp) - ----------
                kp                   kp



h:=k*(k+2)/(4*(N+2))-l^2/(4*N);
c:=(2*N-2)/(N+2);
cm:=k*(k+2)/2; 
cn:=l^2/2;

eqns:={(3*kp-8)*h-c+t*2*N==0, -12*h+2*kp*h*(2*h+1)+t*(cm-cn)==0};

soln=Solve[eqns,{kp,t}]

soln/.{N->3,k->1,l->-1}

                 260       28
Out[25]= {{kp -> ---, t -> ---}}
                 53        477


Simplify[soln/.{k->0,l->1}]

kp1:=4*(N+1)/(N+2);
t1:=4*N^3/(N+2)^3;
kp2:=4*(N+2)/(N+1);
t2:=2^(-4/N)*N^3/((N+2)*(N+1));

Table[{{i,j},soln,kp1,t1,kp2,t2}/.{N->4,k->i,l->j},{i,0,3},{j,0,3}]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                                1    10  32  24  16
Out[75]= {{{{0, 0}, {{kp -> Indeterminate, t -> -}}, --, --, --, --}, 
                                                8    3   27  5   15
 
                       184       25     10  32  24  16
>     {{0, 1}, {{kp -> ---, t -> ---}}, --, --, --, --}, 
                       31        124    3   27  5   15
 
                       52       4    10  32  24  16
>     {{0, 2}, {{kp -> --, t -> -}}, --, --, --, --}, 
                       7        7    3   27  5   15
 
                       248       169    10  32  24  16
>     {{0, 3}, {{kp -> ---, t -> ---}}, --, --, --, --}}, 
                       23        92     3   27  5   15
 
                       144       1     10  32  24  16
>    {{{1, 0}, {{kp -> ---, t -> --}}, --, --, --, --}, 
                       31        31    3   27  5   15
 
                       24       3     10  32  24  16
>     {{1, 1}, {{kp -> --, t -> --}}, --, --, --, --}, 
                       5        40    3   27  5   15
 
                       64       1    10  32  24  16
>     {{1, 2}, {{kp -> --, t -> -}}, --, --, --, --}, 
                       9        3    3   27  5   15
 
                       72       11    10  32  24  16
>     {{1, 3}, {{kp -> --, t -> --}}, --, --, --, --}}, 
                       7        8     3   27  5   15
 
                       39       1     10  32  24  16
>    {{{2, 0}, {{kp -> --, t -> --}}, --, --, --, --}, 
                       11       66    3   27  5   15
 
                       4296        1      10  32  24  16
>     {{2, 1}, {{kp -> ----, t -> ----}}, --, --, --, --}, 
                       1105       1020    3   27  5   15
 
                       84       4     10  32  24  16
>     {{2, 2}, {{kp -> --, t -> --}}, --, --, --, --}, 
                       19       57    3   27  5   15
 
                       6456       529    10  32  24  16
>     {{2, 3}, {{kp -> ----, t -> ---}}, --, --, --, --}}, 
                       671        732    3   27  5   15
 
                       16       1    10  32  24  16
>    {{{3, 0}, {{kp -> --, t -> -}}, --, --, --, --}, 
                       9        3    3   27  5   15
 
                       248       25     10  32  24  16
>     {{3, 1}, {{kp -> ---, t -> ---}}, --, --, --, --}, 
                       117       104    3   27  5   15
 
                       224       1     10  32  24  16
>     {{3, 2}, {{kp -> ---, t -> --}}, --, --, --, --}, 
                       69        23    3   27  5   15
 
                       8       1    10  32  24  16
>     {{3, 3}, {{kp -> -, t -> -}}, --, --, --, --}}}
                       3       8    3   27  5   15

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                        1    10  32  24  16
Out[74]= {{{{{kp -> Indeterminate, t -> -}}, --, --, --, --}, 
                                        8    3   27  5   15
 
               184       25     10  32  24  16
>     {{{kp -> ---, t -> ---}}, --, --, --, --}, 
               31        124    3   27  5   15
 
               52       4    10  32  24  16
>     {{{kp -> --, t -> -}}, --, --, --, --}, 
               7        7    3   27  5   15
 
               248       169    10  32  24  16
>     {{{kp -> ---, t -> ---}}, --, --, --, --}}, 
               23        92     3   27  5   15
 
               144       1     10  32  24  16
>    {{{{kp -> ---, t -> --}}, --, --, --, --}, 
               31        31    3   27  5   15
 
               24       3     10  32  24  16
>     {{{kp -> --, t -> --}}, --, --, --, --}, 
               5        40    3   27  5   15
 
               64       1    10  32  24  16
>     {{{kp -> --, t -> -}}, --, --, --, --}, 
               9        3    3   27  5   15
 
               72       11    10  32  24  16
>     {{{kp -> --, t -> --}}, --, --, --, --}}, 
               7        8     3   27  5   15
 
               39       1     10  32  24  16
>    {{{{kp -> --, t -> --}}, --, --, --, --}, 
               11       66    3   27  5   15
 
               4296        1      10  32  24  16
>     {{{kp -> ----, t -> ----}}, --, --, --, --}, 
               1105       1020    3   27  5   15
 
               84       4     10  32  24  16
>     {{{kp -> --, t -> --}}, --, --, --, --}, 
               19       57    3   27  5   15
 
               6456       529    10  32  24  16
>     {{{kp -> ----, t -> ---}}, --, --, --, --}}, 
               671        732    3   27  5   15
 
               16       1    10  32  24  16
>    {{{{kp -> --, t -> -}}, --, --, --, --}, 
               9        3    3   27  5   15
 
               248       25     10  32  24  16
>     {{{kp -> ---, t -> ---}}, --, --, --, --}, 
               117       104    3   27  5   15
 
               224       1     10  32  24  16
>     {{{kp -> ---, t -> --}}, --, --, --, --}, 
               69        23    3   27  5   15
 
               8       1    10  32  24  16
>     {{{kp -> -, t -> -}}, --, --, --, --}}}
               3       8    3   27  5   15

                  656       27    10  32  24  16
Out[73]= {{{kp -> ---, t -> --}}, --, --, --, --}
                  35        5     3   27  5   15

                  72       11    10  32  24  16
Out[72]= {{{kp -> --, t -> --}}, --, --, --, --}
                  7        8     3   27  5   15

                  64       1    10  32  24  16
Out[71]= {{{kp -> --, t -> -}}, --, --, --, --}
                  9        3    3   27  5   15

                  24       3     10  32  24  16
Out[70]= {{{kp -> --, t -> --}}, --, --, --, --}
                  5        40    3   27  5   15

                  144       1     10  32  24  16
Out[69]= {{{kp -> ---, t -> --}}, --, --, --, --}
                  31        31    3   27  5   15

                           3    10  32  24  16
Out[68]= {{{kp -> -1, t -> -}}, --, --, --, --}
                           2    3   27  5   15

                           25    10  32  24  16
Out[67]= {{{kp -> 19, t -> --}}, --, --, --, --}
                           4     3   27  5   15

                  248       169    10  32  24  16
Out[66]= {{{kp -> ---, t -> ---}}, --, --, --, --}
                  23        92     3   27  5   15

                  72       11    10  32  24  16
Out[65]= {{{kp -> --, t -> --}}, --, --, --, --}
                  7        8     3   27  5   15

                  248       25     10  32  24  16
Out[64]= {{{kp -> ---, t -> ---}}, --, --, --, --}
                  117       104    3   27  5   15

                  84       4     10  32  24  16
Out[63]= {{{kp -> --, t -> --}}, --, --, --, --}
                  19       57    3   27  5   15

                  24       3     10  32  24  16
Out[62]= {{{kp -> --, t -> --}}, --, --, --, --}
                  5        40    3   27  5   15

                  52       4    10  32  24  16
Out[61]= {{{kp -> --, t -> -}}, --, --, --, --}
                  7        7    3   27  5   15

                  184       25     10  32  24  16
Out[60]= {{{kp -> ---, t -> ---}}, --, --, --, --}
                  31        124    3   27  5   15

                  144       1     10  32  24  16
Out[59]= {{{kp -> ---, t -> --}}, --, --, --, --}
                  31        31    3   27  5   15

                  39       1     10  32  24  16
Out[58]= {{{kp -> --, t -> --}}, --, --, --, --}
                  11       66    3   27  5   15

                  84       4     10  32  24  16
Out[57]= {{{kp -> --, t -> --}}, --, --, --, --}
                  19       57    3   27  5   15

                  64       1    10  32  24  16
Out[56]= {{{kp -> --, t -> -}}, --, --, --, --}
                  9        3    3   27  5   15

                  24       3     10  32  24  16
Out[55]= {{{kp -> --, t -> --}}, --, --, --, --}
                  5        40    3   27  5   15

                  144       1     10  32  24  16
Out[54]= {{{kp -> ---, t -> --}}, --, --, --, --}
                  31        31    3   27  5   15


                                  2                      2      3
                 8 (2 + 14 N + 5 N )       8 + 28 N - 2 N  + 8 N
Out[52]= {{kp -> -------------------, t -> ----------------------}}
                 (2 + N) (-1 + 8 N)         2
                                           N  (2 + N) (-1 + 8 N)

                                  2                          2    3
                 2 (8 + 17 N + 5 N )       2 (16 + 20 N + 2 N  + N )
Out[51]= {{kp -> -------------------, t -> -------------------------}}
                 (2 + N) (-1 + 2 N)           2
                                             N  (2 + N) (-1 + 2 N)

                                                          2      3
                 4 (2 + N) (-2 + 5 N)       2 (2 + N - 5 N  + 2 N )
Out[50]= {{kp -> --------------------, t -> -----------------------}}
                                  2           2                 2
                   -10 + 9 N + 4 N           N  (-10 + 9 N + 4 N )

                 8 (2 + N) (-2 + 5 N)           -26 + 8 N
Out[49]= {{kp -> --------------------, t -> -----------------}}
                                  2                         2
                  -18 + 19 N + 8 N          -18 + 19 N + 8 N

                                                      2
                 (2 + N) (-7 + 5 N)        6 - 5 N + N
Out[48]= {{kp -> ------------------, t -> ---------------}}
                               2                   2    3
                   -6 + 3 N + N           -12 + 5 N  + N

                   -48 N       -16    -2 + 2 N
                   ----- - 4 (----- - --------)
                   2 + N      2 + N    2 + N
Out[47]= {{kp -> -(----------------------------), 
                                        4
                             8 N (1 + -----)
                      -24             2 + N
                     ----- + ---------------
                     2 + N        2 + N
 
                     2       3      4
              2 (24 N  - 20 N  + 4 N )
>     t -> ------------------------------}}
            2                          2
           N  (2 + N) (-48 + 24 N + 8 N )

                 10
Out[46]= {{kp -> --, t -> 0}}
                 3

                 520          2
Out[45]= {{kp -> ---, t -> -(---)}}
                 111         111

                 260       28
Out[44]= {{kp -> ---, t -> ---}}
                 53        477

                 208       242
Out[43]= {{kp -> ---, t -> ---}}
                 25        225

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                     2
Out[42]= {{kp -> Indeterminate, t -> --}}
                                     15

                 10
Out[41]= {{kp -> --, t -> 0}}
                 3

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                       8
Out[40]= {{kp -> ComplexInfinity, t -> --}}
                                       15

                 872       14
Out[39]= {{kp -> ---, t -> ---}}
                 275       495

                 8        14
Out[38]= {{kp -> --, t -> --}}
                 15       15

                 38       208
Out[37]= {{kp -> --, t -> ---}}
                 35       315

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                       8
Out[36]= {{kp -> ComplexInfinity, t -> --}}
                                       15

                 344       322
Out[35]= {{kp -> ---, t -> ---}}
                 25        75

                 208       242
Out[34]= {{kp -> ---, t -> ---}}
                 25        225

                 10
Out[33]= {{kp -> --, t -> 0}}
                 3

                 80       14
Out[32]= {{kp -> --, t -> ---}}
                 19       171

                 6680         46
Out[31]= {{kp -> ----, t -> -(---)}}
                 1729         819

                 1160       74
Out[30]= {{kp -> ----, t -> ---}}
                 143        117

                 520          2
Out[29]= {{kp -> ---, t -> -(---)}}
                 111         111

                 712       58
Out[28]= {{kp -> ---, t -> ---}}
                 115       207

                 260       28
Out[27]= {{kp -> ---, t -> ---}}
                 53        477

h/.{N->3,k->1,l->1}

         1
Out[26]= --
         15


         1
Out[24]= --
         15

                 260       28
Out[23]= {{kp -> ---, t -> ---}}
                 53        477

                 460       116
Out[22]= {{kp -> ---, t -> ---}}
                 33        33

                 208       242
Out[21]= {{kp -> ---, t -> ---}}
                 25        225

                 208       242
Out[20]= {{kp -> ---, t -> ---}}
                 25        225

c/.{N->3,k->1,l->1}

         4
Out[18]= -
         5

         1
Out[17]= -
         2



                 16
Out[13]= {{kp -> --, t -> 0}}
                 3

Out[12]= {{kp -> 3, t -> 0}}
