c:=15/2-3*(k+1/k);

Solve[c==7/10,{k}]

               3         5
Out[8]= {{k -> -}, {k -> -}}
               5         3

h:=(2-k)/(2*k);

h/.Out[8]

(6*k-3)/16/.Out[8]

          3   7
Out[12]= {--, --}
          80  16

          7  1
Out[11]= {-, --}
          6  10

Syntax::sntxf: "h." cannot be followed by "/Out[8]".

        15      1
Out[6]= -- - 3 (- + k)
        2       k

Out[4]= True

        15      1
Out[2]= -- - 3 (- + k)
        2       k

Sin[Pi]

Out[13]= 0

Solve[Sin[u]+Sqrt[2]*Sin[3/4*Pi-u]==0,u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                         -2                        2
Out[16]= {{u -> ArcCos[-------]}, {u -> -ArcCos[-------]}}
                       Sqrt[5]                  Sqrt[5]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                -Pi         Pi
Out[15]= {{u -> ---}, {u -> --}}
                 2          2

                                  Pi
Solve::naqs: Sin[u] - Sqrt[2] Sin[-- + u]
                                  4
     is not a quantified system of equations and inequalities.

                                    Pi
Out[14]= Solve[Sin[u] - Sqrt[2] Sin[-- + u], u]
                                    4

Solve[Sin[l+u]/Sin[l]==Sin[2*l-u]/Sin[2*l],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

Out[17]= {{u -> 0}}

Solve[Sin[l+u]/Sin[l]==Sin[3*l-u]/Sin[3*l],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

Out[18]= {{u -> 0}}

Solve[(1+v)/Sqrt[2]==Sin[u+Pi/4],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                               -Sqrt[2] - Sqrt[2] v
                -Pi - 4 ArcSin[--------------------]
                                        2
Out[23]= {{u -> ------------------------------------}, 
                                 4
 
                          Sqrt[2] + Sqrt[2] v
           -Pi + 4 ArcSin[-------------------]
                                   2
>    {u -> -----------------------------------}}
                            4

                                  Pi
                  Sqrt[2] - 2 Sin[-- + u]
                                  4
Out[19]= {{v -> -(-----------------------)}}
                          Sqrt[2]


Solve[(1-v)/Sqrt[2]==Sin[Pi/4-u],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                              Sqrt[2] - Sqrt[2] v
                Pi - 4 ArcSin[-------------------]
                                       2
Out[24]= {{u -> ----------------------------------}, 
                                4
 
                         -Sqrt[2] + Sqrt[2] v
           Pi + 4 ArcSin[--------------------]
                                  2
>    {u -> -----------------------------------}}
                            4

                                Pi
                Sqrt[2] - 2 Sin[-- - u]
                                4
Out[20]= {{v -> -----------------------}}
                        Sqrt[2]

Solve[Sqrt[2]==Sin[u+Pi/4]+Sin[Pi/4-u],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

Out[25]= {{u -> 0}}

Reduce[Reduce[Sqrt[2]==Sin[u+Pi/4]+Sin[Pi/4-u]],u]

Out[29]= C[1] \[Element] Integers && u == 2 Pi C[1]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

Out[28]= {{u -> 0}}

             Pi                       Pi
Out[27]= Sin[-- - u] == Sqrt[2] - Sin[-- + u]
             4                        4

Reduce::naqs: {u -> 0} is not a quantified system of equations and
     inequalities.

Out[26]= Reduce[{{u -> 0}}]

Solve[2==Sqrt[2]*Sin[Pi/4+u]+Sin[u],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                Pi                4
Out[30]= {{u -> --}, {u -> ArcCos[-]}}
                2                 5


Reduce[Sqrt[2]*Sin[Pi/4+u]==Sin[Pi/2-u],u]

Out[32]= C[1] \[Element] Integers && (u == 2 Pi C[1] || u == Pi + 2 Pi C[1])

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

Out[31]= {{u -> 0}}

Solve[Sin[u]==Sqrt[2]*Sin[Pi/4-u],u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                          -2                       2
Out[34]= {{u -> -ArcCos[-------]}, {u -> ArcCos[-------]}}
                        Sqrt[5]                 Sqrt[5]

Out[33]= C[1] \[Element] Integers && 
 
>    (u == -2 ArcTan[2 + Sqrt[5]] + 2 Pi C[1] || 
 
>      u == -2 ArcTan[2 - Sqrt[5]] + 2 Pi C[1])

N[Sin[ArcCos[2/Sqrt[5]]]]

Out[38]= 0.447214

            1
Out[35]= -------
         Sqrt[5]

N[Sqrt[2]*Sin[Pi/4+ArcCos[2/Sqrt[5]]]]

Out[42]= 1.34164

Out[39]= 1.34164

Out[37]= 0.447214

                     Pi             2
Out[36]= Sqrt[2] Sin[-- - ArcCos[-------]]
                     4           Sqrt[5]



N[Sin[Pi/2-ArcCos[2/Sqrt[5]]]]

Out[40]= 0.894427


w1=Sin[Pi/4+u];
w2=Sin[Pi/2-u];
w3=Sin[3*Pi/4-u]/Sin[3*Pi/4];
w4=Sin[Pi/2+u];

w5=Sin[u];
w6=Sqrt[2]*Sin[Pi/4-u];
w7=Sin[u];
w8=Sqrt[2]*Sin[Pi/4-u];

Solve[w1==w2,u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                        -Sqrt[2 + Sqrt[2]]
Out[46]= {{u -> -ArcCos[------------------]}, 
                                2
 
                  Sqrt[2 + Sqrt[2]]
>    {u -> ArcCos[-----------------]}}
                          2

Reduce[w1==w3,u]

Out[57]= C[1] \[Element] Integers && 
 
           -Pi                     3 Pi
>    (u == --- + 2 Pi C[1] || u == ---- + 2 Pi C[1])
            4                       4

Solve[w1==w4,u]

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

{w1,w2,w3,w4,w5,w6,w7,w8}

             Pi                           Pi
Out[7]= {Sin[-- + u], Cos[u], Sqrt[2] Sin[-- + u], Cos[u], Sin[u], 
             4                            4
 
                 Pi                           Pi
>    Sqrt[2] Sin[-- - u], Sin[u], Sqrt[2] Sin[-- - u]}
                 4                            4

Out[6]= {Sin[0.785398 + u], Cos[u], 1.41421 Sin[0.785398 + u], Cos[u], 
 
>    Sin[u], 1.41421 Sin[0.785398 - 1. u], Sin[u], 
 
>    1.41421 Sin[0.785398 - 1. u]}

Out[5]= {0.92388, 0.92388, 1.30656, 0.92388, 0.382683, 0.541196, 0.382683, 
 
>    0.541196}

Out[62]= {0.92388, 1.30656, 1.30656, 0.92388, 0.382683, 0.541196, 0.382683, 
 
>    0.541196}

{w1,w2,w3,w4}

              Pi                   Pi                   Pi
Out[67]= {Sin[-- + u], Sqrt[2] Sin[-- + u], Sqrt[2] Sin[-- + u], Cos[u]}
              4                    4                    4

{w5,w6,w7,w8}

                              Pi                           Pi
Out[68]= {Sin[u], Sqrt[2] Sin[-- - u], Sin[u], Sqrt[2] Sin[-- - u]}
                              4                            4

Out[66]= Cos[u]

                     Pi
Out[64]= Sqrt[2] Sin[-- + u]
                     4

                     Pi
Out[63]= Sqrt[2] Sin[-- + u]
                     4


              Pi               Pi               Pi       Pi       Pi
Out[61]= {Cos[--], Sqrt[2] Cos[--], Sqrt[2] Cos[--], Cos[--], Sin[--], 
              8                8                8        8        8
 
                 Pi       Pi               Pi
>    Sqrt[2] Sin[--], Sin[--], Sqrt[2] Sin[--]}
                 8        8                8

              Pi               Pi               Pi       Pi
Out[59]= {Cos[--], Sqrt[2] Cos[--], Sqrt[2] Cos[--], Cos[--]}
              8                8                8        8

             Pi
Out[58]= Cos[--]
             8

                        -Sqrt[2 + Sqrt[2]]
Out[56]= {{u -> -ArcCos[------------------]}, 
                                2
 
                  Sqrt[2 + Sqrt[2]]
>    {u -> ArcCos[-----------------]}}
                          2

Out[55]= C[1] \[Element] Integers && 
 
                    Sqrt[2] - Sqrt[8 - 4 Sqrt[2]]
>    (u == 2 ArcTan[-----------------------------] + 2 Pi C[1] || 
                            -2 + Sqrt[2]
 
                     Sqrt[2] + Sqrt[8 - 4 Sqrt[2]]
>      u == 2 ArcTan[-----------------------------] + 2 Pi C[1])
                             -2 + Sqrt[2]

Out[54]= C[1] \[Element] Integers && 
 
           -Pi                     3 Pi
>    (u == --- + 2 Pi C[1] || u == ---- + 2 Pi C[1])
            4                       4

Solve::ifun: Inverse functions are being used by Solve, so some solutions may
     not be found; use Reduce for complete solution information.

                -Pi
Out[53]= {{u -> ---}}
                 4

Sin[l-u]/Sin[l]+Sin[u]/Sin[l]*Sin[3*l]/Sin[2*l]

Out[73]= Sin[a + b]
Out[72]= Csc[l] Sin[l - u] + Csc[l] Csc[2 l] Sin[3 l] Sin[u]
                  
Sin[a+b]

Out[74]= Csc[l] Sin[l - u] + Csc[l] Csc[a l] Sin[(1 + a) l] Sin[u] == 
 
>    Csc[a l] Sin[a l + u]


Reduce[Sin[l-u]/Sin[l]+Sin[u]/Sin[l]*Sin[(a+1)*l]/Sin[a*l]==Sin[a*l+u]/Sin[a*l]]

Out[2]= (Sin[l - u] + Csc[a l] Sin[(1 + a) l] Sin[u] != 0 && 
 
                            Csc[a l] Sin[a l + u]
>      Csc[l] == -------------------------------------------) || 
                 Sin[l - u] + Csc[a l] Sin[(1 + a) l] Sin[u]
 
>    (Sin[a l + u] == 0 && Sin[u] == 0 && Sin[l - u] == 0) || 
 
>    (Sin[a l + u] == 0 && Sin[(1 + a) l] Sin[u] != 0 && 
 
>      Csc[a l] == -(Csc[(1 + a) l] Csc[u] Sin[l - u])) || 
 
>    (Sin[a l + u] != 0 && Sin[l - u] == 0 && Csc[a l] == 0) || 
 
>    (Sin[a l + u] == 0 && Sin[l - u] == 0 && Sin[u] != 0 && 
 
>      Sin[(1 + a) l] == 0)

[Calculating...]

Out[76]= (Sin[l - u] + Csc[a l] Sin[(1 + a) l] Sin[u] != 0 && 
 
                            Csc[a l] Sin[a l + u]
>      Csc[l] == -------------------------------------------) || 
                 Sin[l - u] + Csc[a l] Sin[(1 + a) l] Sin[u]
 
>    (Sin[a l + u] == 0 && Sin[u] == 0 && Sin[l - u] == 0) || 
 
>    (Sin[a l + u] == 0 && Sin[(1 + a) l] Sin[u] != 0 && 
 
>      Csc[a l] == -(Csc[(1 + a) l] Csc[u] Sin[l - u])) || 
 
>    (Sin[a l + u] != 0 && Sin[l - u] == 0 && Csc[a l] == 0) || 
 
>    (Sin[a l + u] == 0 && Sin[l - u] == 0 && Sin[u] != 0 && 
 
>      Sin[(1 + a) l] == 0)

Out[75]= Csc[l] Sin[l - u] + Csc[l] Csc[a l] Sin[(1 + a) l] Sin[u] == 
 
>    Csc[a l] Sin[a l + u]


?Integer

Integer is the head used for integers. 

t={{0,1,0,0},{1,0,1,0},{0,1,0,1},{0,0,1,0}}

Last[Sort[Eigenvalues[t]]]

Eigenvectors[t]



           1 + Sqrt[5]      1 + Sqrt[5]  1 + Sqrt[5]
Out[26]= {{-----------, {1, -----------, -----------, 1}}, 
                2                2            2

?ArcCos

t[[1,2]]

Out[29]= 1

Out[28]= {{0, 1, 0, 0}, {1, 0, 1, 0}, {0, 1, 0, 1}, {0, 0, 1, 0}}[{1, 2}]

                                  -1
ArcCos[z] gives the arc cosine cos  (z) of the complex number z. 

W=.

lhs=. removes any rules defined for lhs. 

W[a_Integer,b_Integer,c_Integer,d_Integer,u_]/; (adjMatrix[[a,b]]==1 && adjMatrix[[b,c]]==1 && adjMatrix[[c,d]]==1 && adjMatrix[[d,a]]==1) := a+b+c+d;



W[2,3,2,3]

ww[a_Integer,b_Integer]/;a>0&&a<3:=a*a;

Clear[W]

W[3,2,3,2,u]

Out[78]= 10

Out[77]= W[3, 2, 3, 2]

Out[76]= W[1, 2, 3, 2]

Out[73]= W[1, 2, 3, 2]

Out[46]= W[1, 2, 3, 2]

Out[42]= W[1, 2, 3, 2]

Out[36]= W[1, 2, 3, 2]

Out[33]= W[1, 2, 3, 2]

Out[32]= W[1, 2, 3, 4]

?KroneckerDelta

KroneckerDelta[n , n , ...] gives the Kronecker delta \[Delta]
                1   2                                         n  n  ...
                                                               1  2
    , equal to 1 if all the n  are equal, and 0 otherwise. 
                             i

DifferenceDelta DiracDelta      DiscreteDelta   KroneckerDelta

RSOSModel[adjMatrix__List]:=
	  Module[{es,ev,lm,W},
		 es=Sort[Transpose[Eigensystem[t]],#1[[1]]>#2[[1]]&][[1]];
		 lm=ArcCos[es[[1]]/2];
		 ev=es[[2]];
		 s[r_,u_]:=Sin[r*lm+u]/Sin[lm];
		 W[a_Integer,b_Integer,c_Integer,d_Integer,u_]/; (adjMatrix[[a,b]]==1 && adjMatrix[[b,c]]==1 && adjMatrix[[c,d]]==1 && adjMatrix[[d,a]]==1) := s[1,-u]*KroneckerDelta[a,c]+s[0,u]*Sqrt[ev[[a]]]*Sqrt[ev[[c]]]/ev[[b]]*KroneckerDelta[b,d];
		 W[a_Integer,b_Integer,c_Integer,d_Integer,u_]:=0;
		 W[#1,#2,#3,#4,lm/2]&]

Clear[RSOSModel]

w=RSOSModel[t]

                                lm$182
Out[101]= W$182[#1, #2, #3, #4, ------] & 
                                  2

N[w[1,2,3,4]]

Out[103]= 0.

Out[102]= 0.413304

Out[98]= 0.413304

                               lm$181
Out[96]= W$181[#1, #2, #3, #4, ------] & 
                                 2

                               lm$179
Out[86]= W$179[#1, #2, #3, #4, ------] & 
        
Clear[W]

Table[W[{{d,c},{a,b}}]==w[a,b,c,d],{a,1,4},{b,1,4},{c,1,4},{d,1,4}]

Out[104]= {{{{W[{{1, 1}, {1, 1}}] == 0., W[{{2, 1}, {1, 1}}] == 0., 
 
>       W[{{3, 1}, {1, 1}}] == 0., W[{{4, 1}, {1, 1}}] == 0.}, 
 
>      {W[{{1, 2}, {1, 1}}] == 0., W[{{2, 2}, {1, 1}}] == 0., 
 
>       W[{{3, 2}, {1, 1}}] == 0., W[{{4, 2}, {1, 1}}] == 0.}, 
 
>      {W[{{1, 3}, {1, 1}}] == 0., W[{{2, 3}, {1, 1}}] == 0., 
 
>       W[{{3, 3}, {1, 1}}] == 0., W[{{4, 3}, {1, 1}}] == 0.}, 
 
>      {W[{{1, 4}, {1, 1}}] == 0., W[{{2, 4}, {1, 1}}] == 0., 
 
>       W[{{3, 4}, {1, 1}}] == 0., W[{{4, 4}, {1, 1}}] == 0.}}, 
 
>     {{W[{{1, 1}, {1, 2}}] == 0., W[{{2, 1}, {1, 2}}] == 0.850651, 
 
>       W[{{3, 1}, {1, 2}}] == 0., W[{{4, 1}, {1, 2}}] == 0.}, 
 
>      {W[{{1, 2}, {1, 2}}] == 0., W[{{2, 2}, {1, 2}}] == 0., 
 
>       W[{{3, 2}, {1, 2}}] == 0., W[{{4, 2}, {1, 2}}] == 0.}, 
 
>      {W[{{1, 3}, {1, 2}}] == 0., W[{{2, 3}, {1, 2}}] == 0.413304, 
 
>       W[{{3, 3}, {1, 2}}] == 0., W[{{4, 3}, {1, 2}}] == 0.}, 
 
>      {W[{{1, 4}, {1, 2}}] == 0., W[{{2, 4}, {1, 2}}] == 0., 
 
>       W[{{3, 4}, {1, 2}}] == 0., W[{{4, 4}, {1, 2}}] == 0.}}, 
 
>     {{W[{{1, 1}, {1, 3}}] == 0., W[{{2, 1}, {1, 3}}] == 0., 
 
>       W[{{3, 1}, {1, 3}}] == 0., W[{{4, 1}, {1, 3}}] == 0.}, 
 
>      {W[{{1, 2}, {1, 3}}] == 0., W[{{2, 2}, {1, 3}}] == 0., 
 
>       W[{{3, 2}, {1, 3}}] == 0., W[{{4, 2}, {1, 3}}] == 0.}, 
 
>      {W[{{1, 3}, {1, 3}}] == 0., W[{{2, 3}, {1, 3}}] == 0., 
 
>       W[{{3, 3}, {1, 3}}] == 0., W[{{4, 3}, {1, 3}}] == 0.}, 
 
>      {W[{{1, 4}, {1, 3}}] == 0., W[{{2, 4}, {1, 3}}] == 0., 
 
>       W[{{3, 4}, {1, 3}}] == 0., W[{{4, 4}, {1, 3}}] == 0.}}, 
 
>     {{W[{{1, 1}, {1, 4}}] == 0., W[{{2, 1}, {1, 4}}] == 0., 
 
>       W[{{3, 1}, {1, 4}}] == 0., W[{{4, 1}, {1, 4}}] == 0.}, 
 
>      {W[{{1, 2}, {1, 4}}] == 0., W[{{2, 2}, {1, 4}}] == 0., 
 
>       W[{{3, 2}, {1, 4}}] == 0., W[{{4, 2}, {1, 4}}] == 0.}, 
 
>      {W[{{1, 3}, {1, 4}}] == 0., W[{{2, 3}, {1, 4}}] == 0., 
 
>       W[{{3, 3}, {1, 4}}] == 0., W[{{4, 3}, {1, 4}}] == 0.}, 
 
>      {W[{{1, 4}, {1, 4}}] == 0., W[{{2, 4}, {1, 4}}] == 0., 
 
>       W[{{3, 4}, {1, 4}}] == 0., W[{{4, 4}, {1, 4}}] == 0.}}}, 
 
>    {{{W[{{1, 1}, {2, 1}}] == 0., W[{{2, 1}, {2, 1}}] == 0., 
 
>       W[{{3, 1}, {2, 1}}] == 0., W[{{4, 1}, {2, 1}}] == 0.}, 
 
>      {W[{{1, 2}, {2, 1}}] == 1.37638, W[{{2, 2}, {2, 1}}] == 0., 
 
>       W[{{3, 2}, {2, 1}}] == 0.525731, W[{{4, 2}, {2, 1}}] == 0.}, 
 
>      {W[{{1, 3}, {2, 1}}] == 0., W[{{2, 3}, {2, 1}}] == 0., 
 
>       W[{{3, 3}, {2, 1}}] == 0., W[{{4, 3}, {2, 1}}] == 0.}, 
 
>      {W[{{1, 4}, {2, 1}}] == 0., W[{{2, 4}, {2, 1}}] == 0., 
 
>       W[{{3, 4}, {2, 1}}] == 0., W[{{4, 4}, {2, 1}}] == 0.}}, 
 
>     {{W[{{1, 1}, {2, 2}}] == 0., W[{{2, 1}, {2, 2}}] == 0., 
 
>       W[{{3, 1}, {2, 2}}] == 0., W[{{4, 1}, {2, 2}}] == 0.}, 
 
>      {W[{{1, 2}, {2, 2}}] == 0., W[{{2, 2}, {2, 2}}] == 0., 
 
>       W[{{3, 2}, {2, 2}}] == 0., W[{{4, 2}, {2, 2}}] == 0.}, 
 
>      {W[{{1, 3}, {2, 2}}] == 0., W[{{2, 3}, {2, 2}}] == 0., 
 
>       W[{{3, 3}, {2, 2}}] == 0., W[{{4, 3}, {2, 2}}] == 0.}, 
 
>      {W[{{1, 4}, {2, 2}}] == 0., W[{{2, 4}, {2, 2}}] == 0., 
 
>       W[{{3, 4}, {2, 2}}] == 0., W[{{4, 4}, {2, 2}}] == 0.}}, 
 
>     {{W[{{1, 1}, {2, 3}}] == 0., W[{{2, 1}, {2, 3}}] == 0., 
 
>       W[{{3, 1}, {2, 3}}] == 0., W[{{4, 1}, {2, 3}}] == 0.}, 
 
>      {W[{{1, 2}, {2, 3}}] == 0.525731, W[{{2, 2}, {2, 3}}] == 0., 
 
>       W[{{3, 2}, {2, 3}}] == 1.05146, W[{{4, 2}, {2, 3}}] == 0.}, 
 
>      {W[{{1, 3}, {2, 3}}] == 0., W[{{2, 3}, {2, 3}}] == 0., 
 
>       W[{{3, 3}, {2, 3}}] == 0., W[{{4, 3}, {2, 3}}] == 0.}, 
 
>      {W[{{1, 4}, {2, 3}}] == 0., W[{{2, 4}, {2, 3}}] == 0., 
 
>       W[{{3, 4}, {2, 3}}] == 0.413304, W[{{4, 4}, {2, 3}}] == 0.}}, 
 
>     {{W[{{1, 1}, {2, 4}}] == 0., W[{{2, 1}, {2, 4}}] == 0., 
 
>       W[{{3, 1}, {2, 4}}] == 0., W[{{4, 1}, {2, 4}}] == 0.}, 
 
>      {W[{{1, 2}, {2, 4}}] == 0., W[{{2, 2}, {2, 4}}] == 0., 
 
>       W[{{3, 2}, {2, 4}}] == 0., W[{{4, 2}, {2, 4}}] == 0.}, 
 
>      {W[{{1, 3}, {2, 4}}] == 0., W[{{2, 3}, {2, 4}}] == 0., 
 
>       W[{{3, 3}, {2, 4}}] == 0., W[{{4, 3}, {2, 4}}] == 0.}, 
 
>      {W[{{1, 4}, {2, 4}}] == 0., W[{{2, 4}, {2, 4}}] == 0., 
 
>       W[{{3, 4}, {2, 4}}] == 0., W[{{4, 4}, {2, 4}}] == 0.}}}, 
 
>    {{{W[{{1, 1}, {3, 1}}] == 0., W[{{2, 1}, {3, 1}}] == 0., 
 
>       W[{{3, 1}, {3, 1}}] == 0., W[{{4, 1}, {3, 1}}] == 0.}, 
 
>      {W[{{1, 2}, {3, 1}}] == 0., W[{{2, 2}, {3, 1}}] == 0., 
 
>       W[{{3, 2}, {3, 1}}] == 0., W[{{4, 2}, {3, 1}}] == 0.}, 
 
>      {W[{{1, 3}, {3, 1}}] == 0., W[{{2, 3}, {3, 1}}] == 0., 
 
>       W[{{3, 3}, {3, 1}}] == 0., W[{{4, 3}, {3, 1}}] == 0.}, 
 
>      {W[{{1, 4}, {3, 1}}] == 0., W[{{2, 4}, {3, 1}}] == 0., 
 
>       W[{{3, 4}, {3, 1}}] == 0., W[{{4, 4}, {3, 1}}] == 0.}}, 
 
>     {{W[{{1, 1}, {3, 2}}] == 0., W[{{2, 1}, {3, 2}}] == 0.413304, 
 
>       W[{{3, 1}, {3, 2}}] == 0., W[{{4, 1}, {3, 2}}] == 0.}, 
 
>      {W[{{1, 2}, {3, 2}}] == 0., W[{{2, 2}, {3, 2}}] == 0., 
 
>       W[{{3, 2}, {3, 2}}] == 0., W[{{4, 2}, {3, 2}}] == 0.}, 
 
>      {W[{{1, 3}, {3, 2}}] == 0., W[{{2, 3}, {3, 2}}] == 1.05146, 
 
>       W[{{3, 3}, {3, 2}}] == 0., W[{{4, 3}, {3, 2}}] == 0.525731}, 
 
>      {W[{{1, 4}, {3, 2}}] == 0., W[{{2, 4}, {3, 2}}] == 0., 
 
>       W[{{3, 4}, {3, 2}}] == 0., W[{{4, 4}, {3, 2}}] == 0.}}, 
 
>     {{W[{{1, 1}, {3, 3}}] == 0., W[{{2, 1}, {3, 3}}] == 0., 
 
>       W[{{3, 1}, {3, 3}}] == 0., W[{{4, 1}, {3, 3}}] == 0.}, 
 
>      {W[{{1, 2}, {3, 3}}] == 0., W[{{2, 2}, {3, 3}}] == 0., 
 
>       W[{{3, 2}, {3, 3}}] == 0., W[{{4, 2}, {3, 3}}] == 0.}, 
 
>      {W[{{1, 3}, {3, 3}}] == 0., W[{{2, 3}, {3, 3}}] == 0., 
 
>       W[{{3, 3}, {3, 3}}] == 0., W[{{4, 3}, {3, 3}}] == 0.}, 
 
>      {W[{{1, 4}, {3, 3}}] == 0., W[{{2, 4}, {3, 3}}] == 0., 
 
>       W[{{3, 4}, {3, 3}}] == 0., W[{{4, 4}, {3, 3}}] == 0.}}, 
 
>     {{W[{{1, 1}, {3, 4}}] == 0., W[{{2, 1}, {3, 4}}] == 0., 
 
>       W[{{3, 1}, {3, 4}}] == 0., W[{{4, 1}, {3, 4}}] == 0.}, 
 
>      {W[{{1, 2}, {3, 4}}] == 0., W[{{2, 2}, {3, 4}}] == 0., 
 
>       W[{{3, 2}, {3, 4}}] == 0., W[{{4, 2}, {3, 4}}] == 0.}, 
 
>      {W[{{1, 3}, {3, 4}}] == 0., W[{{2, 3}, {3, 4}}] == 0.525731, 
 
>       W[{{3, 3}, {3, 4}}] == 0., W[{{4, 3}, {3, 4}}] == 1.37638}, 
 
>      {W[{{1, 4}, {3, 4}}] == 0., W[{{2, 4}, {3, 4}}] == 0., 
 
>       W[{{3, 4}, {3, 4}}] == 0., W[{{4, 4}, {3, 4}}] == 0.}}}, 
 
>    {{{W[{{1, 1}, {4, 1}}] == 0., W[{{2, 1}, {4, 1}}] == 0., 
 
>       W[{{3, 1}, {4, 1}}] == 0., W[{{4, 1}, {4, 1}}] == 0.}, 
 
>      {W[{{1, 2}, {4, 1}}] == 0., W[{{2, 2}, {4, 1}}] == 0., 
 
>       W[{{3, 2}, {4, 1}}] == 0., W[{{4, 2}, {4, 1}}] == 0.}, 
 
>      {W[{{1, 3}, {4, 1}}] == 0., W[{{2, 3}, {4, 1}}] == 0., 
 
>       W[{{3, 3}, {4, 1}}] == 0., W[{{4, 3}, {4, 1}}] == 0.}, 
 
>      {W[{{1, 4}, {4, 1}}] == 0., W[{{2, 4}, {4, 1}}] == 0., 
 
>       W[{{3, 4}, {4, 1}}] == 0., W[{{4, 4}, {4, 1}}] == 0.}}, 
 
>     {{W[{{1, 1}, {4, 2}}] == 0., W[{{2, 1}, {4, 2}}] == 0., 
 
>       W[{{3, 1}, {4, 2}}] == 0., W[{{4, 1}, {4, 2}}] == 0.}, 
 
>      {W[{{1, 2}, {4, 2}}] == 0., W[{{2, 2}, {4, 2}}] == 0., 
 
>       W[{{3, 2}, {4, 2}}] == 0., W[{{4, 2}, {4, 2}}] == 0.}, 
 
>      {W[{{1, 3}, {4, 2}}] == 0., W[{{2, 3}, {4, 2}}] == 0., 
 
>       W[{{3, 3}, {4, 2}}] == 0., W[{{4, 3}, {4, 2}}] == 0.}, 
 
>      {W[{{1, 4}, {4, 2}}] == 0., W[{{2, 4}, {4, 2}}] == 0., 
 
>       W[{{3, 4}, {4, 2}}] == 0., W[{{4, 4}, {4, 2}}] == 0.}}, 
 
>     {{W[{{1, 1}, {4, 3}}] == 0., W[{{2, 1}, {4, 3}}] == 0., 
 
>       W[{{3, 1}, {4, 3}}] == 0., W[{{4, 1}, {4, 3}}] == 0.}, 
 
>      {W[{{1, 2}, {4, 3}}] == 0., W[{{2, 2}, {4, 3}}] == 0., 
 
>       W[{{3, 2}, {4, 3}}] == 0.413304, W[{{4, 2}, {4, 3}}] == 0.}, 
 
>      {W[{{1, 3}, {4, 3}}] == 0., W[{{2, 3}, {4, 3}}] == 0., 
 
>       W[{{3, 3}, {4, 3}}] == 0., W[{{4, 3}, {4, 3}}] == 0.}, 
 
>      {W[{{1, 4}, {4, 3}}] == 0., W[{{2, 4}, {4, 3}}] == 0., 
 
>       W[{{3, 4}, {4, 3}}] == 0.850651, W[{{4, 4}, {4, 3}}] == 0.}}, 
 
>     {{W[{{1, 1}, {4, 4}}] == 0., W[{{2, 1}, {4, 4}}] == 0., 
 
>       W[{{3, 1}, {4, 4}}] == 0., W[{{4, 1}, {4, 4}}] == 0.}, 
 
>      {W[{{1, 2}, {4, 4}}] == 0., W[{{2, 2}, {4, 4}}] == 0., 
 
>       W[{{3, 2}, {4, 4}}] == 0., W[{{4, 2}, {4, 4}}] == 0.}, 
 
>      {W[{{1, 3}, {4, 4}}] == 0., W[{{2, 3}, {4, 4}}] == 0., 
 
>       W[{{3, 3}, {4, 4}}] == 0., W[{{4, 3}, {4, 4}}] == 0.}, 
 
>      {W[{{1, 4}, {4, 4}}] == 0., W[{{2, 4}, {4, 4}}] == 0., 
 
>       W[{{3, 4}, {4, 4}}] == 0., W[{{4, 4}, {4, 4}}] == 0.}}}}

                                                             1 + Sqrt[5]
                                                      ArcCos
Out[83]= W$178

Out[82]= W$177

	 

c/(3*k-8)/.{c->7/10, k->5}

          1
Out[147]= --
          10

          7
Out[146]= --
          16

Out[145]= k