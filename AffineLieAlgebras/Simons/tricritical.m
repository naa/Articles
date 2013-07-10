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

w[2,3,2,3]

2*Sin[Pi/10]/Sin[Pi/5]

[Calculating...]

[Calculating...]

             -1 + Sqrt[5]
Out[274]= -------------------
                 5   Sqrt[5]
          2 Sqrt[- - -------]
                 8      8

                       1 + Sqrt[5]
                ArcCos[-----------]
                            4
          2 Sin[-------------------]
                         2
Out[273]= --------------------------
                                 2
                    (1 + Sqrt[5])
           Sqrt[1 - --------------]
                          16

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

2^24

Out[150]= 16777216

Out[149]= 4096

Out[148]= 18446744073709551616


N[Pi]

Out[152]= 3.14159

Out[151]= Pi

Cos[5*Pi/4]

               1
Out[153]= -(-------)
            Sqrt[2]

Cos[4*Pi/5]

          -1 - Sqrt[5]
Out[154]= ------------
               4


N[Sqrt[2*Sqrt[5]+2]]

Out[158]= 2.54404

?Conjugate

Clear[P]

P[theta_,z_]:=1/2*(z+Exp[2*I*theta]*Conjugate[z])

al=Sqrt[2]-1

Out[175]= -1 + Sqrt[2]


P[-3*Pi/8,al^g*Exp[-I/2*ww]]

                         g
           (-1 + Sqrt[2])                  Conjugate[g]
Out[181]= (--------------- + (-1 + Sqrt[2])             
                I/2 ww
               E
 
         (-3 I)/4 Pi + I/2 Conjugate[ww]
>       E                               ) / 2


?Reduce

Reduce[expr, vars] reduces the statement expr
     by solving equations or inequalities for vars
      and eliminating quantifiers. Reduce[expr, vars, dom]

       does the reduction over the domain dom
       . Common choices of dom are Reals, Integers and Complexes. 

N[Reduce[P[-3*Pi/8,al^g*Exp[-I/2*ww]]==P[-3*Pi/8,al^(g+1)*Exp[-Pi*I/4+-I/2*ww]]]/.{g->5,ww->2*Pi}]


Simplify[P[-3*Pi/8,al*Exp[-3*Pi*I/4]]]

Abs[z] gives the absolute value of the real or complex number z. 

Abs[z] gives the absolute value of the real or complex number z. 

                      3   I
                      - + -
                I     2   2
Out[192]= (-1 - -) + -------
                2    Sqrt[2]

                         -1 + Sqrt[2]
          -1 + Sqrt[2] + ------------
                          (3 I)/4 Pi
                         E
Out[191]= ---------------------------
                       2

Out[190]= False
                        5                    5       1/4               5
Out[189]= (-1 + Sqrt[2])  == (-(-1 + Sqrt[2])  - (-1)    (-1 + Sqrt[2])  + 
 
            1/4                       5             3/4
>       (-1)    Sqrt[2] (-1 + Sqrt[2]) ) / (I + (-1)    - I Sqrt[2])
           I/2 ww                       5
Out[188]= E       != 0 && (-1 + Sqrt[2])  == 
 
                       5  I/2 ww + I/2 Conjugate[ww]
>     (-((-1 + Sqrt[2])  E                          ) - 
 
             1/4               5  I/2 ww + I/2 Conjugate[ww]
>        (-1)    (-1 + Sqrt[2])  E                           + 
 
             1/4                       5  I/2 ww + I/2 Conjugate[ww]
>        (-1)    Sqrt[2] (-1 + Sqrt[2])  E                          ) / 
 
                3/4
>      (I + (-1)    - I Sqrt[2])           I/2 ww                       g
Out[187]= E       != 0 && (-1 + Sqrt[2])  == 
 
                       Conjugate[g]  I/2 ww + I/2 Conjugate[ww]
>     (-((-1 + Sqrt[2])             E                          ) - 
 
             1/4               Conjugate[g]  I/2 ww + I/2 Conjugate[ww]
>        (-1)    (-1 + Sqrt[2])             E                           + 
 
             1/4                       Conjugate[g]
>        (-1)    Sqrt[2] (-1 + Sqrt[2])             
 
           I/2 ww + I/2 Conjugate[ww]             3/4
>         E                          ) / (I + (-1)    - I Sqrt[2])

Out[186]= $Aborted

                                        g
                          (-1 + Sqrt[2])
Reduce::nddc: The system (--------------- + 
                               I/2 ww
                              E
                      Conjugate[g]  (-3 I)/4 Pi + I/2 Conjugate[ww]
        (-1 + Sqrt[2])             E                               ) / 2 == 
                                       -I
     <<1>> contains a nonreal constant -- Pi. With the domain Reals
                                       4
     specified, all constants should be real.

                                g
                  (-1 + Sqrt[2])
Out[185]= Reduce[(--------------- + 
                       I/2 ww
                      E
 
                       Conjugate[g]  (-3 I)/4 Pi + I/2 Conjugate[ww]
>        (-1 + Sqrt[2])             E                               ) / 2 == 
 
                     1 + g  -I/4 Pi - I/2 ww
>     ((-1 + Sqrt[2])      E                 + 
 
                       1 + Conjugate[g]  -I/2 Pi + I/2 Conjugate[ww]
>        (-1 + Sqrt[2])                 E                           ) / 2, g, 
 
>    Reals]

           I/2 ww                       g
Out[183]= E       != 0 && (-1 + Sqrt[2])  == 
 
                       Conjugate[g]  I/2 ww + I/2 Conjugate[ww]
>     (-((-1 + Sqrt[2])             E                          ) - 
 
             1/4               Conjugate[g]  I/2 ww + I/2 Conjugate[ww]
>        (-1)    (-1 + Sqrt[2])             E                           + 
 
             1/4                       Conjugate[g]
>        (-1)    Sqrt[2] (-1 + Sqrt[2])             
 
           I/2 ww + I/2 Conjugate[ww]             3/4
>         E                          ) / (I + (-1)    - I Sqrt[2])

                         1 + g  -I/4 Pi - I/2 ww
Out[182]= ((-1 + Sqrt[2])      E                 + 
 
                     1 + Conjugate[g]  -I/2 Pi + I/2 Conjugate[ww]
>      (-1 + Sqrt[2])                 E                           ) / 2

Out[180]= 
 
                   g
     (-1 + Sqrt[2])        1/4               Conjugate[g]  I/2 Conjugate[ww]
     --------------- - (-1)    (-1 + Sqrt[2])             E
          I/2 ww
         E
>    -----------------------------------------------------------------------
                                        2

                         g
           (-1 + Sqrt[2])                  Conjugate[g]
Out[179]= (--------------- + (-1 + Sqrt[2])             
                I/2 ww
               E
 
         (-3 I)/4 Pi + I/2 Conjugate[ww]
>       E                               ) / 2

                             g
               (-1 + Sqrt[2])
Reduce::naqs: (--------------- + 
                    I/2 ww
                   E
                     Conjugate[g]  (-3 I)/4 Pi + I/2 Conjugate[ww]
       (-1 + Sqrt[2])             E                               ) / 2 is not
     a quantified system of equations and inequalities.

                                g
                  (-1 + Sqrt[2])
Out[178]= Reduce[(--------------- + 
                       I/2 ww
                      E
 
                      Conjugate[g]  (-3 I)/4 Pi + I/2 Conjugate[ww]
>       (-1 + Sqrt[2])             E                               ) / 2]


Abs[I]

Out[196]= 1

Out[195]= 1

pr[x_,u_]:=1/2*(x+ Conjugate[x]*u^2/Abs[u]^2)


pr[x*y,u]

pr[x,u]*pr[y,u]


?Reduce

Reduce::ztest: Unable to decide whether numeric quantities 
                                                                      Pi
    {1 - Sqrt[2] + Sqrt[5] - Sqrt[10] - Sqrt[2 (5 - Sqrt[5])] + 4 Cos[--] + 
                                                                      20
            Pi                  Pi                             Pi
      4 Sin[--] - 4 Sqrt[2] Sin[--], 2 + <<6>> - 4 Sqrt[2] Sin[--], 
            20                  20                             20
                                     Pi 2       Pi 2
     -2 + <<7>>, <<2>>, 16 (-1 + Cos[--]  + Sin[--] )} are equal to zero.
                                     20         20
     Assuming they are.

                                                                    Pi
                   -2 Sqrt[2] - Re[x] - Sqrt[5] Re[x] - 4 Re[x] Sin[--]
                                                                    20
Out[247]= Im[x] == ----------------------------------------------------
                              Pi                   2     4
                        4 Cos[--] + Root[80 - 20 #1  + #1  & , 3]
                              20

Reduce[pr[1,Exp[-3*Pi*I/8]]==pr[x*Exp[I*3/2*(Pi/2)],Exp[-3*Pi*I/8]]]


Reduce[pr[x,Exp[-3*Pi*I/8]]==pr[x*Exp[I*3/5*(Pi/2)],Exp[-3*Pi*I/8]]]

                       1 + Sqrt[5]
                ArcCos[-----------]
                            4
          2 Sin[-------------------]
                         2
Out[272]= --------------------------
                                 2
                    (1 + Sqrt[5])
           Sqrt[1 - --------------]
                          16

                       1 + Sqrt[5]
                ArcCos[-----------]
                            4
          2 Sin[-------------------]
                         2
Out[271]= --------------------------
                                 2
                    (1 + Sqrt[5])
           Sqrt[1 - --------------]
                          16

Reduce::ztest: Unable to decide whether numeric quantities 
           Pi      Pi         Pi        Pi
    {-(Cos[--] Cot[--]) + Csc[--] - Sin[--], 
           20      20         20        20
         Pi      Pi        Pi        Pi
     Cos[--] Cot[--] - Csc[--] + Sin[--]} are equal to zero. Assuming they
         20      20        20        20
     are.

                               Pi
                     Re[x] Sin[--]
                               20
Out[270]= Im[x] == -(-------------)
                              Pi
                     -1 + Cos[--]
                              20

Reduce::ztest: Unable to decide whether numeric quantities 
           Pi      Pi         Pi        Pi
    {-(Cos[--] Cot[--]) + Csc[--] - Sin[--], 
           20      20         20        20
         Pi      Pi        Pi        Pi
     Cos[--] Cot[--] - Csc[--] + Sin[--]} are equal to zero. Assuming they
         20      20        20        20
     are.

Solve::fulldim: 
   The solution set contains a full-dimensional component; use Reduce for
    complete solution information.

Out[269]= {{}}

                         x     (3 I)/10 Pi + x
Solve::nddc: The system E  == E                - 
       (-3 I)/4 Pi + Conjugate[x]    (-21 I)/20 Pi + Conjugate[x]
      E                           + E                             contains a
                      3 I
     nonreal constant --- Pi. With the domain Reals
                      10
     specified, all constants should be real.

                 x     (3 I)/10 Pi + x    (-3 I)/4 Pi + Conjugate[x]
Out[268]= Solve[E  == E                - E                           + 
 
        (-21 I)/20 Pi + Conjugate[x]
>      E                            , x, Reals]

Solve::nsmet: This system cannot be solved with the methods available to
    Solve.

                 x     (3 I)/10 Pi + x    (-3 I)/4 Pi + Conjugate[x]
Out[267]= Solve[E  == E                - E                           + 
 
        (-21 I)/20 Pi + Conjugate[x]
>      E                            , x]

           x     (3 I)/10 Pi + x    (-3 I)/4 Pi + Conjugate[x]
Out[266]= E  == E                - E                           + 
 
       (-21 I)/20 Pi + Conjugate[x]
>     E

           x     (-3 I)/10 Pi + x    (-9 I)/20 Pi + Conjugate[x]
Out[265]= E  == E                 + E                            - 
 
       (-3 I)/4 Pi + Conjugate[x]
>     E




Reduce::ztest: Unable to decide whether numeric quantities 
                                                                      Pi
    {1 + Sqrt[2] + Sqrt[5] + Sqrt[10] + Sqrt[2 (5 - Sqrt[5])] - 4 Cos[--] - 
                                                                      20
                    Pi          Pi                                Pi
      4 Sqrt[2] Cos[--] - 4 Sin[--], -1 - Sqrt[2] + <<5>> + 4 Sin[--], <<3>>, 
                    20          20                                20
                  Pi 2       Pi 2
     16 (-1 + Cos[--]  + Sin[--] )} are equal to zero. Assuming they are.
                  20         20

                                                             Pi
                   2 Sqrt[2] + Re[x] + Sqrt[5] Re[x] - 4 Cos[--] Re[x]
                                                             20
Out[264]= Im[x] == ---------------------------------------------------
                                       2     4                Pi
                        Root[80 - 20 #1  + #1  & , 3] + 4 Sin[--]
                                                              20

Reduce::ztest: Unable to decide whether numeric quantities 
                                                                       Pi
    {-1 + Sqrt[2] - Sqrt[5] + Sqrt[10] - Sqrt[2 (5 - Sqrt[5])] - 4 Cos[--] + 
                                                                       20
                    Pi          Pi                               Pi
      4 Sqrt[2] Cos[--] - 4 Sin[--], 1 - Sqrt[2] + <<5>> + 4 Sin[--], <<3>>, 
                    20          20                               20
                  Pi 2       Pi 2
     16 (-1 + Cos[--]  + Sin[--] )} are equal to zero. Assuming they are.
                  20         20

                                                              Pi
                   -2 Sqrt[2] + Re[x] + Sqrt[5] Re[x] + 4 Cos[--] Re[x]
                                                              20
Out[263]= Im[x] == ----------------------------------------------------
                                       2     4                Pi
                        Root[80 - 20 #1  + #1  & , 3] - 4 Sin[--]
                                                              20

Solve[Out[247][[2]][[2]]==0,x]

1/N[-2*Sqrt[2]/(Sqrt[5]+1+4*Sin[Pi/20])]


Pi/N[ArcCos[(1+Sqrt[5])/4]]

Cos[Pi/5]

                   Sqrt[2] + 2 Re[x] + Sqrt[2] Re[x]
Out[262]= Im[x] == ---------------------------------
                                Sqrt[2]

                   Sqrt[2] + 2 Re[x] + Sqrt[2] Re[x]
Out[261]= Im[x] == ---------------------------------
                                Sqrt[2]



          1 + Sqrt[5]
Out[260]= -----------
               4

Out[259]= 5.

Out[258]= 0.628319

                 1 + Sqrt[5]
Out[257]= ArcCos[-----------]
                      4

ArcCos Cos

Information::notfound: Symbol Acos not found.

Out[254]= -1.36535

Out[253]= -0.73241

Solve::fulldim: 
   The solution set contains a full-dimensional component; use Reduce for
    complete solution information.

Out[252]= {{}}

                                                           Pi
Out[251]= -2 Sqrt[2] - Re[x] - Sqrt[5] Re[x] - 4 Re[x] Sin[--]
                                                           20

                              1
Out[250]= -----------------------------------------
                Pi                   2     4
          4 Cos[--] + Root[80 - 20 #1  + #1  & , 3]
                20

                                                           Pi
          -2 Sqrt[2] - Re[x] - Sqrt[5] Re[x] - 4 Re[x] Sin[--]
                                                           20
Out[249]= ----------------------------------------------------
                     Pi                   2     4
               4 Cos[--] + Root[80 - 20 #1  + #1  & , 3]
                     20

[Calculating...]

Reduce[expr, vars] reduces the statement expr
     by solving equations or inequalities for vars
      and eliminating quantifiers. Reduce[expr, vars, dom]

       does the reduction over the domain dom
       . Common choices of dom are Reals, Integers and Complexes. 

Reduce::ztest: Unable to decide whether numeric quantities 
                                                                      Pi
    {1 - Sqrt[2] + Sqrt[5] - Sqrt[10] - Sqrt[2 (5 - Sqrt[5])] + 4 Cos[--] + 
                                                                      20
            Pi                  Pi                             Pi
      4 Sin[--] - 4 Sqrt[2] Sin[--], 2 + <<6>> - 4 Sqrt[2] Sin[--], 
            20                  20                             20
                                     Pi 2       Pi 2
     -2 + <<7>>, <<2>>, 16 (-1 + Cos[--]  + Sin[--] )} are equal to zero.
                                     20         20
     Assuming they are.

                                                                    Pi
                   -2 Sqrt[2] - Re[x] - Sqrt[5] Re[x] - 4 Re[x] Sin[--]
                                                                    20
Out[244]= Im[x] == ----------------------------------------------------
                              Pi                   2     4
                        4 Cos[--] + Root[80 - 20 #1  + #1  & , 3]
                              20

Solve::fulldim: 
   The solution set contains a full-dimensional component; use Reduce for
    complete solution information.

Solve::ztest: Unable to decide whether numeric quantities 
         Pi        (3 I)/10 Pi                (3 I)/10 Pi
    {Cos[--] + Im[E           ] - Sqrt[2] Im[E           ] - 
         20
          (3 I)/10 Pi        Pi                Pi
      Re[E           ] + Sin[--] - Sqrt[2] Sin[--], <<4>>, 
                             20                20
         Pi 2       (3 I)/10 Pi 2       (3 I)/10 Pi 2       Pi 2
     Cos[--]  - Im[E           ]  - Re[E           ]  + Sin[--] } are equal to
         20                                                 20
     zero. Assuming they are.

Out[243]= {{}}

Reduce::ztest: Unable to decide whether numeric quantities 
                                                                      Pi
    {1 - Sqrt[2] + Sqrt[5] - Sqrt[10] - Sqrt[2 (5 - Sqrt[5])] + 4 Cos[--] + 
                                                                      20
            Pi                  Pi                             Pi
      4 Sin[--] - 4 Sqrt[2] Sin[--], 2 + <<6>> - 4 Sqrt[2] Sin[--], 
            20                  20                             20
                                     Pi 2       Pi 2
     -2 + <<7>>, <<2>>, 16 (-1 + Cos[--]  + Sin[--] )} are equal to zero.
                                     20         20
     Assuming they are.

                                                                    Pi
                   -2 Sqrt[2] - Re[x] - Sqrt[5] Re[x] - 4 Re[x] Sin[--]
                                                                    20
Out[242]= Im[x] == ----------------------------------------------------
                              Pi                   2     4
                        4 Cos[--] + Root[80 - 20 #1  + #1  & , 3]
                              20

                   Sqrt[2] + 2 Re[x] + Sqrt[2] Re[x]
Out[241]= Im[x] == ---------------------------------
                                Sqrt[2]

Solve::fulldim: 
   The solution set contains a full-dimensional component; use Reduce for
    complete solution information.

Out[240]= {{}}

Out[239]= {{x -> 
 
               (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>      (-I Im[E                                   ] 
 
               (-3 I)/2 Conjugate[ww]
>          Im[E                      ] + 
 
              (-3 I)/2 Conjugate[ww]      (3 I)/2 ww
>         Im[E                      ] Re[E          ] + 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>         I Im[E                                   ] 
 
               (3 I)/4 (Pi + 2 ww)
>          Re[E                   ] + 
 
              (3 I)/2 ww      (3 I)/4 (Pi + 2 ww)
>         Re[E          ] Re[E                   ] + 
 
              (-3 I)/2 Conjugate[ww]
>         Im[E                      ] 
 
               (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>          Re[E                                   ] + 
 
              (3 I)/4 (Pi + 2 ww)      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>         Re[E                   ] Re[E                                   ] + 
 
              (3 I)/4 (Pi + 2 ww)
>         Im[E                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>          (Im[E                                   ] - 
 
                    (3 I)/2 ww        (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>            I (Re[E          ] + Re[E                                   ]))\
 
                 (3 I)/2 ww       (3 I)/4 (Pi + 2 ww)
>          + Im[E          ] (Im[E                   ] - 
 
                   (-3 I)/2 Conjugate[ww]          (3 I)/4 (Pi + 2 ww)
>            I Im[E                      ] + I Re[E                   ] - 
 
                 (-3 I)/2 Conjugate[ww]
>            Re[E                      ]) - 
 
              (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>         Im[E                                   ] 
 
               (-3 I)/2 Conjugate[ww]
>          Re[E                      ] - 
 
                (3 I)/2 ww      (-3 I)/2 Conjugate[ww]
>         I Re[E          ] Re[E                      ] - 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>         I Re[E                                   ] 
 
               (-3 I)/2 Conjugate[ww]
>          Re[E                      ]) / 
 
             (3 I)/4 (Pi + 2 ww) 2       (-3 I)/2 Conjugate[ww] 2
>       (Im[E                   ]  - Im[E                      ]  + 
 
              (3 I)/4 (Pi + 2 ww) 2       (-3 I)/2 Conjugate[ww] 2
>         Re[E                   ]  - Re[E                      ] )}}

Out[238]= {{x -> 
 
              (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>      (-(Im[E          ] Im[E                   ]) - 
 
               (3 I)/2 (Pi/2 + ww)
>          Im[E                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>           Im[E                                   ] + 
 
               (3 I)/2 ww
>          Im[E          ] Im[E
 
              (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                                                         ] + 
 
               (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>          Im[E                                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>           Im[E                                            ] - 
 
               (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>          Re[E          ] Re[E                   ] - 
 
               (3 I)/2 (Pi/2 + ww)
>          Re[E                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>           Re[E                                   ] + 
 
               (3 I)/2 ww
>          Re[E          ] Re[E
 
              (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                                                         ] + 
 
               (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>          Re[E                                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>           Re[E                                            ]) / 
 
               (3 I)/2 (Pi/2 + ww) 2
>        (-Im[E                   ]  + 
 
               (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>          Im[E                                            ]  - 
 
               (3 I)/2 (Pi/2 + ww) 2
>          Re[E                   ]  + 
 
               (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>          Re[E                                            ] ) + 
 
                 (3 I)/2 ww        (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>       (I (-Im[E          ] - Im[E                                   ] + 
 
                  (3 I)/2 (Pi/2 + ww)
>            (Im[E                   ] 
 
                       (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>               (-(Im[E          ] Im[E                   ]) - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Im[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Im[E                                   ] + 
 
                      (3 I)/2 ww
>                 Im[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Im[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] - 
 
                      (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>                 Re[E          ] Re[E                   ] - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Re[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Re[E                                   ] + 
 
                      (3 I)/2 ww
>                 Re[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Re[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ])) / 
 
                    (3 I)/2 (Pi/2 + ww) 2
>             (-Im[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Im[E                                            ]  - 
 
                    (3 I)/2 (Pi/2 + ww) 2
>               Re[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Re[E                                            ] ) + 
 
                  (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>            (Im[E                                            ] 
 
                       (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>               (-(Im[E          ] Im[E                   ]) - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Im[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Im[E                                   ] + 
 
                      (3 I)/2 ww
>                 Im[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Im[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] - 
 
                      (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>                 Re[E          ] Re[E                   ] - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Re[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Re[E                                   ] + 
 
                      (3 I)/2 ww
>                 Re[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Re[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ])) / 
 
                    (3 I)/2 (Pi/2 + ww) 2
>             (-Im[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Im[E                                            ]  - 
 
                    (3 I)/2 (Pi/2 + ww) 2
>               Re[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Re[E                                            ] ))) / 
 
               (3 I)/2 (Pi/2 + ww)
>        (-Re[E                   ] + 
 
               (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>          Re[E                                            ])}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[237]= {{x -> Indeterminate}}

Out[236]= {{x -> 
 
             Pi 2       Pi      3 Pi        3 Pi      Pi        Pi 2
>      (-Cos[--]  + Cos[--] Cos[----] - Cos[----] Sin[--] + Sin[--]  - 
             16         16       16          16       16        16
 
               Pi      3 Pi          3 Pi      3 Pi        Pi      3 Pi
>          Cos[--] Sin[----] + 2 Cos[----] Sin[----] - Sin[--] Sin[----]) / 
               16       16            16        16         16       16
 
               Pi 2       3 Pi 2       Pi 2       3 Pi 2
>        (-Cos[--]  + Cos[----]  - Sin[--]  + Sin[----] ) + 
               16          16          16          16
 
                Pi        3 Pi
>       (I (Sin[--] - Sin[----] + 
                16         16
 
                  3 Pi        Pi 2       Pi      3 Pi        3 Pi      Pi
>            (Cos[----] (-Cos[--]  + Cos[--] Cos[----] - Cos[----] Sin[--] + 
                   16         16         16       16          16       16
 
                      Pi 2       Pi      3 Pi          3 Pi      3 Pi
>                 Sin[--]  - Cos[--] Sin[----] + 2 Cos[----] Sin[----] - 
                      16         16       16            16        16
 
                      Pi      3 Pi
>                 Sin[--] Sin[----])) / 
                      16       16
 
                    Pi 2       3 Pi 2       Pi 2       3 Pi 2
>             (-Cos[--]  + Cos[----]  - Sin[--]  + Sin[----] ) + 
                    16          16          16          16
 
                  Pi        Pi 2       Pi      3 Pi        3 Pi      Pi
>            (Sin[--] (-Cos[--]  + Cos[--] Cos[----] - Cos[----] Sin[--] + 
                  16        16         16       16          16       16
 
                      Pi 2       Pi      3 Pi          3 Pi      3 Pi
>                 Sin[--]  - Cos[--] Sin[----] + 2 Cos[----] Sin[----] - 
                      16         16       16            16        16
 
                      Pi      3 Pi
>                 Sin[--] Sin[----])) / 
                      16       16
 
                    Pi 2       3 Pi 2       Pi 2       3 Pi 2
>             (-Cos[--]  + Cos[----]  - Sin[--]  + Sin[----] ))) / 
                    16          16          16          16
 
              Pi        3 Pi
>        (Cos[--] + Sin[----])}}
              16         16

                                 1
Power::infy: Infinite expression - encountered.
                                 0

N::meprec: Internal precision limit $MaxExtraPrecision = 50.
                                     Pi 2         Pi      Pi          Pi 2
     reached while evaluating -2 Cos[--]  + 4 Cos[--] Sin[--] + 2 Sin[--] .
                                     8            8       8           8

                                 1
Power::infy: Infinite expression - encountered.
                                 0

N::meprec: Internal precision limit $MaxExtraPrecision = 50.
                                     Pi 2         Pi      Pi          Pi 2
     reached while evaluating -2 Cos[--]  + 4 Cos[--] Sin[--] + 2 Sin[--] .
                                     8            8       8           8

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

N::meprec: Internal precision limit $MaxExtraPrecision = 50.
                                     Pi 2         Pi      Pi          Pi 2
     reached while evaluating -2 Cos[--]  + 4 Cos[--] Sin[--] + 2 Sin[--] .
                                     8            8       8           8

General::stop: Further output of N::meprec
     will be suppressed during this calculation.

Infinity::indet: 
                                 Pi
   Indeterminate expression -Cos[--] + ComplexInfinity + ComplexInfinity - 
                                 8
         Pi
     Sin[--] encountered.
         8

Out[235]= {{x -> Indeterminate}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[234]= {{x -> Indeterminate}}

Out[233]= {{x -> 
 
                            Pi                                  Pi
                     3 (1 + --)                          3 (1 + --)
              3             2     Pi         3   Pi             2     Pi
>      (-(Cos[-] Cos[---------- - --]) + Cos[- - --] Cos[---------- - --] - 
              2          2        4          2   4           2        4
 
                             Pi                            Pi
                      3 (1 + --)                    3 (1 + --)
               3             2          3   Pi             2
>          Cos[-] Cos[----------] + Cos[- - --] Cos[----------] + 
               2          2             2   4           2
 
                             Pi                                 Pi
                      3 (1 + --)                         3 (1 + --)
               3             2     Pi        3   Pi             2     Pi
>          Sin[-] Sin[---------- - --] + Sin[- - --] Sin[---------- - --] - 
               2          2        4         2   4           2        4
 
                             Pi                            Pi
                      3 (1 + --)                    3 (1 + --)
               3             2          3   Pi             2
>          Sin[-] Sin[----------] - Sin[- - --] Sin[----------]) / 
               2          2             2   4           2
 
                     Pi                      Pi                 Pi
              3 (1 + --)              3 (1 + --)         3 (1 + --)
                     2     Pi 2              2   2              2     Pi 2
>        (Cos[---------- - --]  - Cos[----------]  + Sin[---------- - --]  - 
                  2        4              2                  2        4
 
                      Pi
               3 (1 + --)
                      2   2
>          Sin[----------] ) + (I 
                   2
 
                 3        3   Pi
>          (-Sin[-] - Sin[- - --] + 
                 2        2   4
 
                         Pi
                  3 (1 + --)
                         2     Pi
>            (Sin[---------- - --] 
                      2        4
 
                                     Pi
                              3 (1 + --)
                       3             2     Pi
>               (-(Cos[-] Cos[---------- - --]) + 
                       2          2        4
 
                                         Pi                            Pi
                                  3 (1 + --)                    3 (1 + --)
                      3   Pi             2     Pi        3             2
>                 Cos[- - --] Cos[---------- - --] - Cos[-] Cos[----------] + 
                      2   4           2        4         2          2
 
                                         Pi                       Pi
                                  3 (1 + --)               3 (1 + --)
                      3   Pi             2          3             2     Pi
>                 Cos[- - --] Cos[----------] + Sin[-] Sin[---------- - --] + 
                      2   4           2             2          2        4
 
                                         Pi                            Pi
                                  3 (1 + --)                    3 (1 + --)
                      3   Pi             2     Pi        3             2
>                 Sin[- - --] Sin[---------- - --] - Sin[-] Sin[----------] - 
                      2   4           2        4         2          2
 
                                         Pi
                                  3 (1 + --)
                      3   Pi             2
>                 Sin[- - --] Sin[----------])) / 
                      2   4           2
 
                          Pi                      Pi
                   3 (1 + --)              3 (1 + --)
                          2     Pi 2              2   2
>             (Cos[---------- - --]  - Cos[----------]  + 
                       2        4              2
 
                           Pi                      Pi
                    3 (1 + --)              3 (1 + --)
                           2     Pi 2              2   2
>               Sin[---------- - --]  - Sin[----------] ) + 
                        2        4              2
 
                         Pi                        Pi
                  3 (1 + --)                3 (1 + --)
                         2           3             2     Pi
>            (Sin[----------] (-(Cos[-] Cos[---------- - --]) + 
                      2              2          2        4
 
                                         Pi                            Pi
                                  3 (1 + --)                    3 (1 + --)
                      3   Pi             2     Pi        3             2
>                 Cos[- - --] Cos[---------- - --] - Cos[-] Cos[----------] + 
                      2   4           2        4         2          2
 
                                         Pi                       Pi
                                  3 (1 + --)               3 (1 + --)
                      3   Pi             2          3             2     Pi
>                 Cos[- - --] Cos[----------] + Sin[-] Sin[---------- - --] + 
                      2   4           2             2          2        4
 
                                         Pi                            Pi
                                  3 (1 + --)                    3 (1 + --)
                      3   Pi             2     Pi        3             2
>                 Sin[- - --] Sin[---------- - --] - Sin[-] Sin[----------] - 
                      2   4           2        4         2          2
 
                                         Pi
                                  3 (1 + --)
                      3   Pi             2
>                 Sin[- - --] Sin[----------])) / 
                      2   4           2
 
                          Pi                      Pi
                   3 (1 + --)              3 (1 + --)
                          2     Pi 2              2   2
>             (Cos[---------- - --]  - Cos[----------]  + 
                       2        4              2
 
                           Pi                      Pi
                    3 (1 + --)              3 (1 + --)
                           2     Pi 2              2   2
>               Sin[---------- - --]  - Sin[----------] ))) / 
                        2        4              2
 
                      Pi                     Pi
               3 (1 + --)             3 (1 + --)
                      2     Pi               2
>        (-Cos[---------- - --] - Cos[----------])}}
                   2        4             2

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                          0 ComplexInfinity
Infinity::indet: Indeterminate expression ----------------- encountered.
                                               Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.
p
Out[232]= {{x -> Indeterminate}}

Out[231]= {{x -> 
 
              (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>      (-(Im[E          ] Im[E                   ]) - 
 
               (3 I)/2 (Pi/2 + ww)
>          Im[E                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>           Im[E                                   ] + 
 
               (3 I)/2 ww
>          Im[E          ] Im[E
 
              (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                                                         ] + 
 
               (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>          Im[E                                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>           Im[E                                            ] - 
 
               (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>          Re[E          ] Re[E                   ] - 
 
               (3 I)/2 (Pi/2 + ww)
>          Re[E                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>           Re[E                                   ] + 
 
               (3 I)/2 ww
>          Re[E          ] Re[E
 
              (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                                                         ] + 
 
               (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>          Re[E                                   ] 
 
                (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>           Re[E                                            ]) / 
 
               (3 I)/2 (Pi/2 + ww) 2
>        (-Im[E                   ]  + 
 
               (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>          Im[E                                            ]  - 
 
               (3 I)/2 (Pi/2 + ww) 2
>          Re[E                   ]  + 
 
               (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>          Re[E                                            ] ) + 
 
                 (3 I)/2 ww        (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>       (I (-Im[E          ] - Im[E                                   ] + 
 
                  (3 I)/2 (Pi/2 + ww)
>            (Im[E                   ] 
 
                       (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>               (-(Im[E          ] Im[E                   ]) - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Im[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Im[E                                   ] + 
 
                      (3 I)/2 ww
>                 Im[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Im[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] - 
 
                      (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>                 Re[E          ] Re[E                   ] - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Re[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Re[E                                   ] + 
 
                      (3 I)/2 ww
>                 Re[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Re[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ])) / 
 
                    (3 I)/2 (Pi/2 + ww) 2
>             (-Im[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Im[E                                            ]  - 
 
                    (3 I)/2 (Pi/2 + ww) 2
>               Re[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Re[E                                            ] ) + 
 
                  (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>            (Im[E                                            ] 
 
                       (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>               (-(Im[E          ] Im[E                   ]) - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Im[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Im[E                                   ] + 
 
                      (3 I)/2 ww
>                 Im[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Im[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] - 
 
                      (3 I)/2 ww      (3 I)/2 (Pi/2 + ww)
>                 Re[E          ] Re[E                   ] - 
 
                      (3 I)/2 (Pi/2 + ww)
>                 Re[E                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                  Re[E                                   ] + 
 
                      (3 I)/2 ww
>                 Re[E          ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ] + 
 
                      (-3 I)/4 Pi - (3 I)/2 Conjugate[ww]
>                 Re[E                                   ] 
 
                       (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ])) / 
 
                    (3 I)/2 (Pi/2 + ww) 2
>             (-Im[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Im[E                                            ]  - 
 
                    (3 I)/2 (Pi/2 + ww) 2
>               Re[E                   ]  + 
 
                    (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww]) 2
>               Re[E                                            ] ))) / 
 
               (3 I)/2 (Pi/2 + ww)
>        (-Re[E                   ] + 
 
               (-3 I)/4 Pi - (3 I)/2 (Pi/2 + Conjugate[ww])
>          Re[E                                            ])}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[230]= {{x -> Indeterminate}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

                                          0 ComplexInfinity
Infinity::indet: Indeterminate expression ----------------- encountered.
                                               Sqrt[2]

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[229]= {{x -> Indeterminate}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[228]= {{x -> Indeterminate}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[227]= {{x -> Indeterminate}}

2*2

Out[226]= {{x -> 
 
                      (0. - 1.5 I) ww
>      (-1. Im[2.71828               ] 
 
                      (0. - 1.5 I) (1.5708 + ww)
>           Im[2.71828                          ] - 
 
                        (0. - 1.5 I) (1.5708 + ww)
>          1. Im[2.71828                          ] 
 
                      (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>           Im[2.71828                                             ] + 
 
                     (0. - 1.5 I) ww
>          Im[2.71828               ] 
 
>           Im[2.71828
 
              (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                     ] + 
 
                     (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>          Im[2.71828                                             ] 
 
>           Im[2.71828
 
              (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                     ] - 
 
                        (0. - 1.5 I) ww
>          1. Re[2.71828               ] 
 
                      (0. - 1.5 I) (1.5708 + ww)
>           Re[2.71828                          ] - 
 
                        (0. - 1.5 I) (1.5708 + ww)
>          1. Re[2.71828                          ] 
 
                      (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>           Re[2.71828                                             ] + 
 
                     (0. - 1.5 I) ww
>          Re[2.71828               ] 
 
>           Re[2.71828
 
              (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                     ] + 
 
                     (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>          Re[2.71828                                             ] 
 
>           Re[2.71828
 
              (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                     ]) / 
 
                        (0. - 1.5 I) (1.5708 + ww) 2
>        (-1. Im[2.71828                          ]  + 
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>          Im[2.71828                                                        ]
 
            2                (0. - 1.5 I) (1.5708 + ww) 2
>             - 1. Re[2.71828                          ]  + 
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>          Re[2.71828                                                        ]
 
            2                                (0. - 1.5 I) ww
>            ) + ((0. + 1. I) (-1. Im[2.71828               ] - 
 
                          (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>            1. Im[2.71828                                             ] + 
 
                        (0. - 1.5 I) (1.5708 + ww)
>            (Im[2.71828                          ] 
 
                               (0. - 1.5 I) ww
>               (-1. Im[2.71828               ] 
 
                             (0. - 1.5 I) (1.5708 + ww)
>                  Im[2.71828                          ] - 
 
                               (0. - 1.5 I) (1.5708 + ww)
>                 1. Im[2.71828                          ] 
 
                             (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                  Im[2.71828                                             ] + 
 
                            (0. - 1.5 I) ww
>                 Im[2.71828               ] 
 
>                  Im[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
>                   + Im[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                                                                 ] 
 
>                  Im[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
                                   (0. - 1.5 I) ww
>                   - 1. Re[2.71828               ] 
 
                             (0. - 1.5 I) (1.5708 + ww)
>                  Re[2.71828                          ] - 
 
                               (0. - 1.5 I) (1.5708 + ww)
>                 1. Re[2.71828                          ] 
 
                             (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                  Re[2.71828                                             ] + 
 
                            (0. - 1.5 I) ww
>                 Re[2.71828               ] 
 
>                  Re[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
>                   + Re[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                                                                 ] 
 
>                  Re[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
>                 )) / 
 
                             (0. - 1.5 I) (1.5708 + ww) 2
>             (-1. Im[2.71828                          ]  + 
 
>               Im[2.71828
 
                   (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww]) 2
>                                                                          ] \
 
                                (0. - 1.5 I) (1.5708 + ww) 2
>                - 1. Re[2.71828                          ]  + 
 
>               Re[2.71828
 
                   (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww]) 2
>                                                                          ] )
 
>              + (Im[2.71828
 
                  (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                         ] 
 
                               (0. - 1.5 I) ww
>               (-1. Im[2.71828               ] 
 
                             (0. - 1.5 I) (1.5708 + ww)
>                  Im[2.71828                          ] - 
 
                               (0. - 1.5 I) (1.5708 + ww)
>                 1. Im[2.71828                          ] 
 
                             (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                  Im[2.71828                                             ] + 
 
                            (0. - 1.5 I) ww
>                 Im[2.71828               ] 
 
>                  Im[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
>                   + Im[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                                                                 ] 
 
>                  Im[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
                                   (0. - 1.5 I) ww
>                   - 1. Re[2.71828               ] 
 
                             (0. - 1.5 I) (1.5708 + ww)
>                  Re[2.71828                          ] - 
 
                               (0. - 1.5 I) (1.5708 + ww)
>                 1. Re[2.71828                          ] 
 
                             (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                  Re[2.71828                                             ] + 
 
                            (0. - 1.5 I) ww
>                 Re[2.71828               ] 
 
>                  Re[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
>                   + Re[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) Conjugate[ww]
>                                                                 ] 
 
>                  Re[2.71828
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>                                                                            ]
 
>                 )) / 
 
                             (0. - 1.5 I) (1.5708 + ww) 2
>             (-1. Im[2.71828                          ]  + 
 
>               Im[2.71828
 
                   (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww]) 2
>                                                                          ] \
 
                                (0. - 1.5 I) (1.5708 + ww) 2
>                - 1. Re[2.71828                          ]  + 
 
>               Re[2.71828
 
                   (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww]) 2
>                                                                          ] )
 
>            )) / 
 
                        (0. - 1.5 I) (1.5708 + ww)
>        (-1. Re[2.71828                          ] + 
 
                     (0. - 2.35619 I) + (0. + 1.5 I) (1.5708 + Conjugate[ww])
>          Re[2.71828                                                        ]
 
>          )}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[225]= {{x -> Indeterminate}}

Out[224]= 4

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                            0 ComplexInfinity
Infinity::indet: Indeterminate expression -(-----------------) encountered.
                                                 Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[223]= {{x -> Indeterminate}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

                                 1
Power::infy: Infinite expression - encountered.
                                 0

                                          0 ComplexInfinity
Infinity::indet: Indeterminate expression ----------------- encountered.
                                               Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

Infinity::indet: Indeterminate expression 0 ComplexInfinity encountered.

General::stop: Further output of Infinity::indet
     will be suppressed during this calculation.

Out[222]= {{x -> Indeterminate}}

Out[221]= {{x -> -0.375 - 0.0132793 I}}

                                 1
Power::infy: Infinite expression - encountered.
                                 0

N::meprec: Internal precision limit $MaxExtraPrecision = 50.
                                             Pi
                                         Cos[--]
                                   Pi        20        (4 I)/5 Pi
     reached while evaluating -Cos[--] + ------- - Im[E          ] + 
                                   20    Sqrt[2]
                                             Pi
         (4 I)/5 Pi        (4 I)/5 Pi    Sin[--]
     Im[E          ]   Re[E          ]       20
     --------------- - --------------- - -------.
         Sqrt[2]           Sqrt[2]       Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

N::meprec: Internal precision limit $MaxExtraPrecision = 50.
                                             Pi
                                         Cos[--]
                                   Pi        20        (4 I)/5 Pi
     reached while evaluating -Cos[--] + ------- - Im[E          ] + 
                                   20    Sqrt[2]
                                             Pi
         (4 I)/5 Pi        (4 I)/5 Pi    Sin[--]
     Im[E          ]   Re[E          ]       20
     --------------- - --------------- - -------.
         Sqrt[2]           Sqrt[2]       Sqrt[2]

                                 1
Power::infy: Infinite expression - encountered.
                                 0

General::stop: Further output of Power::infy
     will be suppressed during this calculation.

N::meprec: Internal precision limit $MaxExtraPrecision = 50.
                                             Pi
                                         Cos[--]
                                   Pi        20        (4 I)/5 Pi
     reached while evaluating -Cos[--] + ------- - Im[E          ] + 
                                   20    Sqrt[2]
                                             Pi
         (4 I)/5 Pi        (4 I)/5 Pi    Sin[--]
     Im[E          ]   Re[E          ]       20
     --------------- - --------------- - -------.
         Sqrt[2]           Sqrt[2]       Sqrt[2]

General::stop: Further output of N::meprec
     will be suppressed during this calculation.

Infinity::indet: 
                                 Pi
   Indeterminate expression -Cos[--] + ComplexInfinity + ComplexInfinity - 
                                 20
         (4 I)/5 Pi
     Im[E          ] encountered.

Out[220]= {{x -> Indeterminate}}

Out[218]= {{x -> 0. + 0.762757 I}}

[Calculating...]

Out[207]= 0.413304

                                                              1 + Sqrt[5]
                                                       ArcCos[-----------]
                               2                                   4
Out[206]= Sqrt[----------------------------------] Sin[-------------------]
                                               2                2
                                  (1 + Sqrt[5])
               (1 + Sqrt[5]) (1 - --------------)
                                        16

N[w[2,3,4,3]]/N[w[4,3,4,3]]

Out[216]= {{x -> 0. - 1.63928 I}}

Out[215]= {{x -> 0. - 1.63928 I}}

Out[214]= 0.485868

Out[213]= 0.393076

Out[212]= 1.05146

                       1 + Sqrt[5]
                ArcCos[-----------]
                            4
          2 Sin[-------------------]
                         2
Out[211]= --------------------------
                                 2
                    (1 + Sqrt[5])
           Sqrt[1 - --------------]
                          16

Out[210]= 0.413304

Out[209]= 2.41953

Out[208]= 2.41953

Out[205]= {{x -> 0. - 1.63928 I}}

Out[204]= {{x -> 
 
            Pi
        Cos[--]       (-3 I)/10 Pi                            (-3 I)/10 Pi
            20    Im[E            ]       (-3 I)/10 Pi    Re[E            ]
>      (------- + ----------------- - Re[E            ] + ----------------- + 
        Sqrt[2]        Sqrt[2]                                 Sqrt[2]
 
                         Pi
                     Sin[--]
               Pi        20
>          Sin[--] - -------) / 
               20    Sqrt[2]
 
              Pi 2       (-3 I)/10 Pi 2       (-3 I)/10 Pi 2       Pi 2
>        (Cos[--]  - Im[E            ]  - Re[E            ]  + Sin[--] ) + 
              20                                                   20
 
                                    Pi
                                Cos[--]       (-3 I)/10 Pi
               1           Pi       20    Im[E            ]
>       (I (------- - (Cos[--] (------- + ----------------- - 
            Sqrt[2]        20   Sqrt[2]        Sqrt[2]
 
                                                                        Pi
                                          (-3 I)/10 Pi              Sin[--]
                      (-3 I)/10 Pi    Re[E            ]       Pi        20
>                 Re[E            ] + ----------------- + Sin[--] - -------))\
                                           Sqrt[2]            20    Sqrt[2]
 
                      Pi 2       (-3 I)/10 Pi 2       (-3 I)/10 Pi 2
>              / (Cos[--]  - Im[E            ]  - Re[E            ]  + 
                      20
 
                    Pi 2         (-3 I)/10 Pi
>               Sin[--] ) + (Im[E            ] 
                    20
 
                     Pi
                 Cos[--]       (-3 I)/10 Pi
                     20    Im[E            ]       (-3 I)/10 Pi
>               (------- + ----------------- - Re[E            ] + 
                 Sqrt[2]        Sqrt[2]
 
                                                    Pi
                      (-3 I)/10 Pi              Sin[--]
                  Re[E            ]       Pi        20
>                 ----------------- + Sin[--] - -------)) / 
                       Sqrt[2]            20    Sqrt[2]
 
                   Pi 2       (-3 I)/10 Pi 2       (-3 I)/10 Pi 2       Pi 2
>             (Cos[--]  - Im[E            ]  - Re[E            ]  + Sin[--] ))
                   20                                                   20
 
                     (-3 I)/10 Pi        Pi
>          ) / (-Re[E            ] + Sin[--])}}
                                         20

Out[203]= {{x -> 
 
              (-3 I)/5 ww      (-3 I)/5 (Pi/2 + ww)
>      (-(Im[E           ] Im[E                    ]) - 
 
               (-3 I)/5 (Pi/2 + ww)
>          Im[E                    ] 
 
                (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>           Im[E                                   ] + 
 
               (-3 I)/5 ww
>          Im[E           ] Im[E
 
              (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                                                         ] + 
 
               (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>          Im[E                                   ] 
 
                (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>           Im[E                                            ] - 
 
               (-3 I)/5 ww      (-3 I)/5 (Pi/2 + ww)
>          Re[E           ] Re[E                    ] - 
 
               (-3 I)/5 (Pi/2 + ww)
>          Re[E                    ] 
 
                (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>           Re[E                                   ] + 
 
               (-3 I)/5 ww
>          Re[E           ] Re[E
 
              (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                                                         ] + 
 
               (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>          Re[E                                   ] 
 
                (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>           Re[E                                            ]) / 
 
               (-3 I)/5 (Pi/2 + ww) 2
>        (-Im[E                    ]  + 
 
               (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww]) 2
>          Im[E                                            ]  - 
 
               (-3 I)/5 (Pi/2 + ww) 2
>          Re[E                    ]  + 
 
               (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww]) 2
>          Re[E                                            ] ) + 
 
                 (-3 I)/5 ww        (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>       (I (-Im[E           ] - Im[E                                   ] + 
 
                  (-3 I)/5 (Pi/2 + ww)
>            (Im[E                    ] 
 
                       (-3 I)/5 ww      (-3 I)/5 (Pi/2 + ww)
>               (-(Im[E           ] Im[E                    ]) - 
 
                      (-3 I)/5 (Pi/2 + ww)
>                 Im[E                    ] 
 
                       (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                  Im[E                                   ] + 
 
                      (-3 I)/5 ww
>                 Im[E           ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] + 
 
                      (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                 Im[E                                   ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] - 
 
                      (-3 I)/5 ww      (-3 I)/5 (Pi/2 + ww)
>                 Re[E           ] Re[E                    ] - 
 
                      (-3 I)/5 (Pi/2 + ww)
>                 Re[E                    ] 
 
                       (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                  Re[E                                   ] + 
 
                      (-3 I)/5 ww
>                 Re[E           ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ] + 
 
                      (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                 Re[E                                   ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ])) / 
 
                    (-3 I)/5 (Pi/2 + ww) 2
>             (-Im[E                    ]  + 
 
                    (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww]) 2
>               Im[E                                            ]  - 
 
                    (-3 I)/5 (Pi/2 + ww) 2
>               Re[E                    ]  + 
 
                    (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww]) 2
>               Re[E                                            ] ) + 
 
                  (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>            (Im[E                                            ] 
 
                       (-3 I)/5 ww      (-3 I)/5 (Pi/2 + ww)
>               (-(Im[E           ] Im[E                    ]) - 
 
                      (-3 I)/5 (Pi/2 + ww)
>                 Im[E                    ] 
 
                       (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                  Im[E                                   ] + 
 
                      (-3 I)/5 ww
>                 Im[E           ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] + 
 
                      (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                 Im[E                                   ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Im[E                                            ] - 
 
                      (-3 I)/5 ww      (-3 I)/5 (Pi/2 + ww)
>                 Re[E           ] Re[E                    ] - 
 
                      (-3 I)/5 (Pi/2 + ww)
>                 Re[E                    ] 
 
                       (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                  Re[E                                   ] + 
 
                      (-3 I)/5 ww
>                 Re[E           ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ] + 
 
                      (-3 I)/4 Pi + (3 I)/5 Conjugate[ww]
>                 Re[E                                   ] 
 
                       (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>                  Re[E                                            ])) / 
 
                    (-3 I)/5 (Pi/2 + ww) 2
>             (-Im[E                    ]  + 
 
                    (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww]) 2
>               Im[E                                            ]  - 
 
                    (-3 I)/5 (Pi/2 + ww) 2
>               Re[E                    ]  + 
 
                    (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww]) 2
>               Re[E                                            ] ))) / 
 
               (-3 I)/5 (Pi/2 + ww)
>        (-Re[E                    ] + 
 
               (-3 I)/4 Pi + (3 I)/5 (Pi/2 + Conjugate[ww])
>          Re[E                                            ])}}

                2                     2
               u  Conjugate[x]       u  Conjugate[y]
          (x + ---------------) (y + ---------------)
                         2                     2
                   Abs[u]                Abs[u]
Out[202]= -------------------------------------------
                               4

                 2
                u  Conjugate[x y]
          x y + -----------------
                           2
                     Abs[u]
Out[201]= -----------------------
                     2

