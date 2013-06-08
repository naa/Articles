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