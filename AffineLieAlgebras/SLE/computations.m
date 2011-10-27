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

