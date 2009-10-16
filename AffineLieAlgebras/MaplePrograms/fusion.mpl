Psi1:=(p,c,d)->
add(
    add(
        add(
            (-1)^(k-d-c+p-4*(lk+mk)+7)*
            binomial(p,k-1)*binomial(k-1,lk-1)*
            binomial(p-k+1,mk-1)*
            binomial(p-k+1,-d+2*k-5*(lk-1)-2)*
            binomial(k-1,p-c-2*k+2-5*(mk-1)),
            mk=1..p-k+2),
        lk=1..k),
    k=1..p+1);

Gamma:=(p,a,b)->
add(
    add(
        add(
            (-1)^(k-2*(lk+mk)+a+b+4)*
            binomial(p-1,k-1)*binomial(k-1,lk-1)*binomial(p-k,mk-1)*
            binomial(p-k,b+k-3*lk+2)*binomial(k-1,a-k-3*mk+4),
            mk=1..p-k+1),
        lk=1..k),
    k=1..p);

genf:=proc(func,p,fr,t)
    local i,j,aa;
    aa:=[];
    for i from fr to t do
      for j from fr to t do
      	  aa:=[op(aa),[i,j,func(p,i,j)]];
      od;
    od;
end;
