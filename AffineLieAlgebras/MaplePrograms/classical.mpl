if not assigned(maplev_print) then  maplev_print := print fi:
interface(prettyprint=1,verboseproc=2,errorbreak=0,
screenheight=9999,warnlevel=2,errorcursor=false):
kernelopts(printbytes=false):

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSignR;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->weyl['weight_coords'](subs(eps=0,x)),thestar);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 1 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:
Warning, `R` is implicitly declared local to procedure `multiplicities`
Warning, `starSign` is implicitly declared local to procedure `multiplicities`

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->weyl['weight_coords'](subs(eps=0,x)),thestar);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 1 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> weights(B2);
                                  weights(B2)

> weyl['weights'](B2);
                                 e1     e2
                               [---- + ----, e2]
                                 2      2

> weyl['weight_mults'](e1/2+e2/2,40*(e1/2+e2/2),B2);
                                       0

> weyl['weight_mults'](e1/2+e2/2,2*(e1/2+e2/2),B2);
                                       0

> weyl['weight_mults'](e1/2+e2/2,8*(e1/2+e2/2),B2);
                                       0

> with(LinearAlgebra):
> M:=Matrix([[-1, 2, 0,-1,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0],
[ 0,-1, 1, 1, 0,-1,-1, 0, 1, 0, 1,-1, 0, 0, 0],
[ 0, 0,-1, 1, 1, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0],
[ 0, 0, 0,-1, 1, 0, 1, 0,-1, 0,-1, 0, 1, 0, 0],
[ 0, 0, 0, 0,-1, 1, 1, 0, 0,-1,-1, 0, 0, 1, 0],
[ 0, 0, 0, 0, 0,-1, 0, 1, 1, 0, 0,-1, 0,-1, 1],
[ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 1, 0,-1, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1, 0, 0,-1, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1, 0, 0,-1],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 1, 1, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1]]);

M:=Matrix([[-1, 2, 0,-1,-1, 1, 1,-1, 0, 0, 0, 0, 0, 0, 0],
[ 0,-1, 1, 1, 0,-1,-1, 0, 1, 0, 1,-1, 0, 0, 0],
[ 0, 0,-1, 1, 1, 0,-1, 0,-1, 1, 0, 1,-1, 0, 0],
[ 0, 0, 0,-1, 1, 0, 1, 0,-1, 0,-1, 0, 1, 0, 0],
[ 0, 0, 0, 0,-1, 1, 1, 0, 0,-1,-1, 0, 0, 1, 0],
[ 0, 0, 0, 0, 0,-1, 0, 1, 1, 0, 0,-1, 0,-1, 1],
[ 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 1, 0,-1, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1, 0, 0,-1, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 1, 0, 0,-1],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0, 1, 1, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1, 0],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 1],
[ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1]]);
          [-1 , 2 , 0 , -1 , -1 , 1 , 1 , -1 , 0 , 0 , 0 , 0 , 0 , 0 , 0]
          [                                                             ]
          [0 , -1 , 1 , 1 , 0 , -1 , -1 , 0 , 1 , 0 , 1 , -1 , 0 , 0 , 0]
          [                                                             ]
          [0 , 0 , -1 , 1 , 1 , 0 , -1 , 0 , -1 , 1 , 0 , 1 , -1 , 0 , 0]
          [                                                             ]
          [0 ,  0 , 0 , -1 , 1 , 0 , 1 , 0 , -1 , 0 , -1 , 0 , 1 , 0 , 0]
          [                                                             ]
          [0 ,  0 , 0 , 0 , -1 , 1 , 1 , 0 , 0 , -1 , -1 , 0 , 0 , 1 , 0]
          [                                                             ]
          [0 ,  0 , 0 , 0 , 0 , -1 , 0 , 1 , 1 , 0 , 0 , -1 , 0 , -1 , 1]
          [                                                             ]
          [0 ,  0 ,  0 , 0 , 0 , 0 , -1 , 1 , 0 , 0 , 1 , 0 , -1 , 0 , 0]
          [                                                             ]
     M := [0 ,  0 ,  0 , 0 , 0 , 0 , 0 , -1 , 1 , 0 , 1 , 0 , 0 , -1 , 0]
          [                                                             ]
          [0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , -1 , 1 , 0 , 1 , 0 , 0 , -1]
          [                                                             ]
          [0 ,  0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , -1 , 0 , 0 , 1 , 1 , 0]
          [                                                             ]
          [0 ,  0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , 1 , 0 , 0 , 0]
          [                                                             ]
          [0 ,  0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , 1 , 0 , 0]
          [                                                             ]
          [0 ,  0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , 1 , 0]
          [                                                             ]
          [0 ,  0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , 1]
          [                                                             ]
          [0 ,  0 ,  0 ,  0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1]


> MatrixInverse(M);
    [-1 , -2 , -2 , -3 , -4 , -3 , -4 , -6 , -6 , -4 , -5 , -8 , -9 , -8 , -5]

    [0 , -1 , -1 , -2 , -3 , -2 , -3 , -5 , -5 , -3 , -4 , -7 , -8 , -7 , -4]

    [0 , 0 , -1 , -1 , -2 , -2 , -2 , -4 , -4 , -3 , -3 , -6 , -7 , -6 , -4]

    [0 , 0 , 0 , -1 , -1 , -1 , -2 , -3 , -3 , -2 , -3 , -5 , -6 , -5 , -3]

    [0 , 0 , 0 , 0 , -1 , -1 , -1 , -2 , -3 , -2 , -2 , -4 , -5 , -5 , -3]

    [0 , 0 , 0 , 0 , 0 , -1 , 0 , -1 , -2 , -2 , -1 , -2 , -4 , -4 , -3]

    [0 , 0 , 0 , 0 , 0 , 0 , -1 , -1 , -1 , -1 , -2 , -3 , -3 , -3 , -2]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , -1 , -1 , -1 , -2 , -3 , -3 , -2]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , -1 , 0 , -1 , -2 , -3 , -2]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , 0 , 0 , -1 , -2 , -2]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , -1 , -1 , -1 , -1]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , -1 , -1 , -1]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , -1 , -1]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1 , -1]

    [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , -1]

> weyl['weight_mults'](e1/2+e2/2,(e1/2+e2/2),B2);
                                       1

> weyl['weight_mults'](e1/2+e2/2,41*(e1/2+e2/2),B2);
                                       0

> weyl['weight_mults'](e1/2+e2/2,3*(e1/2+e2/2),B2);
                                       0

> weyl['weight_mults']([1,0],[10,0],B2);
Error, (in coxeter/linsolve) invalid input: SolveTools:-Linear expects its 1st
argument, eqs, to be of type {set, list}({algebraic, algebraic = algebraic}),
but received {[9, 0], c[2], c[1]-c[2]}
> weyl['weight_mults'](e1/2+e2/2,4*(e1/2+e2/2),B2);
                                       0

> weyl['weight_mults'](e1/2+e2/2,5*(e1/2+e2/2),B2);
                                       0

> weyl['weights'](B2);
                                 e1     e2
                               [---- + ----, e2]
                                 2      2

> wg:=weyl['weights'](B2);
                                    e1     e2
                            wg := [---- + ----, e2]
                                    2      2

> weyl['weight_mults'](wg[1],wg[1]*5,B2);
                                       0

> weyl['weight_mults'](wg[1]*5,wg[1],B2);
                                       3

> weyl['weight_mults'](wg[1]*555,wg[1],B2);
Interrupted
> weyl['weight_mults'](wg[1]*55,wg[1],B2);
                                      28

> multiplicities(wg[1]*55,wg[1],B2);
Error, (in finite_orbit) cannot determine if this expression is true or false:
0 < iprod(2*e1,1/2*e1+3/2*e2)
> iprod(e1,e1);
                                 iprod(e1, e1)

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2
iprod := (v1,v2) -> coeff(v1,lambda0)*c\

oeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda
oeff(v2,delta)+coeff(v1,delta)*coeff(v2\

,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,d
,lambda0)+coxeter['iprod'](subs({lambda\

0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
iprod := (v1, v2) -> coeff(v1, lambda0) coeff(v2, delta)

     + coeff(v1, delta) coeff(v2, lambda0) + coxeter['iprod'](

    subs({eps = 0, lambda0 = 0, delta = 0}, v1),

    subs({eps = 0, delta = 0, lamda0 = 0}, v2))


weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->weyl['weight_coords'](subs(eps=0,x)),thestar);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 1 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> multiplicities(wg[1]*55,wg[1],B2);
Error, (in unknown) invalid input: weyl/weight_coords uses a 2nd argument, R,
which is missing
> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2
iprod := (v1,v2) -> coeff(v1,lambda0)*c\

oeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda
oeff(v2,delta)+coeff(v1,delta)*coeff(v2\

,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,d
,lambda0)+coxeter['iprod'](subs({lambda\

0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
iprod := (v1, v2) -> coeff(v1, lambda0) coeff(v2, delta)

     + coeff(v1, delta) coeff(v2, lambda0) + coxeter['iprod'](

    subs({eps = 0, lambda0 = 0, delta = 0}, v1),

    subs({eps = 0, delta = 0, lamda0 = 0}, v2))


weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->weyl['weight_coords'](subs(eps=0,x),R),thestar);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 1 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> multiplicities(wg[1]*55,wg[1],B2);
                                       0

> multiplicities(wg[1],wg[1],B2);
                                       1

> multiplicities(wg[1]*3,wg[1],B2);
                                       0

> star(B2);
[eps, e1 - eps, -e1 + e2 - eps, -e1 + 2 e2 + eps, 2 e1 + e2 + eps,

    2 e1 + 2 e2 - eps, 3 e2 - eps, e1 + 3 e2 + eps]

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2
iprod := (v1,v2) -> coeff(v1,lambda0)*c\

oeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda
oeff(v2,delta)+coeff(v1,delta)*coeff(v2\

,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,d
,lambda0)+coxeter['iprod'](subs({lambda\

0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
iprod := (v1, v2) -> coeff(v1, lambda0) coeff(v2, delta)

     + coeff(v1, delta) coeff(v2, lambda0) + coxeter['iprod'](

    subs({eps = 0, lambda0 = 0, delta = 0}, v1),

    subs({eps = 0, delta = 0, lamda0 = 0}, v2))


weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->weyl['weight_coords'](subs(eps=0,x),R),thestar);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 2 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> multiplicities(wg[1]*3,wg[1],B2);
                                       0

> star(B2);
[eps, e1 - eps, -e1 + e2 - eps, -e1 + 2 e2 + eps, 2 e1 + e2 + eps,

    2 e1 + 2 e2 - eps, 3 e2 - eps, e1 + 3 e2 + eps]

>     starSign:=map(x->coeff(x,eps),thestar);
                                 starSign := 0

> thestar:=star(B2);
thestar := [eps, e1 - eps, -e1 + e2 - eps, -e1 + 2 e2 + eps, 2 e1 + e2 + eps,

    2 e1 + 2 e2 - eps, 3 e2 - eps, e1 + 3 e2 + eps]

>     starSign:=map(x->coeff(x,eps),thestar);
                   starSign := [1, -1, -1, 1, 1, -1, -1, 1]

>     starW:=map(x->weyl['weight_coords'](subs(eps=0,x),B2),thestar);
 starW := [[0, 0], [2, -1], [-2, 2], [-2, 3], [4, -1], [4, 0], [0, 3], [2, 2]]

> nops(starW);
                                       8

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2
iprod := (v1,v2) -> coeff(v1,lambda0)*c\

oeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda
oeff(v2,delta)+coeff(v1,delta)*coeff(v2\

,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,d
,lambda0)+coxeter['iprod'](subs({lambda\

0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
iprod := (v1, v2) -> coeff(v1, lambda0) coeff(v2, delta)

     + coeff(v1, delta) coeff(v2, lambda0) + coxeter['iprod'](

    subs({eps = 0, lambda0 = 0, delta = 0}, v1),

    subs({eps = 0, delta = 0, lamda0 = 0}, v2))


weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->weyl['weight_coords'](subs(eps=0,x),R),thestar);

    for i from 2 to n do
        v:=wts[i];
        print(v);
        m:=0;
        for k from 2 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> multiplicities(wg[1]*3,wg[1],B2);
                                   e1    3 e2
                                  ---- + ----
                                   2      2

                                   e1     e2
                                  ---- + ----
                                   2      2

                                       0

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2
iprod := (v1,v2) -> coeff(v1,lambda0)*c\

oeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda
oeff(v2,delta)+coeff(v1,delta)*coeff(v2\

,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,d
,lambda0)+coxeter['iprod'](subs({lambda\

0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
iprod := (v1, v2) -> coeff(v1, lambda0) coeff(v2, delta)

     + coeff(v1, delta) coeff(v2, lambda0) + coxeter['iprod'](

    subs({eps = 0, lambda0 = 0, delta = 0}, v1),

    subs({eps = 0, delta = 0, lamda0 = 0}, v2))


weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->subs(eps=0,x),thestar);

    for i from 2 to n do
        v:=wts[i];
        print(v);
        m:=0;
        for k from 2 to nops(starW) do
            u:=v+starW[k];
            print(u);
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> multiplicities(wg[1]*3,wg[1],B2);
                                   e1    3 e2
                                  ---- + ----
                                   2      2

                                  3 e1   3 e2
                                  ---- + ----
                                   2      2

                                    e1    5 e2
                                 - ---- + ----
                                    2      2

                                    e1    7 e2
                                 - ---- + ----
                                    2      2

                                  5 e1   5 e2
                                  ---- + ----
                                   2      2

                                  5 e1   7 e2
                                  ---- + ----
                                   2      2

                                   e1    9 e2
                                  ---- + ----
                                   2      2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                   e1     e2
                                  ---- + ----
                                   2      2

                                  3 e1    e2
                                  ---- + ----
                                   2      2

                                    e1    3 e2
                                 - ---- + ----
                                    2      2

                                    e1    5 e2
                                 - ---- + ----
                                    2      2

                                  5 e1   3 e2
                                  ---- + ----
                                   2      2

                                  5 e1   5 e2
                                  ---- + ----
                                   2      2

                                   e1    7 e2
                                  ---- + ----
                                   2      2

                                  3 e1   7 e2
                                  ---- + ----
                                   2      2

                                       2

> weyl['weight_mults'](wg[1]*3,wg[1],B2);
                                       2

> weyl['weight_mults'](wg[1]*55,wg[1],B2);
                                      28

> multiplicities(wg[1]*55,wg[1],B2);
                                 53 e1   55 e2
                                 ----- + -----
                                   2       2

                                 55 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   57 e2
                                 ----- + -----
                                   2       2

                                 51 e1   59 e2
                                 ----- + -----
                                   2       2

                                 57 e1   57 e2
                                 ----- + -----
                                   2       2

                                 57 e1   59 e2
                                 ----- + -----
                                   2       2

                                 53 e1   61 e2
                                 ----- + -----
                                   2       2

                                 55 e1   61 e2
                                 ----- + -----
                                   2       2

                                 51 e1   55 e2
                                 ----- + -----
                                   2       2

                                 53 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   57 e2
                                 ----- + -----
                                   2       2

                                 49 e1   59 e2
                                 ----- + -----
                                   2       2

                                 55 e1   57 e2
                                 ----- + -----
                                   2       2

                                 55 e1   59 e2
                                 ----- + -----
                                   2       2

                                 51 e1   61 e2
                                 ----- + -----
                                   2       2

                                 53 e1   61 e2
                                 ----- + -----
                                   2       2

                                 53 e1   53 e2
                                 ----- + -----
                                   2       2

                                 55 e1   53 e2
                                 ----- + -----
                                   2       2

                                 51 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   57 e2
                                 ----- + -----
                                   2       2

                                 57 e1   55 e2
                                 ----- + -----
                                   2       2

                                 57 e1   57 e2
                                 ----- + -----
                                   2       2

                                 53 e1   59 e2
                                 ----- + -----
                                   2       2

                                 55 e1   59 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   57 e2
                                 ----- + -----
                                   2       2

                                 47 e1   59 e2
                                 ----- + -----
                                   2       2

                                 53 e1   57 e2
                                 ----- + -----
                                   2       2

                                 53 e1   59 e2
                                 ----- + -----
                                   2       2

                                 49 e1   61 e2
                                 ----- + -----
                                   2       2

                                 51 e1   61 e2
                                 ----- + -----
                                   2       2

                                 51 e1   53 e2
                                 ----- + -----
                                   2       2

                                 53 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   57 e2
                                 ----- + -----
                                   2       2

                                 55 e1   55 e2
                                 ----- + -----
                                   2       2

                                 55 e1   57 e2
                                 ----- + -----
                                   2       2

                                 51 e1   59 e2
                                 ----- + -----
                                   2       2

                                 53 e1   59 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   57 e2
                                 ----- + -----
                                   2       2

                                 45 e1   59 e2
                                 ----- + -----
                                   2       2

                                 51 e1   57 e2
                                 ----- + -----
                                   2       2

                                 51 e1   59 e2
                                 ----- + -----
                                   2       2

                                 47 e1   61 e2
                                 ----- + -----
                                   2       2

                                 49 e1   61 e2
                                 ----- + -----
                                   2       2

                                 49 e1   53 e2
                                 ----- + -----
                                   2       2

                                 51 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   57 e2
                                 ----- + -----
                                   2       2

                                 53 e1   55 e2
                                 ----- + -----
                                   2       2

                                 53 e1   57 e2
                                 ----- + -----
                                   2       2

                                 49 e1   59 e2
                                 ----- + -----
                                   2       2

                                 51 e1   59 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   57 e2
                                 ----- + -----
                                   2       2

                                 43 e1   59 e2
                                 ----- + -----
                                   2       2

                                 49 e1   57 e2
                                 ----- + -----
                                   2       2

                                 49 e1   59 e2
                                 ----- + -----
                                   2       2

                                 45 e1   61 e2
                                 ----- + -----
                                   2       2

                                 47 e1   61 e2
                                 ----- + -----
                                   2       2

                                 51 e1   51 e2
                                 ----- + -----
                                   2       2

                                 53 e1   51 e2
                                 ----- + -----
                                   2       2

                                 49 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 55 e1   53 e2
                                 ----- + -----
                                   2       2

                                 55 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   57 e2
                                 ----- + -----
                                   2       2

                                 53 e1   57 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   57 e2
                                 ----- + -----
                                   2       2

                                 51 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   57 e2
                                 ----- + -----
                                   2       2

                                 47 e1   59 e2
                                 ----- + -----
                                   2       2

                                 49 e1   59 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   57 e2
                                 ----- + -----
                                   2       2

                                 41 e1   59 e2
                                 ----- + -----
                                   2       2

                                 47 e1   57 e2
                                 ----- + -----
                                   2       2

                                 47 e1   59 e2
                                 ----- + -----
                                   2       2

                                 43 e1   61 e2
                                 ----- + -----
                                   2       2

                                 45 e1   61 e2
                                 ----- + -----
                                   2       2

                                 49 e1   51 e2
                                 ----- + -----
                                   2       2

                                 51 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 53 e1   53 e2
                                 ----- + -----
                                   2       2

                                 53 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   57 e2
                                 ----- + -----
                                   2       2

                                 51 e1   57 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   57 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   57 e2
                                 ----- + -----
                                   2       2

                                 45 e1   59 e2
                                 ----- + -----
                                   2       2

                                 47 e1   59 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   57 e2
                                 ----- + -----
                                   2       2

                                 39 e1   59 e2
                                 ----- + -----
                                   2       2

                                 45 e1   57 e2
                                 ----- + -----
                                   2       2

                                 45 e1   59 e2
                                 ----- + -----
                                   2       2

                                 41 e1   61 e2
                                 ----- + -----
                                   2       2

                                 43 e1   61 e2
                                 ----- + -----
                                   2       2

                                 47 e1   51 e2
                                 ----- + -----
                                   2       2

                                 49 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   53 e2
                                 ----- + -----
                                   2       2

                                 51 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   57 e2
                                 ----- + -----
                                   2       2

                                 49 e1   57 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   57 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   57 e2
                                 ----- + -----
                                   2       2

                                 43 e1   59 e2
                                 ----- + -----
                                   2       2

                                 45 e1   59 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   57 e2
                                 ----- + -----
                                   2       2

                                 37 e1   59 e2
                                 ----- + -----
                                   2       2

                                 43 e1   57 e2
                                 ----- + -----
                                   2       2

                                 43 e1   59 e2
                                 ----- + -----
                                   2       2

                                 39 e1   61 e2
                                 ----- + -----
                                   2       2

                                 41 e1   61 e2
                                 ----- + -----
                                   2       2

                                 49 e1   49 e2
                                 ----- + -----
                                   2       2

                                 51 e1   49 e2
                                 ----- + -----
                                   2       2

                                 47 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 53 e1   51 e2
                                 ----- + -----
                                   2       2

                                 53 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 51 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   57 e2
                                 ----- + -----
                                   2       2

                                 47 e1   57 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   57 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   57 e2
                                 ----- + -----
                                   2       2

                                 41 e1   59 e2
                                 ----- + -----
                                   2       2

                                 43 e1   59 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   57 e2
                                 ----- + -----
                                   2       2

                                 35 e1   59 e2
                                 ----- + -----
                                   2       2

                                 41 e1   57 e2
                                 ----- + -----
                                   2       2

                                 41 e1   59 e2
                                 ----- + -----
                                   2       2

                                 37 e1   61 e2
                                 ----- + -----
                                   2       2

                                 39 e1   61 e2
                                 ----- + -----
                                   2       2

                                 47 e1   49 e2
                                 ----- + -----
                                   2       2

                                 49 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 51 e1   51 e2
                                 ----- + -----
                                   2       2

                                 51 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 49 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   57 e2
                                 ----- + -----
                                   2       2

                                 45 e1   57 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   57 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   57 e2
                                 ----- + -----
                                   2       2

                                 39 e1   59 e2
                                 ----- + -----
                                   2       2

                                 41 e1   59 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   57 e2
                                 ----- + -----
                                   2       2

                                 33 e1   59 e2
                                 ----- + -----
                                   2       2

                                 39 e1   57 e2
                                 ----- + -----
                                   2       2

                                 39 e1   59 e2
                                 ----- + -----
                                   2       2

                                 35 e1   61 e2
                                 ----- + -----
                                   2       2

                                 37 e1   61 e2
                                 ----- + -----
                                   2       2

                                 45 e1   49 e2
                                 ----- + -----
                                   2       2

                                 47 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   51 e2
                                 ----- + -----
                                   2       2

                                 49 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 47 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   57 e2
                                 ----- + -----
                                   2       2

                                 43 e1   57 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   57 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   57 e2
                                 ----- + -----
                                   2       2

                                 37 e1   59 e2
                                 ----- + -----
                                   2       2

                                 39 e1   59 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   57 e2
                                 ----- + -----
                                   2       2

                                 31 e1   59 e2
                                 ----- + -----
                                   2       2

                                 37 e1   57 e2
                                 ----- + -----
                                   2       2

                                 37 e1   59 e2
                                 ----- + -----
                                   2       2

                                 33 e1   61 e2
                                 ----- + -----
                                   2       2

                                 35 e1   61 e2
                                 ----- + -----
                                   2       2

                                 47 e1   47 e2
                                 ----- + -----
                                   2       2

                                 49 e1   47 e2
                                 ----- + -----
                                   2       2

                                 45 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 51 e1   49 e2
                                 ----- + -----
                                   2       2

                                 51 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 49 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 45 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   57 e2
                                 ----- + -----
                                   2       2

                                 41 e1   57 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   57 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   57 e2
                                 ----- + -----
                                   2       2

                                 35 e1   59 e2
                                 ----- + -----
                                   2       2

                                 37 e1   59 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   57 e2
                                 ----- + -----
                                   2       2

                                 29 e1   59 e2
                                 ----- + -----
                                   2       2

                                 35 e1   57 e2
                                 ----- + -----
                                   2       2

                                 35 e1   59 e2
                                 ----- + -----
                                   2       2

                                 31 e1   61 e2
                                 ----- + -----
                                   2       2

                                 33 e1   61 e2
                                 ----- + -----
                                   2       2

                                 45 e1   47 e2
                                 ----- + -----
                                   2       2

                                 47 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 49 e1   49 e2
                                 ----- + -----
                                   2       2

                                 49 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 47 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 43 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   57 e2
                                 ----- + -----
                                   2       2

                                 39 e1   57 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   57 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   57 e2
                                 ----- + -----
                                   2       2

                                 33 e1   59 e2
                                 ----- + -----
                                   2       2

                                 35 e1   59 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   57 e2
                                 ----- + -----
                                   2       2

                                 27 e1   59 e2
                                 ----- + -----
                                   2       2

                                 33 e1   57 e2
                                 ----- + -----
                                   2       2

                                 33 e1   59 e2
                                 ----- + -----
                                   2       2

                                 29 e1   61 e2
                                 ----- + -----
                                   2       2

                                 31 e1   61 e2
                                 ----- + -----
                                   2       2

                                 43 e1   47 e2
                                 ----- + -----
                                   2       2

                                 45 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   49 e2
                                 ----- + -----
                                   2       2

                                 47 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 45 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 41 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   57 e2
                                 ----- + -----
                                   2       2

                                 37 e1   57 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   57 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   57 e2
                                 ----- + -----
                                   2       2

                                 31 e1   59 e2
                                 ----- + -----
                                   2       2

                                 33 e1   59 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   57 e2
                                 ----- + -----
                                   2       2

                                 25 e1   59 e2
                                 ----- + -----
                                   2       2

                                 31 e1   57 e2
                                 ----- + -----
                                   2       2

                                 31 e1   59 e2
                                 ----- + -----
                                   2       2

                                 27 e1   61 e2
                                 ----- + -----
                                   2       2

                                 29 e1   61 e2
                                 ----- + -----
                                   2       2

                                 45 e1   45 e2
                                 ----- + -----
                                   2       2

                                 47 e1   45 e2
                                 ----- + -----
                                   2       2

                                 43 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 49 e1   47 e2
                                 ----- + -----
                                   2       2

                                 49 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 47 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 43 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 39 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   57 e2
                                 ----- + -----
                                   2       2

                                 35 e1   57 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   57 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   57 e2
                                 ----- + -----
                                   2       2

                                 29 e1   59 e2
                                 ----- + -----
                                   2       2

                                 31 e1   59 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   57 e2
                                 ----- + -----
                                   2       2

                                 23 e1   59 e2
                                 ----- + -----
                                   2       2

                                 29 e1   57 e2
                                 ----- + -----
                                   2       2

                                 29 e1   59 e2
                                 ----- + -----
                                   2       2

                                 25 e1   61 e2
                                 ----- + -----
                                   2       2

                                 27 e1   61 e2
                                 ----- + -----
                                   2       2

                                 43 e1   45 e2
                                 ----- + -----
                                   2       2

                                 45 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 47 e1   47 e2
                                 ----- + -----
                                   2       2

                                 47 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 45 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 41 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 37 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   57 e2
                                 ----- + -----
                                   2       2

                                 33 e1   57 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   57 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   57 e2
                                 ----- + -----
                                   2       2

                                 27 e1   59 e2
                                 ----- + -----
                                   2       2

                                 29 e1   59 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   57 e2
                                 ----- + -----
                                   2       2

                                 21 e1   59 e2
                                 ----- + -----
                                   2       2

                                 27 e1   57 e2
                                 ----- + -----
                                   2       2

                                 27 e1   59 e2
                                 ----- + -----
                                   2       2

                                 23 e1   61 e2
                                 ----- + -----
                                   2       2

                                 25 e1   61 e2
                                 ----- + -----
                                   2       2

                                 41 e1   45 e2
                                 ----- + -----
                                   2       2

                                 43 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   47 e2
                                 ----- + -----
                                   2       2

                                 45 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 43 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 39 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 35 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   57 e2
                                 ----- + -----
                                   2       2

                                 31 e1   57 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   57 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   57 e2
                                 ----- + -----
                                   2       2

                                 25 e1   59 e2
                                 ----- + -----
                                   2       2

                                 27 e1   59 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   57 e2
                                 ----- + -----
                                   2       2

                                 19 e1   59 e2
                                 ----- + -----
                                   2       2

                                 25 e1   57 e2
                                 ----- + -----
                                   2       2

                                 25 e1   59 e2
                                 ----- + -----
                                   2       2

                                 21 e1   61 e2
                                 ----- + -----
                                   2       2

                                 23 e1   61 e2
                                 ----- + -----
                                   2       2

                                 43 e1   43 e2
                                 ----- + -----
                                   2       2

                                 45 e1   43 e2
                                 ----- + -----
                                   2       2

                                 41 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 47 e1   45 e2
                                 ----- + -----
                                   2       2

                                 47 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 45 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 41 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 37 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 33 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   57 e2
                                 ----- + -----
                                   2       2

                                 29 e1   57 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   57 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   57 e2
                                 ----- + -----
                                   2       2

                                 23 e1   59 e2
                                 ----- + -----
                                   2       2

                                 25 e1   59 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   57 e2
                                 ----- + -----
                                   2       2

                                 17 e1   59 e2
                                 ----- + -----
                                   2       2

                                 23 e1   57 e2
                                 ----- + -----
                                   2       2

                                 23 e1   59 e2
                                 ----- + -----
                                   2       2

                                 19 e1   61 e2
                                 ----- + -----
                                   2       2

                                 21 e1   61 e2
                                 ----- + -----
                                   2       2

                                 41 e1   43 e2
                                 ----- + -----
                                   2       2

                                 43 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 45 e1   45 e2
                                 ----- + -----
                                   2       2

                                 45 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 43 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 39 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 35 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 31 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   57 e2
                                 ----- + -----
                                   2       2

                                 27 e1   57 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   57 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   57 e2
                                 ----- + -----
                                   2       2

                                 21 e1   59 e2
                                 ----- + -----
                                   2       2

                                 23 e1   59 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   57 e2
                                 ----- + -----
                                   2       2

                                 15 e1   59 e2
                                 ----- + -----
                                   2       2

                                 21 e1   57 e2
                                 ----- + -----
                                   2       2

                                 21 e1   59 e2
                                 ----- + -----
                                   2       2

                                 17 e1   61 e2
                                 ----- + -----
                                   2       2

                                 19 e1   61 e2
                                 ----- + -----
                                   2       2

                                 39 e1   43 e2
                                 ----- + -----
                                   2       2

                                 41 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   45 e2
                                 ----- + -----
                                   2       2

                                 43 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 41 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 37 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 33 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 29 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   57 e2
                                 ----- + -----
                                   2       2

                                 25 e1   57 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   57 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   57 e2
                                 ----- + -----
                                   2       2

                                 19 e1   59 e2
                                 ----- + -----
                                   2       2

                                 21 e1   59 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   57 e2
                                 ----- + -----
                                   2       2

                                 13 e1   59 e2
                                 ----- + -----
                                   2       2

                                 19 e1   57 e2
                                 ----- + -----
                                   2       2

                                 19 e1   59 e2
                                 ----- + -----
                                   2       2

                                 15 e1   61 e2
                                 ----- + -----
                                   2       2

                                 17 e1   61 e2
                                 ----- + -----
                                   2       2

                                 41 e1   41 e2
                                 ----- + -----
                                   2       2

                                 43 e1   41 e2
                                 ----- + -----
                                   2       2

                                 39 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 45 e1   43 e2
                                 ----- + -----
                                   2       2

                                 45 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 43 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 39 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 35 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 31 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 27 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   57 e2
                                 ----- + -----
                                   2       2

                                 23 e1   57 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   57 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   57 e2
                                 ----- + -----
                                   2       2

                                 17 e1   59 e2
                                 ----- + -----
                                   2       2

                                 19 e1   59 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 11 e1   57 e2
                                 ----- + -----
                                   2       2

                                 11 e1   59 e2
                                 ----- + -----
                                   2       2

                                 17 e1   57 e2
                                 ----- + -----
                                   2       2

                                 17 e1   59 e2
                                 ----- + -----
                                   2       2

                                 13 e1   61 e2
                                 ----- + -----
                                   2       2

                                 15 e1   61 e2
                                 ----- + -----
                                   2       2

                                 39 e1   41 e2
                                 ----- + -----
                                   2       2

                                 41 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 43 e1   43 e2
                                 ----- + -----
                                   2       2

                                 43 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 41 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 37 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 33 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 29 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 25 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   57 e2
                                 ----- + -----
                                   2       2

                                 21 e1   57 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   57 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   57 e2
                                 ----- + -----
                                   2       2

                                 15 e1   59 e2
                                 ----- + -----
                                   2       2

                                 17 e1   59 e2
                                 ----- + -----
                                   2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 9 e1   57 e2
                                 ---- + -----
                                  2       2

                                 9 e1   59 e2
                                 ---- + -----
                                  2       2

                                 15 e1   57 e2
                                 ----- + -----
                                   2       2

                                 15 e1   59 e2
                                 ----- + -----
                                   2       2

                                 11 e1   61 e2
                                 ----- + -----
                                   2       2

                                 13 e1   61 e2
                                 ----- + -----
                                   2       2

                                 37 e1   41 e2
                                 ----- + -----
                                   2       2

                                 39 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   43 e2
                                 ----- + -----
                                   2       2

                                 41 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 39 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 35 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 31 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 27 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 23 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   57 e2
                                 ----- + -----
                                   2       2

                                 19 e1   57 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 11 e1   57 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   57 e2
                                 ----- + -----
                                   2       2

                                 13 e1   59 e2
                                 ----- + -----
                                   2       2

                                 15 e1   59 e2
                                 ----- + -----
                                   2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 7 e1   57 e2
                                 ---- + -----
                                  2       2

                                 7 e1   59 e2
                                 ---- + -----
                                  2       2

                                 13 e1   57 e2
                                 ----- + -----
                                   2       2

                                 13 e1   59 e2
                                 ----- + -----
                                   2       2

                                 9 e1   61 e2
                                 ---- + -----
                                  2       2

                                 11 e1   61 e2
                                 ----- + -----
                                   2       2

                                 39 e1   39 e2
                                 ----- + -----
                                   2       2

                                 41 e1   39 e2
                                 ----- + -----
                                   2       2

                                 37 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 43 e1   41 e2
                                 ----- + -----
                                   2       2

                                 43 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 41 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 37 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 33 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 29 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 25 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 21 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   57 e2
                                 ----- + -----
                                   2       2

                                 17 e1   57 e2
                                 ----- + -----
                                   2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 9 e1   57 e2
                                 ---- + -----
                                  2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   57 e2
                                 ----- + -----
                                   2       2

                                 11 e1   59 e2
                                 ----- + -----
                                   2       2

                                 13 e1   59 e2
                                 ----- + -----
                                   2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 5 e1   57 e2
                                 ---- + -----
                                  2       2

                                 5 e1   59 e2
                                 ---- + -----
                                  2       2

                                 11 e1   57 e2
                                 ----- + -----
                                   2       2

                                 11 e1   59 e2
                                 ----- + -----
                                   2       2

                                 7 e1   61 e2
                                 ---- + -----
                                  2       2

                                 9 e1   61 e2
                                 ---- + -----
                                  2       2

                                 37 e1   39 e2
                                 ----- + -----
                                   2       2

                                 39 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 41 e1   41 e2
                                 ----- + -----
                                   2       2

                                 41 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 39 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 35 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 31 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 27 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 23 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 19 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   57 e2
                                 ----- + -----
                                   2       2

                                 15 e1   57 e2
                                 ----- + -----
                                   2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 7 e1   57 e2
                                 ---- + -----
                                  2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   57 e2
                                 ----- + -----
                                   2       2

                                 9 e1   59 e2
                                 ---- + -----
                                  2       2

                                 11 e1   59 e2
                                 ----- + -----
                                   2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 3 e1   57 e2
                                 ---- + -----
                                  2       2

                                 3 e1   59 e2
                                 ---- + -----
                                  2       2

                                 9 e1   57 e2
                                 ---- + -----
                                  2       2

                                 9 e1   59 e2
                                 ---- + -----
                                  2       2

                                 5 e1   61 e2
                                 ---- + -----
                                  2       2

                                 7 e1   61 e2
                                 ---- + -----
                                  2       2

                                 35 e1   39 e2
                                 ----- + -----
                                   2       2

                                 37 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   41 e2
                                 ----- + -----
                                   2       2

                                 39 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 37 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 33 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 29 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 25 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 21 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 17 e1   55 e2
                                 ----- + -----
                                   2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 11 e1   57 e2
                                 ----- + -----
                                   2       2

                                 13 e1   57 e2
                                 ----- + -----
                                   2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                 5 e1   57 e2
                                 ---- + -----
                                  2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 11 e1   57 e2
                                 ----- + -----
                                   2       2

                                 7 e1   59 e2
                                 ---- + -----
                                  2       2

                                 9 e1   59 e2
                                 ---- + -----
                                  2       2

                                 3 e1   55 e2
                                 ---- + -----
                                  2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                  e1    57 e2
                                 ---- + -----
                                  2       2

                                  e1    59 e2
                                 ---- + -----
                                  2       2

                                 7 e1   57 e2
                                 ---- + -----
                                  2       2

                                 7 e1   59 e2
                                 ---- + -----
                                  2       2

                                 3 e1   61 e2
                                 ---- + -----
                                  2       2

                                 5 e1   61 e2
                                 ---- + -----
                                  2       2

                                 37 e1   37 e2
                                 ----- + -----
                                   2       2

                                 39 e1   37 e2
                                 ----- + -----
                                   2       2

                                 35 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 41 e1   39 e2
                                 ----- + -----
                                   2       2

                                 41 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 39 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 35 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 31 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 27 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 23 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 19 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 15 e1   55 e2
                                 ----- + -----
                                   2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 9 e1   57 e2
                                 ---- + -----
                                  2       2

                                 11 e1   57 e2
                                 ----- + -----
                                   2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 3 e1   55 e2
                                 ---- + -----
                                  2       2

                                 3 e1   57 e2
                                 ---- + -----
                                  2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 9 e1   57 e2
                                 ---- + -----
                                  2       2

                                 5 e1   59 e2
                                 ---- + -----
                                  2       2

                                 7 e1   59 e2
                                 ---- + -----
                                  2       2

                                  e1    55 e2
                                 ---- + -----
                                  2       2

                                 3 e1   55 e2
                                 ---- + -----
                                  2       2

                                   e1    57 e2
                                - ---- + -----
                                   2       2

                                   e1    59 e2
                                - ---- + -----
                                   2       2

                                 5 e1   57 e2
                                 ---- + -----
                                  2       2

                                 5 e1   59 e2
                                 ---- + -----
                                  2       2

                                  e1    61 e2
                                 ---- + -----
                                  2       2

                                 3 e1   61 e2
                                 ---- + -----
                                  2       2

                                 35 e1   37 e2
                                 ----- + -----
                                   2       2

                                 37 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 39 e1   39 e2
                                 ----- + -----
                                   2       2

                                 39 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 37 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 33 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 29 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 25 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 21 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 17 e1   53 e2
                                 ----- + -----
                                   2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 13 e1   55 e2
                                 ----- + -----
                                   2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 7 e1   57 e2
                                 ---- + -----
                                  2       2

                                 9 e1   57 e2
                                 ---- + -----
                                  2       2

                                 3 e1   53 e2
                                 ---- + -----
                                  2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                  e1    55 e2
                                 ---- + -----
                                  2       2

                                  e1    57 e2
                                 ---- + -----
                                  2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 7 e1   57 e2
                                 ---- + -----
                                  2       2

                                 3 e1   59 e2
                                 ---- + -----
                                  2       2

                                 5 e1   59 e2
                                 ---- + -----
                                  2       2

                                 33 e1   37 e2
                                 ----- + -----
                                   2       2

                                 35 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   39 e2
                                 ----- + -----
                                   2       2

                                 37 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 35 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 31 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 27 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 23 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 19 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 15 e1   53 e2
                                 ----- + -----
                                   2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 11 e1   55 e2
                                 ----- + -----
                                   2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 3 e1   53 e2
                                 ---- + -----
                                  2       2

                                 3 e1   55 e2
                                 ---- + -----
                                  2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 5 e1   57 e2
                                 ---- + -----
                                  2       2

                                 7 e1   57 e2
                                 ---- + -----
                                  2       2

                                  e1    53 e2
                                 ---- + -----
                                  2       2

                                 3 e1   53 e2
                                 ---- + -----
                                  2       2

                                   e1    55 e2
                                - ---- + -----
                                   2       2

                                   e1    57 e2
                                - ---- + -----
                                   2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                 5 e1   57 e2
                                 ---- + -----
                                  2       2

                                  e1    59 e2
                                 ---- + -----
                                  2       2

                                 3 e1   59 e2
                                 ---- + -----
                                  2       2

                                 35 e1   35 e2
                                 ----- + -----
                                   2       2

                                 37 e1   35 e2
                                 ----- + -----
                                   2       2

                                 33 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 39 e1   37 e2
                                 ----- + -----
                                   2       2

                                 39 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 37 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 33 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 29 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 25 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 21 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 17 e1   51 e2
                                 ----- + -----
                                   2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 13 e1   53 e2
                                 ----- + -----
                                   2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 9 e1   55 e2
                                 ---- + -----
                                  2       2

                                 3 e1   51 e2
                                 ---- + -----
                                  2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                  e1    53 e2
                                 ---- + -----
                                  2       2

                                  e1    55 e2
                                 ---- + -----
                                  2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                 3 e1   57 e2
                                 ---- + -----
                                  2       2

                                 5 e1   57 e2
                                 ---- + -----
                                  2       2

                                 33 e1   35 e2
                                 ----- + -----
                                   2       2

                                 35 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 37 e1   37 e2
                                 ----- + -----
                                   2       2

                                 37 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 35 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 31 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 27 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 23 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 19 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 15 e1   51 e2
                                 ----- + -----
                                   2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 11 e1   53 e2
                                 ----- + -----
                                   2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 3 e1   51 e2
                                 ---- + -----
                                  2       2

                                 3 e1   53 e2
                                 ---- + -----
                                  2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                 7 e1   55 e2
                                 ---- + -----
                                  2       2

                                  e1    51 e2
                                 ---- + -----
                                  2       2

                                 3 e1   51 e2
                                 ---- + -----
                                  2       2

                                   e1    53 e2
                                - ---- + -----
                                   2       2

                                   e1    55 e2
                                - ---- + -----
                                   2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                  e1    57 e2
                                 ---- + -----
                                  2       2

                                 3 e1   57 e2
                                 ---- + -----
                                  2       2

                                 31 e1   35 e2
                                 ----- + -----
                                   2       2

                                 33 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   37 e2
                                 ----- + -----
                                   2       2

                                 35 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 33 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 29 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 25 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 21 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 17 e1   49 e2
                                 ----- + -----
                                   2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 13 e1   51 e2
                                 ----- + -----
                                   2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 9 e1   53 e2
                                 ---- + -----
                                  2       2

                                 3 e1   49 e2
                                 ---- + -----
                                  2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                  e1    51 e2
                                 ---- + -----
                                  2       2

                                  e1    53 e2
                                 ---- + -----
                                  2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                 3 e1   55 e2
                                 ---- + -----
                                  2       2

                                 5 e1   55 e2
                                 ---- + -----
                                  2       2

                                 33 e1   33 e2
                                 ----- + -----
                                   2       2

                                 35 e1   33 e2
                                 ----- + -----
                                   2       2

                                 31 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 37 e1   35 e2
                                 ----- + -----
                                   2       2

                                 37 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 35 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 31 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 27 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 23 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 19 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 15 e1   49 e2
                                 ----- + -----
                                   2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 11 e1   51 e2
                                 ----- + -----
                                   2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 3 e1   49 e2
                                 ---- + -----
                                  2       2

                                 3 e1   51 e2
                                 ---- + -----
                                  2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                 7 e1   53 e2
                                 ---- + -----
                                  2       2

                                  e1    49 e2
                                 ---- + -----
                                  2       2

                                 3 e1   49 e2
                                 ---- + -----
                                  2       2

                                   e1    51 e2
                                - ---- + -----
                                   2       2

                                   e1    53 e2
                                - ---- + -----
                                   2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                  e1    55 e2
                                 ---- + -----
                                  2       2

                                 3 e1   55 e2
                                 ---- + -----
                                  2       2

                                 31 e1   33 e2
                                 ----- + -----
                                   2       2

                                 33 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 35 e1   35 e2
                                 ----- + -----
                                   2       2

                                 35 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 33 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 29 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 25 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 21 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 17 e1   47 e2
                                 ----- + -----
                                   2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 13 e1   49 e2
                                 ----- + -----
                                   2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 9 e1   51 e2
                                 ---- + -----
                                  2       2

                                 3 e1   47 e2
                                 ---- + -----
                                  2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                  e1    49 e2
                                 ---- + -----
                                  2       2

                                  e1    51 e2
                                 ---- + -----
                                  2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                 3 e1   53 e2
                                 ---- + -----
                                  2       2

                                 5 e1   53 e2
                                 ---- + -----
                                  2       2

                                 29 e1   33 e2
                                 ----- + -----
                                   2       2

                                 31 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   35 e2
                                 ----- + -----
                                   2       2

                                 33 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 31 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 27 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 23 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 19 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 15 e1   47 e2
                                 ----- + -----
                                   2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 11 e1   49 e2
                                 ----- + -----
                                   2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 3 e1   47 e2
                                 ---- + -----
                                  2       2

                                 3 e1   49 e2
                                 ---- + -----
                                  2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                 7 e1   51 e2
                                 ---- + -----
                                  2       2

                                  e1    47 e2
                                 ---- + -----
                                  2       2

                                 3 e1   47 e2
                                 ---- + -----
                                  2       2

                                   e1    49 e2
                                - ---- + -----
                                   2       2

                                   e1    51 e2
                                - ---- + -----
                                   2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                  e1    53 e2
                                 ---- + -----
                                  2       2

                                 3 e1   53 e2
                                 ---- + -----
                                  2       2

                                 31 e1   31 e2
                                 ----- + -----
                                   2       2

                                 33 e1   31 e2
                                 ----- + -----
                                   2       2

                                 29 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 35 e1   33 e2
                                 ----- + -----
                                   2       2

                                 35 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 33 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 29 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 25 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 21 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 17 e1   45 e2
                                 ----- + -----
                                   2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 13 e1   47 e2
                                 ----- + -----
                                   2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 9 e1   49 e2
                                 ---- + -----
                                  2       2

                                 3 e1   45 e2
                                 ---- + -----
                                  2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                  e1    47 e2
                                 ---- + -----
                                  2       2

                                  e1    49 e2
                                 ---- + -----
                                  2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                 3 e1   51 e2
                                 ---- + -----
                                  2       2

                                 5 e1   51 e2
                                 ---- + -----
                                  2       2

                                 29 e1   31 e2
                                 ----- + -----
                                   2       2

                                 31 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 33 e1   33 e2
                                 ----- + -----
                                   2       2

                                 33 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 31 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 27 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 23 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 19 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 15 e1   45 e2
                                 ----- + -----
                                   2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 11 e1   47 e2
                                 ----- + -----
                                   2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 3 e1   45 e2
                                 ---- + -----
                                  2       2

                                 3 e1   47 e2
                                 ---- + -----
                                  2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                 7 e1   49 e2
                                 ---- + -----
                                  2       2

                                  e1    45 e2
                                 ---- + -----
                                  2       2

                                 3 e1   45 e2
                                 ---- + -----
                                  2       2

                                   e1    47 e2
                                - ---- + -----
                                   2       2

                                   e1    49 e2
                                - ---- + -----
                                   2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                  e1    51 e2
                                 ---- + -----
                                  2       2

                                 3 e1   51 e2
                                 ---- + -----
                                  2       2

                                 27 e1   31 e2
                                 ----- + -----
                                   2       2

                                 29 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   33 e2
                                 ----- + -----
                                   2       2

                                 31 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 29 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 25 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 21 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 17 e1   43 e2
                                 ----- + -----
                                   2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 13 e1   45 e2
                                 ----- + -----
                                   2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 9 e1   47 e2
                                 ---- + -----
                                  2       2

                                 3 e1   43 e2
                                 ---- + -----
                                  2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                  e1    45 e2
                                 ---- + -----
                                  2       2

                                  e1    47 e2
                                 ---- + -----
                                  2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                 3 e1   49 e2
                                 ---- + -----
                                  2       2

                                 5 e1   49 e2
                                 ---- + -----
                                  2       2

                                 29 e1   29 e2
                                 ----- + -----
                                   2       2

                                 31 e1   29 e2
                                 ----- + -----
                                   2       2

                                 27 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 33 e1   31 e2
                                 ----- + -----
                                   2       2

                                 33 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 31 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 27 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 23 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 19 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 15 e1   43 e2
                                 ----- + -----
                                   2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 11 e1   45 e2
                                 ----- + -----
                                   2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 3 e1   43 e2
                                 ---- + -----
                                  2       2

                                 3 e1   45 e2
                                 ---- + -----
                                  2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                 7 e1   47 e2
                                 ---- + -----
                                  2       2

                                  e1    43 e2
                                 ---- + -----
                                  2       2

                                 3 e1   43 e2
                                 ---- + -----
                                  2       2

                                   e1    45 e2
                                - ---- + -----
                                   2       2

                                   e1    47 e2
                                - ---- + -----
                                   2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                  e1    49 e2
                                 ---- + -----
                                  2       2

                                 3 e1   49 e2
                                 ---- + -----
                                  2       2

                                 27 e1   29 e2
                                 ----- + -----
                                   2       2

                                 29 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 31 e1   31 e2
                                 ----- + -----
                                   2       2

                                 31 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 29 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 25 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 21 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 17 e1   41 e2
                                 ----- + -----
                                   2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 13 e1   43 e2
                                 ----- + -----
                                   2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 9 e1   45 e2
                                 ---- + -----
                                  2       2

                                 3 e1   41 e2
                                 ---- + -----
                                  2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                  e1    43 e2
                                 ---- + -----
                                  2       2

                                  e1    45 e2
                                 ---- + -----
                                  2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                 3 e1   47 e2
                                 ---- + -----
                                  2       2

                                 5 e1   47 e2
                                 ---- + -----
                                  2       2

                                 25 e1   29 e2
                                 ----- + -----
                                   2       2

                                 27 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   31 e2
                                 ----- + -----
                                   2       2

                                 29 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 27 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 23 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 19 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 15 e1   41 e2
                                 ----- + -----
                                   2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 11 e1   43 e2
                                 ----- + -----
                                   2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 3 e1   41 e2
                                 ---- + -----
                                  2       2

                                 3 e1   43 e2
                                 ---- + -----
                                  2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                 7 e1   45 e2
                                 ---- + -----
                                  2       2

                                  e1    41 e2
                                 ---- + -----
                                  2       2

                                 3 e1   41 e2
                                 ---- + -----
                                  2       2

                                   e1    43 e2
                                - ---- + -----
                                   2       2

                                   e1    45 e2
                                - ---- + -----
                                   2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                  e1    47 e2
                                 ---- + -----
                                  2       2

                                 3 e1   47 e2
                                 ---- + -----
                                  2       2

                                 27 e1   27 e2
                                 ----- + -----
                                   2       2

                                 29 e1   27 e2
                                 ----- + -----
                                   2       2

                                 25 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 31 e1   29 e2
                                 ----- + -----
                                   2       2

                                 31 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 29 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 25 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 21 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 17 e1   39 e2
                                 ----- + -----
                                   2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 13 e1   41 e2
                                 ----- + -----
                                   2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 9 e1   43 e2
                                 ---- + -----
                                  2       2

                                 3 e1   39 e2
                                 ---- + -----
                                  2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                  e1    41 e2
                                 ---- + -----
                                  2       2

                                  e1    43 e2
                                 ---- + -----
                                  2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                 3 e1   45 e2
                                 ---- + -----
                                  2       2

                                 5 e1   45 e2
                                 ---- + -----
                                  2       2

                                 25 e1   27 e2
                                 ----- + -----
                                   2       2

                                 27 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 29 e1   29 e2
                                 ----- + -----
                                   2       2

                                 29 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 27 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 23 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 19 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 15 e1   39 e2
                                 ----- + -----
                                   2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 11 e1   41 e2
                                 ----- + -----
                                   2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 3 e1   39 e2
                                 ---- + -----
                                  2       2

                                 3 e1   41 e2
                                 ---- + -----
                                  2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                 7 e1   43 e2
                                 ---- + -----
                                  2       2

                                  e1    39 e2
                                 ---- + -----
                                  2       2

                                 3 e1   39 e2
                                 ---- + -----
                                  2       2

                                   e1    41 e2
                                - ---- + -----
                                   2       2

                                   e1    43 e2
                                - ---- + -----
                                   2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                  e1    45 e2
                                 ---- + -----
                                  2       2

                                 3 e1   45 e2
                                 ---- + -----
                                  2       2

                                 23 e1   27 e2
                                 ----- + -----
                                   2       2

                                 25 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   29 e2
                                 ----- + -----
                                   2       2

                                 27 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 25 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 21 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 17 e1   37 e2
                                 ----- + -----
                                   2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 13 e1   39 e2
                                 ----- + -----
                                   2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 9 e1   41 e2
                                 ---- + -----
                                  2       2

                                 3 e1   37 e2
                                 ---- + -----
                                  2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                  e1    39 e2
                                 ---- + -----
                                  2       2

                                  e1    41 e2
                                 ---- + -----
                                  2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                 3 e1   43 e2
                                 ---- + -----
                                  2       2

                                 5 e1   43 e2
                                 ---- + -----
                                  2       2

                                 25 e1   25 e2
                                 ----- + -----
                                   2       2

                                 27 e1   25 e2
                                 ----- + -----
                                   2       2

                                 23 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 29 e1   27 e2
                                 ----- + -----
                                   2       2

                                 29 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 27 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 23 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 19 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 15 e1   37 e2
                                 ----- + -----
                                   2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 11 e1   39 e2
                                 ----- + -----
                                   2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 3 e1   37 e2
                                 ---- + -----
                                  2       2

                                 3 e1   39 e2
                                 ---- + -----
                                  2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                 7 e1   41 e2
                                 ---- + -----
                                  2       2

                                  e1    37 e2
                                 ---- + -----
                                  2       2

                                 3 e1   37 e2
                                 ---- + -----
                                  2       2

                                   e1    39 e2
                                - ---- + -----
                                   2       2

                                   e1    41 e2
                                - ---- + -----
                                   2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                  e1    43 e2
                                 ---- + -----
                                  2       2

                                 3 e1   43 e2
                                 ---- + -----
                                  2       2

                                 23 e1   25 e2
                                 ----- + -----
                                   2       2

                                 25 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 27 e1   27 e2
                                 ----- + -----
                                   2       2

                                 27 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 25 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 21 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 17 e1   35 e2
                                 ----- + -----
                                   2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 13 e1   37 e2
                                 ----- + -----
                                   2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 9 e1   39 e2
                                 ---- + -----
                                  2       2

                                 3 e1   35 e2
                                 ---- + -----
                                  2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                  e1    37 e2
                                 ---- + -----
                                  2       2

                                  e1    39 e2
                                 ---- + -----
                                  2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                 3 e1   41 e2
                                 ---- + -----
                                  2       2

                                 5 e1   41 e2
                                 ---- + -----
                                  2       2

                                 21 e1   25 e2
                                 ----- + -----
                                   2       2

                                 23 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   27 e2
                                 ----- + -----
                                   2       2

                                 25 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 23 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 19 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 15 e1   35 e2
                                 ----- + -----
                                   2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 11 e1   37 e2
                                 ----- + -----
                                   2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 3 e1   35 e2
                                 ---- + -----
                                  2       2

                                 3 e1   37 e2
                                 ---- + -----
                                  2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                 7 e1   39 e2
                                 ---- + -----
                                  2       2

                                  e1    35 e2
                                 ---- + -----
                                  2       2

                                 3 e1   35 e2
                                 ---- + -----
                                  2       2

                                   e1    37 e2
                                - ---- + -----
                                   2       2

                                   e1    39 e2
                                - ---- + -----
                                   2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                  e1    41 e2
                                 ---- + -----
                                  2       2

                                 3 e1   41 e2
                                 ---- + -----
                                  2       2

                                 23 e1   23 e2
                                 ----- + -----
                                   2       2

                                 25 e1   23 e2
                                 ----- + -----
                                   2       2

                                 21 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 27 e1   25 e2
                                 ----- + -----
                                   2       2

                                 27 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 25 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 21 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 17 e1   33 e2
                                 ----- + -----
                                   2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 13 e1   35 e2
                                 ----- + -----
                                   2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 9 e1   37 e2
                                 ---- + -----
                                  2       2

                                 3 e1   33 e2
                                 ---- + -----
                                  2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                  e1    35 e2
                                 ---- + -----
                                  2       2

                                  e1    37 e2
                                 ---- + -----
                                  2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                 3 e1   39 e2
                                 ---- + -----
                                  2       2

                                 5 e1   39 e2
                                 ---- + -----
                                  2       2

                                 21 e1   23 e2
                                 ----- + -----
                                   2       2

                                 23 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 25 e1   25 e2
                                 ----- + -----
                                   2       2

                                 25 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 23 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 19 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 15 e1   33 e2
                                 ----- + -----
                                   2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 11 e1   35 e2
                                 ----- + -----
                                   2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 3 e1   33 e2
                                 ---- + -----
                                  2       2

                                 3 e1   35 e2
                                 ---- + -----
                                  2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                 7 e1   37 e2
                                 ---- + -----
                                  2       2

                                  e1    33 e2
                                 ---- + -----
                                  2       2

                                 3 e1   33 e2
                                 ---- + -----
                                  2       2

                                   e1    35 e2
                                - ---- + -----
                                   2       2

                                   e1    37 e2
                                - ---- + -----
                                   2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                  e1    39 e2
                                 ---- + -----
                                  2       2

                                 3 e1   39 e2
                                 ---- + -----
                                  2       2

                                 19 e1   23 e2
                                 ----- + -----
                                   2       2

                                 21 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   25 e2
                                 ----- + -----
                                   2       2

                                 23 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 21 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 17 e1   31 e2
                                 ----- + -----
                                   2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 13 e1   33 e2
                                 ----- + -----
                                   2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 9 e1   35 e2
                                 ---- + -----
                                  2       2

                                 3 e1   31 e2
                                 ---- + -----
                                  2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                  e1    33 e2
                                 ---- + -----
                                  2       2

                                  e1    35 e2
                                 ---- + -----
                                  2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                 3 e1   37 e2
                                 ---- + -----
                                  2       2

                                 5 e1   37 e2
                                 ---- + -----
                                  2       2

                                 21 e1   21 e2
                                 ----- + -----
                                   2       2

                                 23 e1   21 e2
                                 ----- + -----
                                   2       2

                                 19 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 25 e1   23 e2
                                 ----- + -----
                                   2       2

                                 25 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 23 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   23 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 19 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 15 e1   31 e2
                                 ----- + -----
                                   2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 11 e1   33 e2
                                 ----- + -----
                                   2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 3 e1   31 e2
                                 ---- + -----
                                  2       2

                                 3 e1   33 e2
                                 ---- + -----
                                  2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                 7 e1   35 e2
                                 ---- + -----
                                  2       2

                                  e1    31 e2
                                 ---- + -----
                                  2       2

                                 3 e1   31 e2
                                 ---- + -----
                                  2       2

                                   e1    33 e2
                                - ---- + -----
                                   2       2

                                   e1    35 e2
                                - ---- + -----
                                   2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                  e1    37 e2
                                 ---- + -----
                                  2       2

                                 3 e1   37 e2
                                 ---- + -----
                                  2       2

                                 19 e1   21 e2
                                 ----- + -----
                                   2       2

                                 21 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 23 e1   23 e2
                                 ----- + -----
                                   2       2

                                 23 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 21 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 17 e1   29 e2
                                 ----- + -----
                                   2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 13 e1   31 e2
                                 ----- + -----
                                   2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 9 e1   33 e2
                                 ---- + -----
                                  2       2

                                 3 e1   29 e2
                                 ---- + -----
                                  2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                  e1    31 e2
                                 ---- + -----
                                  2       2

                                  e1    33 e2
                                 ---- + -----
                                  2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                 3 e1   35 e2
                                 ---- + -----
                                  2       2

                                 5 e1   35 e2
                                 ---- + -----
                                  2       2

                                 17 e1   21 e2
                                 ----- + -----
                                   2       2

                                 19 e1   21 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   23 e2
                                 ----- + -----
                                   2       2

                                 21 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 19 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 15 e1   29 e2
                                 ----- + -----
                                   2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 11 e1   31 e2
                                 ----- + -----
                                   2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 3 e1   29 e2
                                 ---- + -----
                                  2       2

                                 3 e1   31 e2
                                 ---- + -----
                                  2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                 7 e1   33 e2
                                 ---- + -----
                                  2       2

                                  e1    29 e2
                                 ---- + -----
                                  2       2

                                 3 e1   29 e2
                                 ---- + -----
                                  2       2

                                   e1    31 e2
                                - ---- + -----
                                   2       2

                                   e1    33 e2
                                - ---- + -----
                                   2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                  e1    35 e2
                                 ---- + -----
                                  2       2

                                 3 e1   35 e2
                                 ---- + -----
                                  2       2

                                 19 e1   19 e2
                                 ----- + -----
                                   2       2

                                 21 e1   19 e2
                                 ----- + -----
                                   2       2

                                 17 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 23 e1   21 e2
                                 ----- + -----
                                   2       2

                                 23 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 21 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   21 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 17 e1   27 e2
                                 ----- + -----
                                   2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 13 e1   29 e2
                                 ----- + -----
                                   2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 9 e1   31 e2
                                 ---- + -----
                                  2       2

                                 3 e1   27 e2
                                 ---- + -----
                                  2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                  e1    29 e2
                                 ---- + -----
                                  2       2

                                  e1    31 e2
                                 ---- + -----
                                  2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                 3 e1   33 e2
                                 ---- + -----
                                  2       2

                                 5 e1   33 e2
                                 ---- + -----
                                  2       2

                                 17 e1   19 e2
                                 ----- + -----
                                   2       2

                                 19 e1   19 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 21 e1   21 e2
                                 ----- + -----
                                   2       2

                                 21 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 19 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 15 e1   27 e2
                                 ----- + -----
                                   2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 11 e1   29 e2
                                 ----- + -----
                                   2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 3 e1   27 e2
                                 ---- + -----
                                  2       2

                                 3 e1   29 e2
                                 ---- + -----
                                  2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                 7 e1   31 e2
                                 ---- + -----
                                  2       2

                                  e1    27 e2
                                 ---- + -----
                                  2       2

                                 3 e1   27 e2
                                 ---- + -----
                                  2       2

                                   e1    29 e2
                                - ---- + -----
                                   2       2

                                   e1    31 e2
                                - ---- + -----
                                   2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                  e1    33 e2
                                 ---- + -----
                                  2       2

                                 3 e1   33 e2
                                 ---- + -----
                                  2       2

                                 15 e1   19 e2
                                 ----- + -----
                                   2       2

                                 17 e1   19 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   21 e2
                                 ----- + -----
                                   2       2

                                 19 e1   23 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 17 e1   25 e2
                                 ----- + -----
                                   2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 13 e1   27 e2
                                 ----- + -----
                                   2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 9 e1   29 e2
                                 ---- + -----
                                  2       2

                                 3 e1   25 e2
                                 ---- + -----
                                  2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                  e1    27 e2
                                 ---- + -----
                                  2       2

                                  e1    29 e2
                                 ---- + -----
                                  2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                 3 e1   31 e2
                                 ---- + -----
                                  2       2

                                 5 e1   31 e2
                                 ---- + -----
                                  2       2

                                 17 e1   17 e2
                                 ----- + -----
                                   2       2

                                 19 e1   17 e2
                                 ----- + -----
                                   2       2

                                 15 e1   19 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 21 e1   19 e2
                                 ----- + -----
                                   2       2

                                 21 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 19 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 15 e1   19 e2
                                 ----- + -----
                                   2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 15 e1   25 e2
                                 ----- + -----
                                   2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 11 e1   27 e2
                                 ----- + -----
                                   2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 3 e1   25 e2
                                 ---- + -----
                                  2       2

                                 3 e1   27 e2
                                 ---- + -----
                                  2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                 7 e1   29 e2
                                 ---- + -----
                                  2       2

                                  e1    25 e2
                                 ---- + -----
                                  2       2

                                 3 e1   25 e2
                                 ---- + -----
                                  2       2

                                   e1    27 e2
                                - ---- + -----
                                   2       2

                                   e1    29 e2
                                - ---- + -----
                                   2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                  e1    31 e2
                                 ---- + -----
                                  2       2

                                 3 e1   31 e2
                                 ---- + -----
                                  2       2

                                 15 e1   17 e2
                                 ----- + -----
                                   2       2

                                 17 e1   17 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 19 e1   19 e2
                                 ----- + -----
                                   2       2

                                 19 e1   21 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 17 e1   23 e2
                                 ----- + -----
                                   2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 13 e1   25 e2
                                 ----- + -----
                                   2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 9 e1   27 e2
                                 ---- + -----
                                  2       2

                                 3 e1   23 e2
                                 ---- + -----
                                  2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                  e1    25 e2
                                 ---- + -----
                                  2       2

                                  e1    27 e2
                                 ---- + -----
                                  2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                 3 e1   29 e2
                                 ---- + -----
                                  2       2

                                 5 e1   29 e2
                                 ---- + -----
                                  2       2

                                 13 e1   17 e2
                                 ----- + -----
                                   2       2

                                 15 e1   17 e2
                                 ----- + -----
                                   2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   19 e2
                                 ----- + -----
                                   2       2

                                 17 e1   21 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 15 e1   23 e2
                                 ----- + -----
                                   2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 11 e1   25 e2
                                 ----- + -----
                                   2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 3 e1   23 e2
                                 ---- + -----
                                  2       2

                                 3 e1   25 e2
                                 ---- + -----
                                  2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                 7 e1   27 e2
                                 ---- + -----
                                  2       2

                                  e1    23 e2
                                 ---- + -----
                                  2       2

                                 3 e1   23 e2
                                 ---- + -----
                                  2       2

                                   e1    25 e2
                                - ---- + -----
                                   2       2

                                   e1    27 e2
                                - ---- + -----
                                   2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                  e1    29 e2
                                 ---- + -----
                                  2       2

                                 3 e1   29 e2
                                 ---- + -----
                                  2       2

                                 15 e1   15 e2
                                 ----- + -----
                                   2       2

                                 17 e1   15 e2
                                 ----- + -----
                                   2       2

                                 13 e1   17 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 19 e1   17 e2
                                 ----- + -----
                                   2       2

                                 19 e1   19 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 17 e1   21 e2
                                 ----- + -----
                                   2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 13 e1   17 e2
                                 ----- + -----
                                   2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 15 e1   19 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 13 e1   23 e2
                                 ----- + -----
                                   2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 9 e1   25 e2
                                 ---- + -----
                                  2       2

                                 3 e1   21 e2
                                 ---- + -----
                                  2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                  e1    23 e2
                                 ---- + -----
                                  2       2

                                  e1    25 e2
                                 ---- + -----
                                  2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                 3 e1   27 e2
                                 ---- + -----
                                  2       2

                                 5 e1   27 e2
                                 ---- + -----
                                  2       2

                                 13 e1   15 e2
                                 ----- + -----
                                   2       2

                                 15 e1   15 e2
                                 ----- + -----
                                   2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 17 e1   17 e2
                                 ----- + -----
                                   2       2

                                 17 e1   19 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 15 e1   21 e2
                                 ----- + -----
                                   2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 11 e1   23 e2
                                 ----- + -----
                                   2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 3 e1   21 e2
                                 ---- + -----
                                  2       2

                                 3 e1   23 e2
                                 ---- + -----
                                  2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                 7 e1   25 e2
                                 ---- + -----
                                  2       2

                                  e1    21 e2
                                 ---- + -----
                                  2       2

                                 3 e1   21 e2
                                 ---- + -----
                                  2       2

                                   e1    23 e2
                                - ---- + -----
                                   2       2

                                   e1    25 e2
                                - ---- + -----
                                   2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                  e1    27 e2
                                 ---- + -----
                                  2       2

                                 3 e1   27 e2
                                 ---- + -----
                                  2       2

                                 11 e1   15 e2
                                 ----- + -----
                                   2       2

                                 13 e1   15 e2
                                 ----- + -----
                                   2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 15 e1   17 e2
                                 ----- + -----
                                   2       2

                                 15 e1   19 e2
                                 ----- + -----
                                   2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 13 e1   21 e2
                                 ----- + -----
                                   2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 9 e1   23 e2
                                 ---- + -----
                                  2       2

                                 3 e1   19 e2
                                 ---- + -----
                                  2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                  e1    21 e2
                                 ---- + -----
                                  2       2

                                  e1    23 e2
                                 ---- + -----
                                  2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                 3 e1   25 e2
                                 ---- + -----
                                  2       2

                                 5 e1   25 e2
                                 ---- + -----
                                  2       2

                                 13 e1   13 e2
                                 ----- + -----
                                   2       2

                                 15 e1   13 e2
                                 ----- + -----
                                   2       2

                                 11 e1   15 e2
                                 ----- + -----
                                   2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 17 e1   15 e2
                                 ----- + -----
                                   2       2

                                 17 e1   17 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 15 e1   19 e2
                                 ----- + -----
                                   2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 11 e1   15 e2
                                 ----- + -----
                                   2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 13 e1   17 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 11 e1   21 e2
                                 ----- + -----
                                   2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 3 e1   19 e2
                                 ---- + -----
                                  2       2

                                 3 e1   21 e2
                                 ---- + -----
                                  2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                 7 e1   23 e2
                                 ---- + -----
                                  2       2

                                  e1    19 e2
                                 ---- + -----
                                  2       2

                                 3 e1   19 e2
                                 ---- + -----
                                  2       2

                                   e1    21 e2
                                - ---- + -----
                                   2       2

                                   e1    23 e2
                                - ---- + -----
                                   2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                  e1    25 e2
                                 ---- + -----
                                  2       2

                                 3 e1   25 e2
                                 ---- + -----
                                  2       2

                                 11 e1   13 e2
                                 ----- + -----
                                   2       2

                                 13 e1   13 e2
                                 ----- + -----
                                   2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 15 e1   15 e2
                                 ----- + -----
                                   2       2

                                 15 e1   17 e2
                                 ----- + -----
                                   2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 13 e1   19 e2
                                 ----- + -----
                                   2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 9 e1   21 e2
                                 ---- + -----
                                  2       2

                                 3 e1   17 e2
                                 ---- + -----
                                  2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                  e1    19 e2
                                 ---- + -----
                                  2       2

                                  e1    21 e2
                                 ---- + -----
                                  2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                 3 e1   23 e2
                                 ---- + -----
                                  2       2

                                 5 e1   23 e2
                                 ---- + -----
                                  2       2

                                 9 e1   13 e2
                                 ---- + -----
                                  2       2

                                 11 e1   13 e2
                                 ----- + -----
                                   2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 13 e1   15 e2
                                 ----- + -----
                                   2       2

                                 13 e1   17 e2
                                 ----- + -----
                                   2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 11 e1   19 e2
                                 ----- + -----
                                   2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 3 e1   17 e2
                                 ---- + -----
                                  2       2

                                 3 e1   19 e2
                                 ---- + -----
                                  2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                 7 e1   21 e2
                                 ---- + -----
                                  2       2

                                  e1    17 e2
                                 ---- + -----
                                  2       2

                                 3 e1   17 e2
                                 ---- + -----
                                  2       2

                                   e1    19 e2
                                - ---- + -----
                                   2       2

                                   e1    21 e2
                                - ---- + -----
                                   2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                  e1    23 e2
                                 ---- + -----
                                  2       2

                                 3 e1   23 e2
                                 ---- + -----
                                  2       2

                                 11 e1   11 e2
                                 ----- + -----
                                   2       2

                                 13 e1   11 e2
                                 ----- + -----
                                   2       2

                                 9 e1   13 e2
                                 ---- + -----
                                  2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 15 e1   13 e2
                                 ----- + -----
                                   2       2

                                 15 e1   15 e2
                                 ----- + -----
                                   2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 13 e1   17 e2
                                 ----- + -----
                                   2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 9 e1   13 e2
                                 ---- + -----
                                  2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                 11 e1   15 e2
                                 ----- + -----
                                   2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 9 e1   19 e2
                                 ---- + -----
                                  2       2

                                 3 e1   15 e2
                                 ---- + -----
                                  2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                  e1    17 e2
                                 ---- + -----
                                  2       2

                                  e1    19 e2
                                 ---- + -----
                                  2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                 3 e1   21 e2
                                 ---- + -----
                                  2       2

                                 5 e1   21 e2
                                 ---- + -----
                                  2       2

                                 9 e1   11 e2
                                 ---- + -----
                                  2       2

                                 11 e1   11 e2
                                 ----- + -----
                                   2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 13 e1   13 e2
                                 ----- + -----
                                   2       2

                                 13 e1   15 e2
                                 ----- + -----
                                   2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 11 e1   17 e2
                                 ----- + -----
                                   2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 3 e1   15 e2
                                 ---- + -----
                                  2       2

                                 3 e1   17 e2
                                 ---- + -----
                                  2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                 7 e1   19 e2
                                 ---- + -----
                                  2       2

                                  e1    15 e2
                                 ---- + -----
                                  2       2

                                 3 e1   15 e2
                                 ---- + -----
                                  2       2

                                   e1    17 e2
                                - ---- + -----
                                   2       2

                                   e1    19 e2
                                - ---- + -----
                                   2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                  e1    21 e2
                                 ---- + -----
                                  2       2

                                 3 e1   21 e2
                                 ---- + -----
                                  2       2

                                 7 e1   11 e2
                                 ---- + -----
                                  2       2

                                 9 e1   11 e2
                                 ---- + -----
                                  2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                 11 e1   13 e2
                                 ----- + -----
                                   2       2

                                 11 e1   15 e2
                                 ----- + -----
                                   2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 9 e1   17 e2
                                 ---- + -----
                                  2       2

                                 3 e1   13 e2
                                 ---- + -----
                                  2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                  e1    15 e2
                                 ---- + -----
                                  2       2

                                  e1    17 e2
                                 ---- + -----
                                  2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                 3 e1   19 e2
                                 ---- + -----
                                  2       2

                                 5 e1   19 e2
                                 ---- + -----
                                  2       2

                                  9 e1   9 e2
                                  ---- + ----
                                   2      2

                                 11 e1   9 e2
                                 ----- + ----
                                   2      2

                                 7 e1   11 e2
                                 ---- + -----
                                  2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 13 e1   11 e2
                                 ----- + -----
                                   2       2

                                 13 e1   13 e2
                                 ----- + -----
                                   2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 11 e1   15 e2
                                 ----- + -----
                                   2       2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                 7 e1   11 e2
                                 ---- + -----
                                  2       2

                                 3 e1   13 e2
                                 ---- + -----
                                  2       2

                                 3 e1   15 e2
                                 ---- + -----
                                  2       2

                                 9 e1   13 e2
                                 ---- + -----
                                  2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                 7 e1   17 e2
                                 ---- + -----
                                  2       2

                                  e1    13 e2
                                 ---- + -----
                                  2       2

                                 3 e1   13 e2
                                 ---- + -----
                                  2       2

                                   e1    15 e2
                                - ---- + -----
                                   2       2

                                   e1    17 e2
                                - ---- + -----
                                   2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                  e1    19 e2
                                 ---- + -----
                                  2       2

                                 3 e1   19 e2
                                 ---- + -----
                                  2       2

                                  7 e1   9 e2
                                  ---- + ----
                                   2      2

                                  9 e1   9 e2
                                  ---- + ----
                                   2      2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                 11 e1   11 e2
                                 ----- + -----
                                   2       2

                                 11 e1   13 e2
                                 ----- + -----
                                   2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 9 e1   15 e2
                                 ---- + -----
                                  2       2

                                 3 e1   11 e2
                                 ---- + -----
                                  2       2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                  e1    13 e2
                                 ---- + -----
                                  2       2

                                  e1    15 e2
                                 ---- + -----
                                  2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                 3 e1   17 e2
                                 ---- + -----
                                  2       2

                                 5 e1   17 e2
                                 ---- + -----
                                  2       2

                                  5 e1   9 e2
                                  ---- + ----
                                   2      2

                                  7 e1   9 e2
                                  ---- + ----
                                   2      2

                                 3 e1   11 e2
                                 ---- + -----
                                  2       2

                                 3 e1   13 e2
                                 ---- + -----
                                  2       2

                                 9 e1   11 e2
                                 ---- + -----
                                  2       2

                                 9 e1   13 e2
                                 ---- + -----
                                  2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                 7 e1   15 e2
                                 ---- + -----
                                  2       2

                                  e1    11 e2
                                 ---- + -----
                                  2       2

                                 3 e1   11 e2
                                 ---- + -----
                                  2       2

                                   e1    13 e2
                                - ---- + -----
                                   2       2

                                   e1    15 e2
                                - ---- + -----
                                   2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                  e1    17 e2
                                 ---- + -----
                                  2       2

                                 3 e1   17 e2
                                 ---- + -----
                                  2       2

                                  7 e1   7 e2
                                  ---- + ----
                                   2      2

                                  9 e1   7 e2
                                  ---- + ----
                                   2      2

                                  5 e1   9 e2
                                  ---- + ----
                                   2      2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                 11 e1   9 e2
                                 ----- + ----
                                   2      2

                                 11 e1   11 e2
                                 ----- + -----
                                   2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 9 e1   13 e2
                                 ---- + -----
                                  2       2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                  5 e1   9 e2
                                  ---- + ----
                                   2      2

                                  e1    11 e2
                                 ---- + -----
                                  2       2

                                  e1    13 e2
                                 ---- + -----
                                  2       2

                                 7 e1   11 e2
                                 ---- + -----
                                  2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                 3 e1   15 e2
                                 ---- + -----
                                  2       2

                                 5 e1   15 e2
                                 ---- + -----
                                  2       2

                                  5 e1   7 e2
                                  ---- + ----
                                   2      2

                                  7 e1   7 e2
                                  ---- + ----
                                   2      2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                 3 e1   11 e2
                                 ---- + -----
                                  2       2

                                  9 e1   9 e2
                                  ---- + ----
                                   2      2

                                 9 e1   11 e2
                                 ---- + -----
                                  2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                 7 e1   13 e2
                                 ---- + -----
                                  2       2

                                   e1    9 e2
                                  ---- + ----
                                   2      2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                   e1    11 e2
                                - ---- + -----
                                   2       2

                                   e1    13 e2
                                - ---- + -----
                                   2       2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                  e1    15 e2
                                 ---- + -----
                                  2       2

                                 3 e1   15 e2
                                 ---- + -----
                                  2       2

                                  3 e1   7 e2
                                  ---- + ----
                                   2      2

                                  5 e1   7 e2
                                  ---- + ----
                                   2      2

                                   e1    9 e2
                                  ---- + ----
                                   2      2

                                  e1    11 e2
                                 ---- + -----
                                  2       2

                                  7 e1   9 e2
                                  ---- + ----
                                   2      2

                                 7 e1   11 e2
                                 ---- + -----
                                  2       2

                                 3 e1   13 e2
                                 ---- + -----
                                  2       2

                                 5 e1   13 e2
                                 ---- + -----
                                  2       2

                                  5 e1   5 e2
                                  ---- + ----
                                   2      2

                                  7 e1   5 e2
                                  ---- + ----
                                   2      2

                                  3 e1   7 e2
                                  ---- + ----
                                   2      2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                  9 e1   7 e2
                                  ---- + ----
                                   2      2

                                  9 e1   9 e2
                                  ---- + ----
                                   2      2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                 7 e1   11 e2
                                 ---- + -----
                                  2       2

                                   e1    7 e2
                                  ---- + ----
                                   2      2

                                  3 e1   7 e2
                                  ---- + ----
                                   2      2

                                    e1    9 e2
                                 - ---- + ----
                                    2      2

                                   e1    11 e2
                                - ---- + -----
                                   2       2

                                  5 e1   9 e2
                                  ---- + ----
                                   2      2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                  e1    13 e2
                                 ---- + -----
                                  2       2

                                 3 e1   13 e2
                                 ---- + -----
                                  2       2

                                  3 e1   5 e2
                                  ---- + ----
                                   2      2

                                  5 e1   5 e2
                                  ---- + ----
                                   2      2

                                   e1    7 e2
                                  ---- + ----
                                   2      2

                                   e1    9 e2
                                  ---- + ----
                                   2      2

                                  7 e1   7 e2
                                  ---- + ----
                                   2      2

                                  7 e1   9 e2
                                  ---- + ----
                                   2      2

                                 3 e1   11 e2
                                 ---- + -----
                                  2       2

                                 5 e1   11 e2
                                 ---- + -----
                                  2       2

                                   e1    5 e2
                                  ---- + ----
                                   2      2

                                  3 e1   5 e2
                                  ---- + ----
                                   2      2

                                    e1    7 e2
                                 - ---- + ----
                                    2      2

                                    e1    9 e2
                                 - ---- + ----
                                    2      2

                                  5 e1   7 e2
                                  ---- + ----
                                   2      2

                                  5 e1   9 e2
                                  ---- + ----
                                   2      2

                                  e1    11 e2
                                 ---- + -----
                                  2       2

                                 3 e1   11 e2
                                 ---- + -----
                                  2       2

                                  3 e1   3 e2
                                  ---- + ----
                                   2      2

                                  5 e1   3 e2
                                  ---- + ----
                                   2      2

                                   e1    5 e2
                                  ---- + ----
                                   2      2

                                   e1    7 e2
                                  ---- + ----
                                   2      2

                                  7 e1   5 e2
                                  ---- + ----
                                   2      2

                                  7 e1   7 e2
                                  ---- + ----
                                   2      2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                  5 e1   9 e2
                                  ---- + ----
                                   2      2

                                   e1    3 e2
                                  ---- + ----
                                   2      2

                                  3 e1   3 e2
                                  ---- + ----
                                   2      2

                                    e1    5 e2
                                 - ---- + ----
                                    2      2

                                    e1    7 e2
                                 - ---- + ----
                                    2      2

                                  5 e1   5 e2
                                  ---- + ----
                                   2      2

                                  5 e1   7 e2
                                  ---- + ----
                                   2      2

                                   e1    9 e2
                                  ---- + ----
                                   2      2

                                  3 e1   9 e2
                                  ---- + ----
                                   2      2

                                   e1     e2
                                  ---- + ----
                                   2      2

                                  3 e1    e2
                                  ---- + ----
                                   2      2

                                    e1    3 e2
                                 - ---- + ----
                                    2      2

                                    e1    5 e2
                                 - ---- + ----
                                    2      2

                                  5 e1   3 e2
                                  ---- + ----
                                   2      2

                                  5 e1   5 e2
                                  ---- + ----
                                   2      2

                                   e1    7 e2
                                  ---- + ----
                                   2      2

                                  3 e1   7 e2
                                  ---- + ----
                                   2      2

                                      28

> read("coxeter.mpl");
coxeter and weyl 2.4v loaded.
Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2
iprod := (v1,v2) -> coeff(v1,lambda0)*c\

oeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda
oeff(v2,delta)+coeff(v1,delta)*coeff(v2\

,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,d
,lambda0)+coxeter['iprod'](subs({lambda\

0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
iprod := (v1, v2) -> coeff(v1, lambda0) coeff(v2, delta)

     + coeff(v1, delta) coeff(v2, lambda0) + coxeter['iprod'](

    subs({eps = 0, lambda0 = 0, delta = 0}, v1),

    subs({eps = 0, delta = 0, lamda0 = 0}, v2))


weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;
    weyl_vector := proc(R) return convert(weyl['weights'](R), `+`) end proc


finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,
finite_orbit:=proc(v0,R) local S,coS,EP\

S,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
    S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
    if type(v0,'set') or type(v0,'list') then orb:=op(v0)
    else orb:=coxeter['vec2fc'](v0,S)
    fi;
    if type([op(S),orb],'list'('polynom'('rational'))) then
        EPS:=0; is_mem:=(x,y,z)->member(x,y)
    else
        EPS:=`coxeter/default`['epsilon'];
        is_mem:=`coxeter/orbit/f`
    fi;
    new:=[orb];
    nparity:=[seq(-1,i=1..nops(new))];
    pres:=nparity;
    while nops(new)>0 do
        old:=new; new:=[];
        parity:=nparity;
        nparity:=[];
        for j to nops(old) do
            u:=old[j];
            for i to nops(S) do
                v:=iprod(coS[i],u);#coxeter['iprod'](coS[i],u);
                if v>EPS then
                    v:=u-v*S[i];
                    if not is_mem(v,new,EPS) then
                        new:=[op(new),v];
                        nparity:=[op(nparity),parity[j]*(-1)];
                    fi
                fi
            od
        od;
        orb:=orb,op(new);
        pres:=[op(pres),op(nparity)];
    od;
    return zip((x,y)->x+y*eps,[orb],pres);
end:

star:=proc(R)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,finite_orbit(rho,R));
end proc;
star := proc(R)
local rho;
    rho := weyl_vector(R); return map(x -> rho - x, finite_orbit(rho, R))
end proc


multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j\

,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R;
    if nargs<>3 then
        R:=args[2];
        S:=coxeter['base'](args[2])
    else
        R:=args[3];
        S:=coxeter['base'](args[3]);
        u:=coxeter['vec2fc'](args[2],S);
        v:=coxeter['root_coords'](v0-u,S);
        if not type(v,'list'('nonnegint')) then RETURN(0) fi
    fi;
    wts:=weyl['weight_sys'](v0,S,'wl');
    if nargs>3 then pr:=args[3]; prwc:=args[4] else
        pr:=coxeter['pos_roots'](S);
        coS:=coxeter['co_base'](S);
        prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
    fi;
    r0:=convert(pr,`+`)/2; mults[1]:=1;
    if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
    thestar:=star(R);
    starSign:=map(x->coeff(x,eps),thestar);
    starW:=map(x->subs(eps=0,x),thestar);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 2 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next fi;
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then mults[n] else
        convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:

> weyl['weight_mults'](wg[1]*155,wg[1],B2);
                                      78

> multiplicities(wg[1]*155,wg[1],B2);
                                      78

> time(multiplicities(wg[1]*155,wg[1],B2));
                                     2.920

> time(multiplicities(wg[1]*201,wg[1],B2));
                                     5.844

> time(multiplicities(wg[1]*401,wg[1],B2));
                                    49.935

> xs:=[1..100];
                               xs := [1 .. 100]

> xs[3];
Error, invalid subscript selector
> xs:=[seq(i,i=1..100)];
xs := [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,

    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,

    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,

    59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77,

    78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96,

    97, 98, 99, 100]

> ts:=map(x->time(multiplicities(wg[1]*(2*x+1),wg[1],B2)),xs);
ts := [0.008, 0.008, 0.012, 0.016, 0.032, 0.024, 0.096, 0.032, 0.036, 0.056,

    0.048, 0.100, 0.064, 0.116, 0.080, 0.132, 0.104, 0.148, 0.168, 0.180,

    0.192, 0.204, 0.220, 0.276, 0.252, 0.308, 0.284, 0.348, 0.360, 0.376,

    0.400, 0.464, 0.436, 0.500, 0.528, 0.544, 0.572, 0.648, 0.664, 0.652,

    0.836, 0.844, 0.860, 0.852, 0.888, 0.908, 0.992, 1.012, 1.056, 1.128,

    1.172, 1.200, 1.288, 1.332, 1.344, 1.436, 1.532, 1.516, 1.596, 1.860,

    1.748, 1.784, 1.948, 1.940, 2.028, 2.128, 2.156, 2.252, 2.300, 2.416,

    2.512, 2.568, 2.672, 2.776, 2.836, 2.936, 3.044, 3.156, 3.220, 3.388,

    3.464, 3.564, 3.624, 3.780, 3.896, 3.968, 4.136, 4.248, 4.372, 4.528,

    4.728, 4.852, 4.968, 5.476, 5.696, 5.796, 5.748, 5.680, 5.924, 6.300]

> plotsetup(X11);
> rr:=zip((x,y)->[x,y],xs,ts);
rr := [[1, 0.008], [2, 0.008], [3, 0.012], [4, 0.016], [5, 0.032], [6, 0.024],

    [7, 0.096], [8, 0.032], [9, 0.036], [10, 0.056], [11, 0.048], [12, 0.100],

    [13, 0.064], [14, 0.116], [15, 0.080], [16, 0.132], [17, 0.104],

    [18, 0.148], [19, 0.168], [20, 0.180], [21, 0.192], [22, 0.204],

    [23, 0.220], [24, 0.276], [25, 0.252], [26, 0.308], [27, 0.284],

    [28, 0.348], [29, 0.360], [30, 0.376], [31, 0.400], [32, 0.464],

    [33, 0.436], [34, 0.500], [35, 0.528], [36, 0.544], [37, 0.572],

    [38, 0.648], [39, 0.664], [40, 0.652], [41, 0.836], [42, 0.844],

    [43, 0.860], [44, 0.852], [45, 0.888], [46, 0.908], [47, 0.992],

    [48, 1.012], [49, 1.056], [50, 1.128], [51, 1.172], [52, 1.200],

    [53, 1.288], [54, 1.332], [55, 1.344], [56, 1.436], [57, 1.532],

    [58, 1.516], [59, 1.596], [60, 1.860], [61, 1.748], [62, 1.784],

    [63, 1.948], [64, 1.940], [65, 2.028], [66, 2.128], [67, 2.156],

    [68, 2.252], [69, 2.300], [70, 2.416], [71, 2.512], [72, 2.568],

    [73, 2.672], [74, 2.776], [75, 2.836], [76, 2.936], [77, 3.044],

    [78, 3.156], [79, 3.220], [80, 3.388], [81, 3.464], [82, 3.564],

    [83, 3.624], [84, 3.780], [85, 3.896], [86, 3.968], [87, 4.136],

    [88, 4.248], [89, 4.372], [90, 4.528], [91, 4.728], [92, 4.852],

    [93, 4.968], [94, 5.476], [95, 5.696], [96, 5.796], [97, 5.748],

    [98, 5.680], [99, 5.924], [100, 6.300]]

> plots[pointplot](rr);
>     weyl['weight_mults'](10*e1,e1,A1);
                                       1

>     weyl['weight_mults'](100*e1,e1,A1);
                                       1

>     weyl['weight_mults'](1000*e1,e1,A1);
Interrupted
> xs:=[seq(i,i=1..50)];
xs := [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,

    21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39,

    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50]

> ts:=map(x->time(weyl['weight_mults'](wg[1]*(2*x+1),wg[1],B2)),xs);
ts := [0.008, 0.016, 0.024, 0.044, 0.056, 0.112, 0.092, 0.160, 0.196, 0.236,

    0.320, 0.380, 0.460, 0.568, 0.648, 0.780, 0.880, 1.040, 1.200, 1.404,

    1.588, 1.752, 2.172, 2.260, 2.576, 2.792, 3.156, 3.544, 3.940, 4.320,

    4.772, 5.060, 5.724, 6.068, 6.600, 7.196, 7.708, 8.396, 9.068, 9.796,

    10.412, 10.956, 11.524, 12.484, 13.768, 14.528, 15.640, 16.505, 17.913,

    18.665]

> rr:=zip((x,y)->[x,y],xs,ts);
rr := [[1, 0.008], [2, 0.016], [3, 0.024], [4, 0.044], [5, 0.056], [6, 0.112],

    [7, 0.092], [8, 0.160], [9, 0.196], [10, 0.236], [11, 0.320], [12, 0.380],

    [13, 0.460], [14, 0.568], [15, 0.648], [16, 0.780], [17, 0.880],

    [18, 1.040], [19, 1.200], [20, 1.404], [21, 1.588], [22, 1.752],

    [23, 2.172], [24, 2.260], [25, 2.576], [26, 2.792], [27, 3.156],

    [28, 3.544], [29, 3.940], [30, 4.320], [31, 4.772], [32, 5.060],

    [33, 5.724], [34, 6.068], [35, 6.600], [36, 7.196], [37, 7.708],

    [38, 8.396], [39, 9.068], [40, 9.796], [41, 10.412], [42, 10.956],

    [43, 11.524], [44, 12.484], [45, 13.768], [46, 14.528], [47, 15.640],

    [48, 16.505], [49, 17.913], [50, 18.665]]

> plots[pointplot](rr);
> weights(C4);
                                  weights(C4)

> weyl['weights'](C4);
                [e1 + e2 + e3 + e4, e2 + e3 + e4, e3 + e4, e4]

> wg:=weyl['weights'](C4);
             wg := [e1 + e2 + e3 + e4, e2 + e3 + e4, e3 + e4, e4]

>     weyl['weight_mults'](10*wg[1],wg[1],C4);
                                     2485

> time(weyl['weight_mults'](10*wg[1],wg[1],C4));
                                    10.788

>     multiplicities(10*wg[1],wg[1],C4);
                                     2485

>time(multiplicities(10*wg[1],wg[1],C4));
                                    52.275

>
