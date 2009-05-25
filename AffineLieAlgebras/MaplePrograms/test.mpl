read("coxeter.mpl");

iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));

weyl_vector:=proc(R)
         return convert(weyl['weights'](R),`+`);
     end proc;

finite_orbit:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v,parity,u1,j,nparity,pres;
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
    local rho,res,cr0;
    rho:=weyl_vector(R);
    res:=map(x->rho-x,finite_orbit(rho,R));
    res:=map(x->[x,iprod(x,rho)],res);
    res:=sort(res,(x,y)->evalb(x[2]<y[2]));
    return map(x->x[1],res);
end proc;

multiplicities:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0,thestar,starW,starSign,R,MMatrix;
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
    MMatrix:=Matrix(n,n);

    for i from 2 to n do
        v:=wts[i];
        m:=0;
        for k from 2 to nops(starW) do
            u:=v+starW[k];
            u0:=coxeter['vec2fc'](u,S);
            if not member(u0,wts,'j') then next; fi;
            MMatrix[i,j]:=MMatrix[i,j]+starSign[k];
            m:=m-starSign[k]*mults[j];
        od;
        mults[i]:=m;
    od;
    if nargs=3 then [mults[n],MMatrix] else
        [convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`),MMatrix] fi
end:
