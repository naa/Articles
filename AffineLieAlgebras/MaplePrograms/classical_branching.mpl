read("coxeter.mpl");
iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));
pairing:=(v1,v2) -> 2*iprod(v1,v2)/iprod(v2,v2);
get_indices:=r->map(x->x[1],[indices(r)]);

projection:=proc(wg,real_roots)
    local v,t,res,weights,icm;
    icm:=linalg[inverse](coxeter['cartan_matrix']([e3-e4,e4]));
    res:=table();
    weights:=wg;
    for v in weights do
        t:=convert(map(x->x*pairing(v,x)/2,real_roots),`+`);
        if assigned(res[t]) then
            res[t]:=res[t]+coeff(v,eps);
        else
            res[t]:=coeff(v,eps);
        fi;
    od;
    return map(x->x+res[x]*eps,get_indices(res));
end proc:


positive_roots:=proc(R)
             map(x->x+eps,coxeter['pos_roots'](R));
         end;

mult_inj:=proc(r,roots)
         if assigned(roots[subs(eps=0,r)]) then
             return roots[subs(eps=0,r)];
         else
             return 0;
         fi;
     end proc:

fan_table:=proc(pos_roots,inj_roots)
local t,n,t2,res,p,j;
    n:=coeff(pos_roots[1],eps)-mult_inj(pos_roots[1],inj_roots);
    res:=table();
    p:=subs(eps=0,pos_roots[1]);
    if nops(pos_roots)=1 then
        for j from 1 to n do
                res[j*p]:=(-1)^(j) * binomial(n,j);
        od;
        return res;
    else
        t:=fan_table(pos_roots[2..-1],inj_roots);
        t2:=fan_table([pos_roots[1]],inj_roots);
        for x in get_indices(t) do
            for j from 1 to n do
                if not(assigned(res[p*j+x])) then res[p*j+x]:=0; fi;
                res[p*j+x]:=res[p*j+x]+t[x]*(-1)^(j) * binomial(n,j);
            od;
            if not(assigned(res[x])) then res[x]:=0; fi;
            res[x]:=res[x]+t[x];
        od;
        for x in get_indices(t2) do
            if not(assigned(res[x])) then res[x]:=0; fi;
            res[x]:=res[x]+t2[x];
        od;
        return res;
    fi;
end proc:

fan:=proc(pos_roots,inj_roots)
    local t;
    t:=fan_table(pos_roots,inj_roots);
    if assigned(t[0]) then t[0]:=t[0]+1;
    else
        t[0]:=1;
    end;
    select(x->coeff(x,eps)<>0,map(x->x+t[x]*eps,get_indices(t)));
end proc:

inj_roots:=proc(roots)
        local t;
        t:=table();
        for x in roots do
            t[subs(eps=0,x)]:=coeff(x,eps);
        od;
        return t;
    end proc:

calculate_branching_coefficient:=proc(xi,fan,gamma0,ap,res,is_in_borders)
                    local xi0,a;
                        if not is_in_borders(xi) then return 0; fi;
                        xi0:=subs(eps=0,xi);
                        if assigned(res[xi0]) then
                            return res[xi0];
                        fi;
                        res[xi0]:=convert(
                            map(x->coeff(x,eps)*calculate_branching_coefficient(xi+x-gamma0,fan,gamma0,ap,res,is_in_borders),
                                select(x->subs(eps=0,x)<>subs(eps=0,gamma0),fan)),`+`);
                        a:= select(x->subs(eps=0,x)=subs(eps=0,xi-gamma0),ap);
                        if nops(a)>0 then
                            res[xi0]:=res[xi0]-convert(map(x->coeff(x,eps),a),`+`);
                        fi;
                        res[xi0]:=-1/coeff(gamma0,eps)*res[xi0];
                        return res[xi0];
                    end proc:

#inject:=(sal,al,injection)->subs({(alpha[i]=al[i])$i=1..nops(al)-1},subs(solve(injection,{e||(i$i=1..(nops(sal)-1))}),sal));

inject:=(weights,subalgebra,injection)->map(x->
                                            convert(
                                                zip((y,z)->y*z,
                                                    coxeter['root_coords'](subs(eps=0,x),subalgebra),
                                                    injection),
                                                `+`)+coeff(x,eps)*eps,weights);

plot_projection:=proc(tt,real_roots)
         plots[textplot](map(x->[pairing(x,real_roots[1])/2,pairing(x,real_roots[2])/2,coeff(x,eps)],[op(tt)]));
     end proc:

get_fan:=proc(algebra,subalgebra,injection)
    local i_i_b,pr,spr,al,sal,sal_inj,t,sing,wv,i,sl,mi,mv;
        al:=coxeter['base'](algebra);
        sal:=coxeter['base'](subalgebra);

        sal_inj:=inject(sal,subalgebra,injection);
        pr:=projection(positive_roots(algebra),sal_inj);
        pr:=select(x->subs(eps=0,x)<>0,pr);

        spr:=inject(positive_roots(subalgebra),subalgebra,injection);
        f:=[eps,op(fan(pr,inj_roots(spr)))];
        return f;
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

anomalous_points:=proc(weight,R)
          local rho;
              rho:=`weyl/rho`(R);
              return map(x->subs(eps=0,x)-rho-coeff(x,eps)*eps,finite_orbit(weight+rho,R));
          end proc;

external_border:=proc(points)
             local p,p1,res,ind,c1,c2;
             res:=0;
             for p1 in points do
                 p:=subs({delta=0,eps=0,lambda0=0},p1);
                 for ind in indets(p) do
                     c1:=coeff(p,ind);
                     c2:=coeff(res,ind);
                     if c2=0 then
                         res:=res+[c1,c1]*ind;
                     elif c1<c2[1] then
                         res:=subs(ind=0,res)+[c1,c2[2]]*ind;
                     elif c1>c2[2] then
                         res:=subs(ind=0,res)+[c2[1],c1]*ind;
                     fi;
                 od;
             od;
             return res;
         end proc;

is_in_borders:=proc(p,borders)
          local ind,c1,c2;
          for ind in indets(borders) do
              c1:=coeff(p,ind);
              c2:=coeff(borders,ind);
              if c1<c2[1] or c1>c2[2] then
                  return false;
              fi;
          od;
          return true;
      end proc;

branching:=proc(weight,algebra,subalgebra,injection)
local i_i_b,pr,spr,al,sal,sal_inj,t,sing,wv,i,sl,mi,mv;
    al:=coxeter['base'](algebra);
    sal:=coxeter['base'](subalgebra);

    sal_inj:=inject(sal,subalgebra,injection);
    pr:=projection(positive_roots(algebra),sal_inj);
    pr:=select(x->subs(eps=0,x)<>0,pr);

    spr:=inject(positive_roots(subalgebra),subalgebra,injection);
    f:=[eps,op(fan(pr,inj_roots(spr)))];
    print(f);
    wv:=op(inject([`weyl/rho`(subalgebra)],subalgebra,injection));
    sing:=projection(anomalous_points(weight,algebra),sal_inj[1..-1]):
    print(sing);
    borders:=external_border(sing);
    print(borders);
    i_i_b:=rcurry(is_in_borders,borders);

    t:=table();
    t[op(projection([weight],sal_inj))]:=1;


    sl:=select(x->coeff(x,delta)=0,f);
    mi:=1;
    mv:=iprod(sl[1],wv);
    for i from 2 to nops(sl) do
        if iprod(sl[i],wv)<mv then
            mv:=iprod(sl[i],wv);
            mi:=i;
        fi;
    od;
    calculate_branching_coefficient(0,f,sl[mi],sing,t,i_i_b);

    return t;
end proc:

branching_rules:=proc(hweight,algebra,subalgebra,injection,max_grade)
              local t,dom_weights,dw,res,inds,al,dw1;
              al:=algebra_roots(algebra);
              t:=branching(hweight,algebra,subalgebra,injection,max_grade);
              dom_weights:=dominant_weights(subalgebra,level(hweight));
              res:=table();
              inds:=get_indices(t);
              for dw in dom_weights do
                  dw1:=inject([dw],subalgebra,injection);
                  res[dw]:=sort(
                      convert(map(x->t[x]*q^(-coeff(x,delta)),
                                  select(y->subs({eps=0,delta=0},y)=dw1,inds))
                              ,`+`),
                      q,ascending);
              od;
              return res;
          end proc;

maximal_regular_subalgebras:=proc(R)
                local cm,cm1,R1,i,j,k,result,perm,names,alg,al_roots,sub_roots,r,es,ls;
                    R1:=finite_dimensional_root_system(R);
                    cm:=matrix(coxeter['cartan_matrix'](R1));
#                    cm:=convert(cartan_matrix(R),'matrix');
                    result:=table();
                    al_roots:=coxeter['base'](R1);
                    al_roots:=algebra_roots(R);
                    for i from 1 to linalg[rowdim](cm) do
                        cm1:=linalg[delrows](linalg[delcols](cm,i..i),i..i);
                        names:=coxeter['names_of'](cm1,'perm');
                        for alg in names do
                            sub_roots:=[];
                            for j from 1 to coxeter['rank'](alg) do
                                sub_roots:=[op(sub_roots),al_roots[perm[1]]];
                                perm:=perm[2..-1];
                            end do:
                            ls:=labels(alg);
                            sub_roots:=[op(sub_roots),delta-add(sub_roots[k]*ls[k],k=1..nops(sub_roots))];
                            if assigned(result[names][alg]) then
                                result[names][alg]:=[op(result[names][alg]),sub_roots];
                            else
                                result[names][alg]:=[sub_roots];
                            fi;
                        end do:
                    end do;
                    return result;
                end proc:


character:=proc(mu,R)
    local anom_points, pr,rh,denominator,numerator;
    rh:=weyl['rho'](R);
    anom_points:=map(x->x-rh,finite_orbit(mu+rh,R));
    pr:=coxeter['pos_roots'](R);
    denominator:=convert(map(x->(1-exp(-x)),pr),`*`);
    numerator:=convert(map(x->coeff(x,eps)*exp(subs(eps=0,x)),anom_points),`+`):
    return numerator/denominator;
end;
