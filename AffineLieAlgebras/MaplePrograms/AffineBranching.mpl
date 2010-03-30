#!/opt/maple11/bin/maple
read("affine.mpl"):
plotsetup(X11);
get_indices:=r->map(x->x[1],[indices(r)]);
iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)
+coeff(v1,delta)*coeff(v2,lambda0)+
coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),
                 subs({lamda0=0,delta=0,eps=0},v2));
pairing:=(v1,v2) -> 2*iprod(v1,v2)/iprod(v2,v2);

projection:=proc(weights,base)
local cm;
    if nops(base)<=0 then
        return [];
    fi;
    cm:=linalg[inverse](
        array([seq([seq(iprod(r1,s),s=base)],r1=[seq(2*r/iprod(r,r),r=base)])]));
    if type(weights,'list') then
        map(x->convert(zip(`*`,
                           convert(linalg[multiply](cm,array(map(y->pairing(x,y),base))),'list'),
                           base)
                       ,`+`),weights);
    else
        convert(zip(`*`,
                    convert(linalg[multiply](cm,array(map(y->pairing(weights,y),base))),'list'),
                    base),`+`);
    end;
end;

projection1:=proc(weights,roots)
    map(x->projection(subs(eps=0,x),roots[1..-2])+delta*coeff(x,delta)+eps*coeff(x,eps),weights); end;

affine_projection:=proc(w,roots)
           if type(w,'list') or type(w,'set') then
               map(x->projection(subs(eps=0,x),roots[1..-2])+delta*coeff(x,delta)+eps*coeff(x,eps),w);
           else
               projection(subs(eps=0,w),roots[1..-2])+delta*coeff(w,delta)+eps*coeff(w,eps);
           end;
       end;

embedding_level:=proc(embedded_roots,algebra_name)
          local pa,s;
              pa:=projection(highest_root(algebra_name),embedded_roots[1..-2]);
              s:=coxeter['highest_root'](embedded_roots[1..-2]);
              return iprod(pa,pa)/iprod(s,s);
          end;

calculate_branching_coefficient:=proc(xi,fan,gamma0,ap,res,is_in_borders)
                    local xi0,a;
                        if not is_in_borders(xi) then return 0; fi;
                        xi0:=subs(eps=0,xi);
                        if assigned(res[xi0]) then
                            return res[xi0];
                        fi;
                        print(xi);

                        res[xi0]:=convert(
                            map(gam->coeff(gam,eps)*
                                calculate_branching_coefficient(xi+subs(eps=0,gam),
                                                                fan,gamma0,ap,res,
                                                                is_in_borders),fan)
                            ,`+`);
                        if assigned(ap[xi0-subs(eps=0,gamma0)]) then
                            res[xi0]:=res[xi0]+ap[xi0-subs(eps=0,gamma0)];
                        end;
                        res[xi0]:=-1/coeff(gamma0,eps)*res[xi0];
                        return res[xi0];
                    end proc:

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

external_border:=proc(points)
         local p,p1,res,ind,c1,c2;
             res:=0;
             for p1 in points do
                 p:=subs({eps=0,lambda0=0},p1);
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
            t[subs(eps=0,x)]:=1;
        od;
        return t;
    end proc:

dimension:=proc(v,pos_roots) local f,r0;
    r0:=convert(pos_roots,`+`)/2;
    f:=[seq(1+coxeter['iprod'](r,v)/coxeter['iprod'](r,r0),r=pos_roots)];
    convert(f,`*`)
end;

branching:=proc(highest_weight,subalgebra_roots,subalgebra_pos_roots,algebra_name,max_grade)
local algebra_pos_roots, algebra_simple_roots,
    rho,
    Abot_roots,Abot_rho,
    selected_points,
    ppts,
    f, inj_roots,
    Gamma, gamma0,
    k;
    algebra_simple_roots:=algebra_roots(algebra_name);
    #print(algebra_simple_roots);
    algebra_pos_roots:=positive_roots(algebra_name,max_grade);
    rho:=weyl_vector(algebra_name);
    #print(rho);
    # Next line is incorrect
    # I don't know how to fix it.
    #subalgebra_pos_roots:=map(x->-x,
    #                          coxeter['orbit'](map(x->-x,subalgebra_roots[1..-2]),
    #                                           subalgebra_roots[1..-2]));
    #print(subalgebra_pos_roots);
    Abot_roots:=select(x->affine_projection(subs(eps=0,x),subalgebra_roots)=0,algebra_pos_roots);
    Abot_rho:=1/2*convert(Abot_roots,`+`);
    #print(Abot_roots);
    anom_points:=[op(anomalous_points(highest_weight,algebra_name,max_grade))]:
    selected_points:=select(x->
                            andmap(y->iprod(x+rho-projection(x+rho,subalgebra_roots),y)>=0,
                                   Abot_roots),
                            anom_points):
    #print(nops(selected_points));

    ppts:=map(x->
              affine_projection(subs(eps=0,x),subalgebra_roots)+eps*coeff(x,eps)*
              dimension(projection(x+rho,Abot_roots)-Abot_rho,Abot_roots),
              selected_points):
    mult_table:=table();

    for r in affine_projection([op({op(algebra_pos_roots)}
                                   minus {op(Abot_roots)})],subalgebra_roots)
    do
        #        print(r);
        if assigned(mult_table[subs(eps=0,r)]) then
            mult_table[subs(eps=0,r)]:=mult_table[subs(eps=0,r)]+coeff(r,eps);
        else
            mult_table[subs(eps=0,r)]:=coeff(r,eps);
        end;
    end;


    pppr:=map(x->x+eps*mult_table[x],get_indices(mult_table));
    print(pppr);
    f:=fan(pppr,inj_roots(subalgebra_pos_roots));
    f1:=select(x->coeff(x,delta)<=max_grade,f);

    gamma0:=f1[1];
    subalgebra_rho:=subs(delta=lambda0,convert(select(x->coeff(x,delta)=0,subalgebra_roots),`+`)/2);

    for r in f1 do
        if iprod(r,subalgebra_rho)<iprod(gamma0,subalgebra_rho) and coeff(r,delta)<=coeff(gamma0,delta) then
            gamma0:=r;
        end;
    end;
    Gamma:=select(y->subs(eps=0,y)<>0,map(x->x-subs(eps=0,gamma0),f1));

    borders:=external_border([op(ppts),1/2*delta,-max_grade*delta]);
    print(borders);
    i_i_b:=rcurry(is_in_borders,borders);

    sing_table:=table();
    for v in ppts do
        if assigned(sing_table[subs(eps=0,v)]) then
            sing_table[subs(eps=0,v)]:=sing_table[subs(eps=0,v)]+coeff(v,eps);
        else

            sing_table[subs(eps=0,v)]:=coeff(v,eps);
        fi;
    end;
    print(gamma0);
#    mu0:=affine_projection(-16*(e1-e3)/2-max_grade*delta,subalgebra_roots);
    mu0:=affine_projection(subs(lambda0=0,highest_weight)-max_grade*delta,subalgebra_roots);
    t:=table();
#    calculate_branching_coefficient(-mu0,Gamma,gamma0,sing_table,t,i_i_b);

    pppp:=dominant_weights(subalgebra_roots,embedding_level(subalgebra_roots,algebra_name));

    print(pppp);
    map(x->calculate_branching_coefficient(subs(lambda0=0,x)-max_grade*delta,Gamma,gamma0,sing_table,t,i_i_b),
        map(y->subs(eps=0,y),finite_orbit(pppp,finite_dimensional_root_system(subalgebra_roots))));

    return [t,sing_table,Gamma,gamma0];
end;

projection_with_mults:=proc(weights,subalgebra_roots)
                local r, mult_table;
                    mult_table:=table();
                    for r in affine_projection(weights,subalgebra_roots)
                    do
                        if assigned(mult_table[subs(eps=0,r)]) then
                            mult_table[subs(eps=0,r)]:=mult_table[subs(eps=0,r)]+coeff(r,eps);
                        else
                            mult_table[subs(eps=0,r)]:=coeff(r,eps);
                        end;
                    end;
                    map(x->x+eps*mult_table[x],get_indices(mult_table));
                end proc:
