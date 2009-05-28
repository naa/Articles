read("affine.mpl"):

projection:=proc(wg,real_roots,img_root)
    local v,t,res,imaginary_root,weights;
    res:=table();
    imaginary_root:=img_root;
    weights:=wg;
    if grade(imaginary_root)=0 then
        weights:=subs(delta=0,weights);
        imaginary_root:=imaginary_root+delta;
    fi;
    for v in weights do
        t:=convert(map(x->x*pairing(v,x)/2,real_roots),`+`)+delta*grade(v)/grade(imaginary_root)+level(v)*lambda0;
        if assigned(res[t]) then
            res[t]:=res[t]+coeff(v,eps);
        else
            res[t]:=coeff(v,eps);
        fi;
    od;
    return map(x->x+res[x]*eps,get_indices(res));
end proc:


finite_positive_roots:=proc(R) local S;
                    S:=finite_roots(R);
                    map(x->-x,coxeter['orbit'](map(x->-x,S),S))
                end:


tA2_positive_roots:=n->[e1+i*delta+eps$i=0..n,2*e1+(2*i+1)*delta+eps$i=0..n/2,-2*e1+(2*i+1)*delta+eps$i=0..n/2,-e1+i*delta+eps$i=1..n,i*delta+eps$i=1..n];

positive_roots:=proc(R,max_grade)
         local cpr,l,p,r;
             if R=tA2 then
                 return tA2_positive_roots(max_grade);
             fi;
             if is_twisted(R) then
                 error("Twisted algebras are unsupported yet");
             fi;
             cpr:=finite_positive_roots(R);
             r:=coxeter['rank'](finite_dimensional_root_system(R));
             [op(map(x->x+l*delta+eps$l=0..max_grade,cpr)),op(map(x->-x+p*delta+eps$p=1..max_grade,cpr)),p*delta+r*eps$p=1..max_grade];
         end proc:

mult_inj:=proc(r,roots)
         if assigned(roots[subs(eps=0,r)]) then
             return roots[subs(eps=0,r)];
         else
             return 0;
         fi;
     end proc:

fan_table:=proc(pos_roots,inj_roots,max_grade)
local t,n,t2,res,p,j;
    n:=coeff(pos_roots[1],eps)-mult_inj(pos_roots[1],inj_roots);
    res:=table();
    p:=subs(eps=0,pos_roots[1]);
    if nops(pos_roots)=1 then
        for j from 1 to n do
            if coeff(j*p,delta)<=max_grade then
                res[j*p]:=(-1)^(j) * binomial(n,j);
            fi;
        od;
        return res;
    else
        t:=fan_table(pos_roots[2..-1],inj_roots,max_grade);
        t2:=fan_table([pos_roots[1]],inj_roots,max_grade);
        for x in get_indices(t) do
            for j from 1 to n do
                if coeff(p*j+x,delta)<=max_grade then
                    if not(assigned(res[p*j+x])) then res[p*j+x]:=0; fi;
                    res[p*j+x]:=res[p*j+x]+t[x]*(-1)^(j) * binomial(n,j);
                fi;
            od;
            if coeff(x,delta)<=max_grade then
                if not(assigned(res[x])) then res[x]:=0; fi;
                res[x]:=res[x]+t[x];
            fi;
        od;
        for x in get_indices(t2) do
            if coeff(x,delta)<=max_grade then
                if not(assigned(res[x])) then res[x]:=0; fi;
                res[x]:=res[x]+t2[x];
            fi;
        od;
        return res;
    fi;
end proc:

fan:=proc(pos_roots,inj_roots,max_grade)
    local t;
    t:=fan_table(pos_roots,inj_roots,max_grade);
    map(x->x+t[x]*eps,get_indices(t));
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
                        if coeff(xi,delta)>0 or not is_in_borders(xi) then return 0; fi;
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

inject:=(sal,al,injection)->subs({(alpha[i]=al[i])$i=1..nops(al)-1},subs(solve(injection,{e||(i$i=1..(nops(sal)-1))}),sal));

plot_projection:=proc(tt,real_roots,imaginary_root)
         plots[textplot](map(x->[coeff(x,delta)/coeff(imaginary_root,delta),pairing(x,real_roots[1])/2,coeff(x,eps)],[op(tt)]));
     end proc:

get_fan:=proc(weight,algebra,subalgebra,injection,max_grade)
    local i_i_b,pr,spr,al,sal,sal_inj,t,sing,wv,i,sl,mi,mv;
        al:=algebra_roots(algebra);
        sal:=algebra_roots(subalgebra);

        sal_inj:=inject(sal,al,injection);
        pr:=projection(positive_roots(algebra,max_grade),sal_inj[1..-2],delta);
        pr:=select(x->subs(eps=0,x)<>0,pr);

        spr:=inject(positive_roots(subalgebra,max_grade),al,injection);
        f:=[eps,op(fan(pr,inj_roots(spr),max_grade))];
        return f;
    end proc;

branching:=proc(weight,algebra,subalgebra,injection,max_grade)
local i_i_b,pr,spr,al,sal,sal_inj,t,sing,wv,i,sl,mi,mv;
    al:=algebra_roots(algebra);
    sal:=algebra_roots(subalgebra);

    sal_inj:=inject(sal,al,injection);
    pr:=projection(positive_roots(algebra,max_grade),sal_inj[1..-2],delta);
    pr:=select(x->subs(eps=0,x)<>0,pr);

    spr:=inject(positive_roots(subalgebra,max_grade),al,injection);
    f:=[eps,op(fan(pr,inj_roots(spr),max_grade))];

    wv:=inject(weyl_vector(subalgebra),al,injection);
    sing:=projection(anomalous_points(weight,algebra,max_grade),sal_inj[1..-2],delta):
    borders:=external_border(sing);
    i_i_b:=rcurry(is_in_borders,borders);

    t:=table();
    t[weight]:=1;


    sl:=select(x->coeff(x,delta)=0,f);
    mi:=1;
    mv:=iprod(sl[1],wv);
    for i from 2 to nops(sl) do
        if iprod(sl[i],wv)<mv then
            mv:=iprod(sl[i],wv);
            mi:=i;
        fi;
    od;
#    print(sl[mi]);
    calculate_branching_coefficient(-max_grade*delta+weight,f,sl[mi],sing,t,i_i_b);

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
                  dw1:=inject(dw,al,injection);
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
