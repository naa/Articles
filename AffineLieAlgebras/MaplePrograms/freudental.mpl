read("affine.mpl"):

test:=proc(x,R)
    local iw,pr,p,w;
    iw := map(w->subs(lambda0=0,w)+level(x)*lambda0,dominant_weights(R,level(x)));
    pr:=coxeter['pos_roots'](finite_dimensional_root_system(R));
    map(pos_r->[pos_r,map(w->[w,isolve(iprod(x,x)+2*p*iprod(x,pos_r)+p*p*iprod(pos_r,pos_r)-iprod(w,w)+2*q*iprod(w,delta))],iw)
        ],pr);
end proc;

calculateMultiplicity:=proc(weight,highest_weight,dominant_weights,positive_roots,rho,result_table)
local bare_weight,root,dom_weight,solutions,res;
    if coeff(weight,delta)>0 then return 0; fi;
    bare_weight:=subs(eps=0,weight);
    if not assigned(result_table[bare_weight]) then
        res:=0;
        for root in positive_roots do
            for dom_weight in dominant_weights do
                solutions:='solutions';
                solutions:=isolve(iprod(weight,weight)+2*p*iprod(weight,root)+p*p*iprod(root,root)-iprod(dom_weight,dom_weight)+2*q*iprod(dom_weight,delta));
                if assigned(solutions) and solutions<>NULL then
                    eq1:=select(x->rhs(x)=_Z1,solve(subs(solutions,{p>0,q>=0})))[1];
                    eq2:=select(x->lhs(x)=_Z1,solve(subs(solutions,{p>0,q>=0})))[1];
                    for i from ceil(lhs(eq1)) to floor(rhs(eq2)) do
                        concrete_solution:=subs(_Z1=i,solutions);
#                        print(concrete_solution);
                        r:=rhs(select(x->lhs(x)=q,concrete_solution)[1]);
                        k:=rhs(select(x->lhs(x)=p,concrete_solution)[1]);
                        if r<=abs(coeff(weight,delta)) and k>0 then
                            res:=res+iprod(weight+k*root,root)*calculateMultiplicity(dom_weight-r*delta,highest_weight,dominant_weights,positive_roots,rho,result_table)
                        end if:
                    end do:
                end if:
            end do:
        end do:
        result_table[bare_weight]:=res*2/(iprod(highest_weight+rho,highest_weight+rho)-iprod(weight+rho,weight+rho));
    fi;
    return result_table[bare_weight];
end proc:

stringFunction:=proc(weight,highest_weight,max_grade,R)
    rt:=table();
    rt[highest_weight]:=1;
    R0:=finite_dimensional_root_system(R);
    dominant_weights:=dominant_weights(R,coeff(weight,lambda0));
    positive_roots:=coxeter['pos_roots'](R0);
    for i from 1 to coxeter['rank'](R0) do
        positive_roots:=[op(positive_roots),delta];
    end do:
#    print(positive_roots);
    rho:=weyl_vector(R);
    return map(x->calculateMultiplicity(weight-x*delta,highest_weight,dominant_weights,positive_roots,rho,rt),[$0..max_grade]);
end proc:
