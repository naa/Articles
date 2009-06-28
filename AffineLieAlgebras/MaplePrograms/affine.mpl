#AffineLieAlgebra:=module()
#export iprod,roots,weights,labels,co_labels,weyl_vector,highest_root,star0;
#local t_b_action,weyl;
#option package;
#    use LinearAlgebra, ListTools in

read("coxeter.mpl");

# Inner product for root space


# scalar product in weight space
iprod := (v1,v2) -> coeff(v1,lambda0)*coeff(v2,delta)+coeff(v1,delta)*coeff(v2,lambda0)+coxeter['iprod'](subs({lambda0=0,delta=0,eps=0},v1),subs({lamda0=0,delta=0,eps=0},v2));

pairing:=(v1,v2) -> 2*iprod(v1,v2)/iprod(v2,v2);

# checks, whether given affine algebra is twisted. So to specify
# algebra D_3^{(2)} we should use the notation tD3
is_twisted := proc(R)
       if substring(R,1..1)='t' then
           return true;
       else return false;
       end if;
   end proc;

# Gets root system of finite dimensional subalgebra
# For example, it will be root system of A1 if specified algebra is A1
finite_dimensional_root_system := proc(R)
                        local r,S,L;
                            if type(R,'indexed') then
                                L:=op(0,R)
                            else
                                L:=substring(R,1..1)
                            fi;
                            if is_twisted(R) then
                                L:=substring(R,2..2);
                                r:=parse(substring(R,3..-1));
                                if L='A' then
                                    return C||(floor(r/2));
                                elif L='D' then
                                    return B||(r-1);
                                end if
                            else
                                return R;
                            end if;
                        end proc;

# Returns arrayof finite roots
finite_roots:=proc(R1)
           local R;
           R:=finite_dimensional_root_system(R1);
           if is_twisted(R1) then
               return coxeter['co_base'](R);
           else
               return coxeter['base'](R);
           fi;
       end proc;

# All roots of affine algebra: finite roots + \alpha_0. \alpha_0 is
# the last element of the array of roots
algebra_roots := proc(R1)
        local al,ls,R;
            al:=finite_roots(R1);
            ls:=labels(R1);
            al:=[op(al),delta-add(al[i]*ls[i],i=1..nops(al))];
            return al;
        end proc;

# weights of algebra, finite first, then lambda_0
weights:=proc(R1)
local wg,al,i,R,lb,clb;
    R:=finite_dimensional_root_system(R1);
    al:=algebra_roots(R1);
    wg:=weyl['weights'](R);
    lb:=labels(R1);
    for i to nops(wg) do
        wg[i]:=wg[i]+lambda0;
    end do;
    wg:=[op(wg),lambda0];
    return wg;
end proc;

dominant_weights := proc(R,level)
           local res,wg,clb,w,addon,term,k,wi,n;
               clb:=co_labels(R);
               wg:=weights(R);
               res:={[0 $ i=1..nops(wg)]};
               for w from 1 to nops(wg)-1 do
                   addon:={};
                   k:=0;
                   term:=false;
                   while not term do
                       k:=k+1;
                       term:=true;
                       for wi in res do
                           nwi:=wi;
                           nwi[w]:=nwi[w]+k;
                           n:=add(nwi[i]*clb[i],i=1..nops(nwi)-1);
                           if n <= level then
                               addon:={op(addon),nwi};
#                               print(addon);
                               term:=false;
                           fi;
                       od;
                   od;
                   res:={op(res),op(addon)};
               od;
               return map(x->subs(lambda0=0,add(x[i]*wg[i],i=1..nops(x)))+level*lambda0,res);
           end proc;

level:=x->coeff(x,lambda0);

grade:=x->coeff(x,delta);

# if vector is in main chamber of algebra
is_in_main_chamber:=proc(vector,R)
               local al,a;
               al:=algebra_roots(R);
               for a in al do
                   if iprod(vector,a)<0 then
                       return false;
                   fi;
               od;
               return true;
end proc;

# labels
labels:=proc(R)
local r,S,L;

    if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
    if is_twisted(R) then
        L:=substring(R,2..2);
        r:=parse(substring(R,3..-1));
        if L='A' then
            if r = 2 then
                S:=[2,1];
            elif r mod 2 = 0 then
                S:=[seq(2,j=1..(r/2-1)),1,2];
            else
                S:=[1,seq(2,j=1..(floor(r/2)-2)),1,1];
            end if;
        elif L='D' then
            S:=[seq(1,j=1..r)];
        else
            error "Unsupported";
        end if;
    else
        r:=coxeter['rank'](R);
        if L='A' then
            S:=[seq(1,j=1..r+1)];
        elif L='B' then
            S:=[1,seq(2,j=1..r-1),1];
        elif L='C' then
            S:=[seq(2,j=1..r-1),1,1];
        elif L='D' then
            S:=[1,seq(2,j=1..r-3),1,1,1];
        elif L='E' then
            if r=6 then
                S:=[1,2,3,2,2,1,1];
            elif r=7 then
                S:=[1,2,3,4,2,3,2,1];
            elif r=8 then
                S:=[2,3,4,5,6,3,4,2,1];
            fi;
        elif L='F' and r=4 then
            S:=[2,3,4,2,1];
        elif L='G' then
            S:=[2,3,1];
        fi;
    end if;
    return S;
end proc;

algebra_co_roots:=proc(R)
               local r;
               [seq(2*r/iprod(r,r),r=algebra_roots(R))]
           end proc;
#colabels
co_labels:=proc(R1)
   local l,r;
       r:=algebra_roots(R1);
       l:=labels(R1);
       return [seq(l[i]*iprod(r[i],r[i])/2,i=1..nops(r))];
   end proc;

cartan_matrix2:=proc(R)
           local S,coS,r,s;
           S:=algebra_roots(R);
           coS:=algebra_co_roots(R);
           array([seq([seq(iprod(r,s),s=S)],r=coS)]);
       end proc;

cartan_matrix:=proc(R)
           local m,al,ls,cls,i,j;
           al:=algebra_roots(R);
           ls:=labels(R);
           cls:=co_labels(R);
           m:=Matrix(nops(al),nops(al));
           for i from 1 to nops(al) do
               for j from 1 to nops(al) do
                   m[j,i]:=ls[j]/cls[j]*iprod(al[i],al[j]);
               od;
           od;
           return m;
       end proc;

#highest root? Unused
highest_root:=proc(R)
        local ls,rt;
            rt:=algebra_roots(R);
            ls:=labels(R);
            return add(rt[i]*ls[i],i=1..nops(rt)-1);
        end proc;

# weyl vector - sum of all weights
weyl_vector:=proc(R)
         return convert(weights(R),`+`);
     end proc;


t_b_action := proc(b,lambda)
#        if b=0 then return lambda; fi;
        return lambda+iprod(lambda,delta)*b-(iprod(b,b)/2*iprod(lambda,delta)+iprod(lambda,b))*delta;
    end proc;

translations:=proc(weights,R,max_grade)
    if type(weights,sequential) then
        return {op(map(x->op(translations_correct(x,R,max_grade)),weights))};
    else
        return {op(translations_correct(weights,R,max_grade))};
    fi;
end proc:

translations_correct:=proc(w,R,max_grade)
             local al,A,K,C,D,L,low,high,v;
#                 print(w);
                 al:=coxeter['base'](finite_dimensional_root_system(R));
                 A:=Matrix(map(x->map(y->1/2*iprod(x,y)*iprod(w,delta),al),al));
#                 print(A);
                 K:=Vector(map(x->iprod(w,x),al));
#                 print(K);
                 C:=max_grade+coeff(w,delta);
#                 print(C);
                 L,D:=LinearAlgebra[Eigenvectors](A);
#                 print(L);
#                 print(D);
#                 print(D.A.LinearAlgebra[MatrixInverse](D));
                 dim:=LinearAlgebra[Dimension](L);
                 for i from 1 to dim do
                     for j from 1 to dim do
                         v:=floor(evalf(-sqrt((C+1/4*(LinearAlgebra[Transpose](K).A.K))/L[i])*D[i,j]-
                                        (LinearAlgebra[MatrixInverse](A).K)[j]/2));
                         if not(assigned(low[j])) or v<low[j] then
                             low[j]:=v;
                         fi;
                         v:=ceil(evalf(sqrt((C+1/4*(LinearAlgebra[Transpose](K).A.K))/L[i])*D[i,j]-
                                        (LinearAlgebra[MatrixInverse](A).K)[j]/2));
                         if not(assigned(high[j])) or v>high[j] then
                             high[j]:=v;
                         fi;
                     od;
                 od;
#                 print(low);
#                 print(high);
                 res:= translate(w,coxeter['base'](finite_dimensional_root_system(R)),convert(low,list),convert(high,list),1,0,max_grade);
#                 print(nops(res));
                 return res;
             end proc:

translate:=proc(weight,roots,low,high,j,vec,max_grade)
local k,v,res;
    if j>nops(low) then
        v:= t_b_action(vec,weight);
        if abs(coeff(v,delta))<= max_grade then
            return [v];
        else return [];
        fi;
    fi:
    v:= t_b_action(vec,weight);

    if abs(coeff(v,delta))<= max_grade then
        res:=[v];
    else
        res:=[];
    fi;
    for k from low[j] to high[j] do
        addon:=translate(weight,roots,low,high,j+1,vec+k*roots[j],max_grade);
#        print(addon);
        res:=[op(res),op(addon)];
    od:
#    print(res);
    return res;
end proc:

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

orbit:=proc(w,R,size)
    return translations(finite_orbit(w,finite_dimensional_root_system(R)),R,size);
end proc;

external_contour:=proc(weight,R,size)
             return orbit(weight,R,size);
         end proc;

anomalous_points:=proc(weight,R,size)
          local rho;
              rho:=weyl_vector(R);
              return map(x->subs(eps=0,x)-rho-coeff(x,eps)*eps,translations(finite_orbit(weight+rho,finite_dimensional_root_system(R)),R,size));
          end proc;

star:=proc(R,max_grade)
    local rho;
    rho:=weyl_vector(R);
    return map(x->rho-x,orbit(rho,R,max_grade));
end proc;

multiplicities := proc (weight,R,md)
    local anomPoints,theStar,highestPoints,p,resTable,borders,maximum,max_delta;
    max_delta:=md;
    anomPoints:=anomalous_points(weight,R,max_delta);
    theStar:=select(x->subs(eps=0,x)<>0,star(R,max_delta));
#    highestPoints:= select(p->abs(coeff(p,delta))=max_delta,anomPoints);
#    borders:=external_border(highestPoints);

#    borders:=external_border(anomPoints);
    resTable:=table();
    maximum:=min(op(map(x->coeff(x,delta),anomPoints)));
    borders:=coeff(t_b_action(-finite_part(weight)/coeff(weight,lambda0),weight),delta);
#    calc_mult(weight-max_delta*delta,theStar,resTable,anomPoints,[weight]);
#    calc_mult(weight-max_delta*delta,theStar,resTable,anomPoints,[op(finite_orbit(weight,finite_dimensional_root_system(R)))],R,max_delta);
    map(x->calc_mult(x,theStar,resTable,anomPoints,borders,R,max_delta),
        map(x->subs(delta=0,x)-max_delta*delta,orbit(weight,R,max_delta)));
    return resTable;
end proc;

multiplicity:=proc(weight,highest_weight,R)
    local max_grade,r;
    max_grade:=-coeff(weight,delta);
    r:=multiplicities(highest_weight,R,max_grade);
    return r[weight];
end proc;

fusion_coefficient:=proc(mu,lambda,nu,R)
           convert(map(x->coeff(x,eps)*multiplicity(subs(eps=0,x-subs(lambda0=0,lambda)),subs(eps=0,mu),R),orbit(nu,R,3)),`+`);
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

add_star:= (w,star)-> map(x->x+w,star);

detw:=(w,anomPoints)->
coeff(
    apply(t->if nops(t)=0 then 0; else t[1]; fi,
          select(x->subs(eps=0,x)=subs(eps=0,w),anomPoints))
    ,eps);

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

is_in_borders2:=proc(p,ws)
          local kn,hweight;
          for hweight in ws do
              kn:=coeff(hweight,lambda0);
              if coeff(p,delta)<=-iprod(hweight,p-hweight)/kn-iprod(p-hweight,p-hweight)/2/kn then
                  return true;
              fi;
          od;
          return false;
      end proc;

is_in_borders3:=proc(p,hweight,star)
          local kn;
          kn:=coeff(hweight,lambda0);
          for st in star do
              if coeff(p+st,delta)<=-iprod(hweight,p+st-hweight)/kn-iprod(p+st-hweight,p+st-hweight)/2/kn then
                  return true;
              fi;
          od;
          return false;
      end proc;

is_in_borders4:=proc(p,borders,R)
          return is_in_borders(p,borders) and is_in_main_chamber(p,R);
      end proc;

vector_length:=v->foldl((x,y)->x+y^2,0,coeffs(v));

is_in_borders5:=(p,hw)->
    ormap(
    hw->evalb(
        vector_length(finite_part(p))<vector_length(finite_part(hw)) or
        coeff(p,delta)<=coeff(t_b_action(subs({delta=0,lambda=0},p-hw),hw)+hw,delta)),hw);
#        coeff(p,delta)<=coeff(t_b_action(subs({delta=0,lambda=0},p-hw),lambda0*coeff(hw,lambda0))+hw,delta)),hw);

is_in_borders6:=(p,hw)->
    ormap(
    hw->evalb(
        vector_length(finite_part(p))<vector_length(finite_part(hw)) or
        coeff(p,delta)<=coeff(t_b_action(finite_part(p)-finite_part(hw),hw),delta)),
        hw);


#        coeff(p,delta)<=coeff(t_b_action(subs({delta=0,lambda=0},p-hw),lambda0*coeff(hw,lambda0))+hw,delta)),hw);

is_in_borders7:=(p,hw,R)->
evalb(iprod(p,lambda0+foldl((x,y)->x+y,0,op(zip((x,y)->x*y,co_labels(R),finite_roots(R)))))<=
      iprod(hw,lambda0+foldl((x,y)->x+y,0,op(zip((x,y)->x*y,co_labels(R),finite_roots(R))))));


is_in_borders8:=(p,borders)->
evalb(coeff(t_b_action(-finite_part(p)/coeff(p,lambda0),p),delta)<=borders);

calc_mult:=proc (w0,star,tab,anomPoints,borders,R,size)
         local res,w;
         w:=subs(eps=0,w0);
         if coeff(w,delta)>0 or not is_in_borders8(w0,borders) then return 0; fi;
         if assigned(tab[w]) then return tab[w];
         else
             res:=detw(w,anomPoints);
             res:=res-add(coeff(x,eps)*calc_mult(w+x,star,tab,anomPoints,borders,R,size),x in star);

#             res:=add(coeff(x,eps)*calc_mult(w+x,star,tab,anomPoints,borders,R,size),x in star);
#             res:=res+detw(w,anomPoints);
             tab[w]:=res;
             return res;
         fi;
     end proc;

get_indices:=r->map(x->x[1],[indices(r)]);


draw2d := proc(r)
    plots[textplot](map(x->[coeff(x,delta),coeff(x,e1),-r[x]],get_indices(r)));
end proc;

draw3d := proc(r)
    plots[textplot3d](map(x->[coeff(x,delta),coeff(x,e1),coeff(x,e2),r[x]],get_indices(r)));
end proc;

ap_draw2d := proc(ap)
    plots[textplot](map(x->[coeff(x,delta),coeff(x,e1),coeff(x,eps)],ap));
end proc;

sub_draw2d := proc(points,r)
    plots[textplot](map(x->[coeff(x,delta),coeff(x,e1),r[x]],points));
end proc;

# get weights from result table
extract_weights:=get_indices;

# classical part of weight vector (put delta=0, lambda0=0)
finite_part:=x->subs({lambda0=0,delta=0,eps=0},x);

# This function selects string function coefficients from all the multiplicities
filter_string_function:=(weight,resTable)->
sort(
    map(x->resTable[x],
        select(x->finite_part(x)=finite_part(weight),
               extract_weights(resTable))));

fold_weight:=proc(weight0,R,mg)
     local fp,a,al,wg,ex,weights,weight,p,p1;
         weights:=finite_orbit(subs(eps=0,weight0),finite_dimensional_root_system(R));
         al:=algebra_roots(R);

         for weight in weights do
             wg:=weight;
             ex:=false;
             while (not ex) and (not is_in_main_chamber(wg,R)) do
                 ex:=true;
                 for a in al[1..-2] do
#                     p:=trunc(pairing(wg,a)/2/coeff(weight0,lambda0));
#                     p1:=0;
#                     while pairing(t_b_action(-a,wg),a)>=1 do
#                         wg:=t_b_action(-a,wg);
#                         p1:=p1+1;
#                         ex:=false;
#                     od;
#                     print(p,p1);

                     p:=trunc(pairing(wg,a)/2/coeff(weight0,lambda0));
                     if p>0 then
                         wg:=t_b_action(-p*a,wg);
                         ex:=false;
                     end if;
                     if coeff(wg,delta)<0 then
                         break;
                     fi;
                 od;
             od;
             if is_in_main_chamber(wg,R) and coeff(wg,delta)>=coeff(weight0,delta) then
                 return subs(eps=0,wg)+eps*coeff(weight0,eps);
             fi;
         end do;
#         print(weight0);
         return 0;

     end;

fold_weight_strange:=proc(weight0,R,mg)
            local fp,a,al,wg,ex,weights,weight,p,p1;
                weights:=finite_orbit(subs(eps=0,weight0),finite_dimensional_root_system(R));
                al:=algebra_roots(R);

                for weight in weights do
                    wg:=t_b_action(convert(zip((x,y)->-x/2/coeff(weight,lambda0)*y,
                                               al[1..-2],root_coeffs(weight,R)),`+`),weight);
                    if is_in_main_chamber(wg,R) and coeff(wg,delta)>=coeff(weight0,delta) then
                        return subs(eps=0,wg)+eps*coeff(weight0,eps);
                    fi;
                od;
                return 0;
            end proc;

fold_weight3:=proc(weight,R,mg)
         local tr,al,rc,level,wg,orb,actor;
         al:=algebra_roots(R);
         level:=coeff(weight,lambda0);
         tr:=weight;
         while true do
             rc:=root_coeffs(tr,R);
             actor:=convert(zip((x,y)->-(trunc(x/2/level)+sign(frac(x/2/level)))*y,rc[1..-2],al[1..-2]),`+`);
             print(actor);
             if actor=0 or is_in_main_chamber(tr,R) then
                 break;
             fi;
             tr:=t_b_action(actor,tr);
             print(root_coeffs(tr,R));
         end do;
         orb:=finite_orbit(subs(eps=0,tr),finite_dimensional_root_system(R));
         print(map(x->root_coeffs(x,R),orb));
         wg:=select(x->is_in_main_chamber(x,R),
                    orb);
         if nops(wg)>0 then
             print(root_coeffs(wg[1],R));
#             return subs(eps=0,wg[1])+coeff(weight,eps)*eps;
             return wg[1];
         else
             return 0;
         end if;
     end proc;

fold_weight_simple:=proc (weight0,R,max_grade)
            local fp,a,al,wg,ex,weights,weight;
                weights:=finite_orbit(subs(eps=0,weight0),finite_dimensional_root_system(R));
                al:=algebra_roots(R);
                for weight in weights do
                    wg:=select(x->is_in_main_chamber(x,R) and coeff(x,delta)>=coeff(weight0,delta),
                           translations(weight,R,max_grade));
                    if nops(wg)>0 then
                        return subs(eps=0,wg[1])+coeff(weight0,eps)*eps;
                    end if;
                end do;
                return 0;
            end;

fold_star := (weight,star,R,max_grade)->select(x->x[2]<>0,# and coeff(x,delta)<=0,
                                                map(x->fold_weight_simple(x,R,max_grade),
                                                    map(x->weight+x,[op(star)])));


calc_m := proc(weight,star,iw,folded_stars,resTable)
         local w,w1;
         w:=subs(eps=0,weight);
         if coeff(w,delta)>0 then return 0; fi;
         if not assigned(resTable[w]) then
             w1:=subs(delta=0,w);
             if not assigned(folded_stars[w1]) then
#                 folded_stars[w1] := fold_star(w1,star,R,max_grade);
                 folded_stars[w1] := select(x->subs(eps=0,w1)<>subs(eps=0,x),fold_star_short(w1,star,iw));
#                 print(w1);
#                 print(folded_stars[w1]);
             fi;
             resTable[w]:=-add(coeff(x,eps)*calc_m(x+coeff(w,delta)*delta,star,iw,folded_stars,resTable),x in folded_stars[w1]);
         fi;
         return resTable[w];
     end;

string_function:=proc(weight,hweight,R,max_grade)
           local rt,k,theStar,fs,iw;
           rt:=table();
           fs:=table();
           rt[subs(eps=0,hweight)]:=1;
           iw:=map(x->subs(lambda0=0,x)+coeff(hweight,lambda0)*lambda0,dominant_weights(R,coeff(weight,lambda0)));
           theStar:=select(x->subs(eps=0,x)<>0,star(R,max_grade));
           return map(x->calc_m(weight-x*delta,theStar,iw,fs,rt),[$0..max_grade]);
       end;


fold_weight_short:=proc(weight,hweight)
         return hweight+(iprod(weight,weight)-iprod(hweight,hweight))/2/coeff(weight,lambda0)*delta+eps*coeff(weight,eps);
     end proc;

fold_star_short := (weight,star,iw)->select(x->x<>0,# and subs(eps=0,x)<>subs(eps=0,weight),# and coeff(x,delta)<=0,

                                                    map(w->
                                                        select(x->type(coeff(x,delta),integer),
                                                               map(x->
                                                                   fold_weight_short(w,subs(lambda0=0,x)+
                                                                                     coeff(weight,lambda0)*lambda0),iw))[1],
                                                        map(x->weight+x,[op(star)])));

root_coeffs:=(weight,R)->map(x->pairing(weight,x),algebra_roots(R));
coroot_coeffs:=(weight,R)->map(x->pairing(weight,x),algebra_co_roots(R));


weight_coeff:=(weight,R)->map(x->iprod(weight,x),weights(R));


to_simple_roots_base:=(w,R)->if type(w,'sequential') then
                                map(x->[LinearAlgebra[MatrixInverse](Matrix(coxeter['cartan_matrix'](finite_dimensional_root_system(R)))).Vector(root_coeffs(x,R)[1..-2]),
                                                                     coeff(x,delta),
                                        coeff(x,eps)],w);
                            else
                                [LinearAlgebra[MatrixInverse](Matrix(coxeter['cartan_matrix'](finite_dimensional_root_system(R)))).Vector(root_coeffs(w,R)[1..-2]),
                                                                     coeff(w,delta),
                                        coeff(w,eps)];
                            fi;


