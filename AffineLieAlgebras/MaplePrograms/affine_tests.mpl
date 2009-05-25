read("affine.mpl"):


test_star_A2:={ [0,1,0,0,1] ,[ 2,1,0,0,-1] ,[1,0,0,0,1] ,[ 1,2,0,0,-1] ,[ 2,2,0,0,1] ,[ 3,1,0,1,1] ,[-1,1,0,1,-1] ,[ 1,3,0,1,1] ,[ 1,-1,0,1,-1] ,[ 3,3,0,1,-1],[ -1,-1,0,1,1] ,[ 3,4,0,2,1] ,[ 0,-2,0,2,1] ,[ 2,4,0,2,-1],[ -1,-2,0,2,-1] ,[ 4,3,0,2,1] ,[ -2,0,0,2,1] ,[ 4,2,0,2,-1],[ -2,-1,0,2,-1] ,[ 0,3,0,2,-1] ,[ 3,0,0,2,-1] ,[ -1,2,0,2,1],[ 2,-1,0,2,1] ,[ 0,4,0,4,1] ,[ -3,-2,0,4,1] ,[ 5,4,0,4,-1],[ 2,-2,0,4,-1] ,[ 4,0,0,4,1] ,[ -2,-3,0,4,1] ,[ 4,5,0,4,-1] ,[ -2,2,0,4,-1] ,[ -3,0,0,4,-1] ,[ 1,-3,0,5,1] ,[ 5,1,0,5,-1],[ 5,5,0,5,1] ,[ 1,5,0,5,-1] ,[ 0,-3,0,4,-1] ,[ 2,5,0,4,1],[ 5,2,0,4,1] , [ -3,-3,0,5,-1] ,[ -3,1,0,5,1] ,[ 6,4,0,6,1],[ 3,-2,0,6,1] ,[ -1,4,0,6,-1] ,[ -4,-2,0,6,-1] ,[ 4,6,0,6,1],[ -2,3,0,6,1] ,[ 4,-1,0,6,-1] ,[ -2,-4,0,6,-1] ,[ 3,6,0,6,-1],[ 6,3,0,6,-1] ,[ -4,-1,0,6,1] ,[ -1,-4,0,6,1] ,[ 6,1,0,8,1],[ -4,1,0,8,-1] ,[ 1,6,0,8,1] ,[ 1,-4,0,8,-1] ,[ 6,6,0,8,-1],[ -4,-4,0,8,1] ,[ 3,7,0,9,1] ,[ -3,-5,0,9,1] ,[ 5,7,0,9,-1],[ -1,-5,0,9,-1] ,[ 7,3,0,9,1] ,[ -5,-3,0,9,1] ,[ 7,5,0,9,-1],[ -5,-1,0,9,-1] ,[ -3,3,0,9,-1],[ 3,-3,0,9,-1] ,[ -1,5,0,9,1] ,[ 5,-1,0,9,1]}:

folded_star_A2_level1:={
 [ 0,0,0,1] ,[ 0,0,1,-2] ,[0,0,2,-1] ,[ 0,0,3,2] ,[ 0,0,4,1] ,[0,0,5,2] ,
[ 0,0,7,-2] ,[ 0,0,8,-2] ,[ 0,0,9,1],[ 0,0,10,1] ,[ 0,0,13,2] ,[ 0,0,14,3] ,
[ 0,0,15,-2] ,[ 0,0,16,2] ,[ 0,0,19,-2],[ 0,0,20,-2]
}:

# to get folded star in the same form from fold_star_short use
# fs:=fold_star_short(2*lambda0,test_star_A2_expanded,[iw[1],iw[-1]]);
# map(x->[LinearAlgebra[MatrixInverse](m).Vector(root_coeffs(x,A2)[1..2]),coeff(x,delta),coeff(x,eps)],fs);


roots_A2:=algebra_roots(A2);
test_star_A2_expanded:=map(x->x[1]*roots_A2[1]+x[2]*roots_A2[2]+x[3]*lambda0+x[4]*delta-x[5]*eps,test_star_A2):

run_tests:=proc()
        st:=map(x->[op(root_coeffs(x,A2)[1..2]),level(x),coeff(x,delta),coeff(x,eps)],star(A2,9));
        print(evalb(st=test_star_A2));
        r:=multiplicities(lambda0,A1,10);
        print("Testing multiplicities");
        print(evalb(filter_string_function(lambda0,r)[1..10]=[1,1,2,3,5,7,11,15,22,30]));
#string_function(_

    end proc:


translation_with_coeff := proc(lambda,b,k)
        return lambda+iprod(lambda,delta)*b*k-(k*k*iprod(b,b)/2*iprod(lambda,delta)+k*iprod(lambda,b))*delta;
    end proc;


test_string_f := proc(weight,hweight)
local st,rt,w;
                R:=A2;
                max_grade:=9;
    rt:=table();
    rt[subs(eps=0,hweight)]:=1;
    w:=fold_weight_short(weight,hweight);
    st:=fold_star_short(hweight,hweight,select(x->subs(eps=0,x)<>0,test_star_A2_expanded));
    return map(x->calc(weight-x*delta,st,rt),[$0..max_grade]);
end proc;


test_string_func:=proc(weight,hweight)
            local rt,k,theStar,fs,iw;
                R:=A2;
                max_grade:=9;
           rt:=table();
           fs:=table();
           rt[subs(eps=0,hweight)]:=1;
           iw:=map(x->subs(lambda0=0,x)+coeff(hweight,lambda0)*lambda0,dominant_weights(R,coeff(weight,lambda0)));
           theStar:=select(x->subs(eps=0,x)<>0,test_star_A2_expanded);
           return map(x->calc_m(weight-x*delta,theStar,iw,fs,rt),[$0..max_grade]);
       end;

