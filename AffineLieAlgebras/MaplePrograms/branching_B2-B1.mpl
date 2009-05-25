read("branching.mpl"):

plotsetup(X11);

# Let's calculate branching rules for twisted subalgebra A^{(2)}_2 of
# algebra A^{(1)}_2
# Simple roots of A^{(1)}_2 are

B2_roots:=algebra_roots(B2);

# The injection is following:

B1_roots:=algebra_roots(B1);

# then classical root of A^{(2)}_2 is equal to half of first classical root of
# A^{(1)}_2, and imaginary root of A^{(2)}_2 is equal to imaginary root of
# A^{(1)}_2

B1_roots_injected:=subs({e1=e2-e1},B1_roots);

# Then, we should calculate positive roots of algebra A_2, for
# example to 1-st grade

pos_roots:=positive_roots(B2,2):
B1_pos_roots:=positive_roots(B1,2):

pos_roots_projected:=projection(pos_roots,B1_roots_injected[1..-2],delta):

B1_pos_roots_injected:=projection(subs({e1=e2-e1},B1_pos_roots),B1_roots_injected[1..-2],delta):

# Then we get projection of singular weights

sing:=projection(anomalous_points(lambda0,B2,5),B1_roots_injected[1..-2],delta):
#sing:=[lambda0+eps,sing[]];
# we could easily plot projection to check it
plot_projection:=proc(tt,real_roots,imaginary_root)
         plots[textplot](map(x->[coeff(x,delta)/coeff(imaginary_root,delta),pairing(x,real_roots[1])/2,coeff(x,eps)],[op(tt)]));
     end proc:

plot_projection(sing,B1_roots_injected,delta);

# now we can construct a fan of the injection:
#ft:=fan_table(pos_roots_projected,inj_roots(B12_pos_roots_injected)):
f:=fan(select(x->subs(eps=0,x)<>0,pos_roots_projected),inj_roots(B1_pos_roots_injected)):

t:=table();
t[lambda0]:=1;
f1:=[eps,op(select(x->coeff(x,delta)<=5,f))];

plot_projection(projection(f1,B1_roots_injected[1..-2],delta),B1_roots_injected,delta);

# Something wrong in two following lines
# I really don't know how to find borders of injection diagram
#borders:=coeff(t_b_action(-finite_part(weyl_vector(tB1))/coeff(weyl_vector(tB1),lambda0),weyl_vector(tB1)),delta):
borders:=coeff(t_b_action(-finite_part(weyl_vector(B2)),weyl_vector(B2)),delta):
#borders:=0;
calculate_branching_coefficient(-5*delta+lambda0,f1,-(e2-e1)/2-eps,sing,t,rcurry(is_in_borders8,borders));

#Here we have some random limits
#iib:=x->coeff(x,e1)<10 and coeff(x,e2)<10;
#calculate_branching_coefficient(-10*delta,f1,-(e2-e1)/2-eps,sing,t,iib);

plot_projection(projection(map(x->x+eps*t[x],get_indices(t)),B1_roots_injected[1..-2],delta),B1_roots_injected,delta);
