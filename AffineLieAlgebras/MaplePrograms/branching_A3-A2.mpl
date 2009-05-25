read("branching.mpl"):

plotsetup(X11);

# Let's calculate branching rules for twisted subalgebra A^{(2)}_2 of
# algebra A^{(1)}_2
# Simple roots of A^{(1)}_2 are

A3_roots:=algebra_roots(A3);

# The injection is following:

A2_roots:=algebra_roots(A2);

# then classical root of A^{(2)}_2 is equal to half of first classical root of
# A^{(1)}_2, and imaginary root of A^{(2)}_2 is equal to imaginary root of
# A^{(1)}_2

A2_roots_injected:=A2_roots;

# Then, we should calculate positive roots of algebra A_2, for
# example to 1-st grade

pos_roots:=positive_roots(A3,10):
A2_pos_roots:=positive_roots(A2,10):

pos_roots_projected:=projection(pos_roots,A2_roots_injected[1..-2],delta):

A2_pos_roots_injected:=A2_pos_roots;

# Then we get projection of singular weights

sing:=projection(anomalous_points(lambda0,A3,5),A2_roots_injected[1..-2],delta):

# we could easily plot projection to check it
plot_projection:=proc(tt,real_roots,imaginary_root)
         plots[textplot3d](map(x->[coeff(x,delta)/coeff(imaginary_root,delta),pairing(x,real_roots[1])/2,pairing(x,real_roots[2])/2,coeff(x,eps)],[op(tt)]));
     end proc:

plot_projection(sing,A2_roots_injected,delta);

# now we can construct a fan of the injection:
#ft:=fan_table(pos_roots_projected,inj_roots(A22_pos_roots_injected)):
f:=fan(pos_roots_projected,inj_roots(A2_pos_roots_injected)):

t:=table();

f1:=[eps,op(select(x->coeff(x,delta)<=5,f))];

plot_projection(projection(f1,A2_roots_injected[1..-2],delta),A2_roots_injected,delta);

# Something wrong in two following lines
# I really don't know how to find borders of injection diagram
#borders:=coeff(t_b_action(-finite_part(weyl_vector(tA2))/coeff(weyl_vector(tA2),lambda0),weyl_vector(tA2)),delta):
borders:=-coeff(t_b_action(finite_part(weyl_vector(A3)),weyl_vector(A3)),delta):
#borders:=0;
calculate_branching_coefficient(-5*delta+lambda0,f1,-(e2-e1)/2-eps,sing,t,rcurry(is_in_borders8,borders));

#Here we have some random limits
#iib:=x->coeff(x,e1)<10 and coeff(x,e2)<10;
#calculate_branching_coefficient(-10*delta,f1,-(e2-e1)/2-eps,sing,t,iib);

plot_projection(projection(map(x->x+eps*t[x],get_indices(t)),A2_roots_injected[1..-2],delta),A2_roots_injected,delta);
