read("branching.mpl"):

plotsetup(X11);

# Let's calculate branching rules for twisted subalgebra A^{(2)}_2 of
# algebra A^{(1)}_2
# Simple roots of A^{(1)}_2 are

A12_roots:=algebra_roots(A2);

# The injection is following:

A22_roots:=algebra_roots(tA2);

# then classical root of A^{(2)}_2 is equal to half of first classical root of
# A^{(1)}_2, and imaginary root of A^{(2)}_2 is equal to imaginary root of
# A^{(1)}_2

A22_roots_injected:=subs({e1=A12_roots[1]/2},A22_roots);

# Then, we should calculate positive roots of algebra A_2, for
# example to 1-st grade

pos_roots:=positive_roots(A2,25):
i:='i';
tA2_positive_roots:=n->[e1+i*delta+eps$i=0..n,2*e1+(2*i+1)*delta+eps$i=0..n/2,-2*e1+(2*i+1)*delta+eps$i=0..n/2,-e1+i*delta+eps$i=1..n,i*delta+eps$i=1..n];
A22_pos_roots:=tA2_positive_roots(25);

# We should project positive roots of A^{(1)}_2

pos_roots_projected:=projection(pos_roots,A22_roots_injected[1..-2],delta):

A22_pos_roots_injected:=subs({e1=A12_roots[1]/2,delta=delta},A22_pos_roots);

# Then we get projection of singular weights


sing:=projection(anomalous_points(lambda0,A2,25),A22_roots_injected[1..-2],delta):

# we could easily plot projection to check it
plot_projection:=proc(tt,real_roots,imaginary_root)
         plots[textplot](map(x->[coeff(x,delta)/coeff(imaginary_root,delta),pairing(x,real_roots[1])/2,coeff(x,eps)],[op(tt)]));
     end proc:

plot_projection(sing,A22_roots_injected,delta);


# now we can construct a fan of the injection:
#ft:=fan_table(pos_roots_projected,inj_roots(A22_pos_roots_injected)):

f:=fan(pos_roots_projected,inj_roots(A22_pos_roots_injected),25):


t:=table();

t[lambda0]:=1;

f1:=[eps,op(select(x->coeff(x,delta)<=25,f))];

plot_projection(projection(f1,A22_roots_injected[1..-2],delta),A22_roots_injected,delta);

# Something wrong in two following lines
# I really don't know how to find borders of injection diagram
#borders:=coeff(t_b_action(-finite_part(weyl_vector(tA2))/coeff(weyl_vector(tA2),lambda0),weyl_vector(tA2)),delta):

#calculate_branching_coefficient(-25*delta,f1,-(e2-e1)/2-eps,sing,t,rcurry(is_in_borders8,borders));

#iib:=x->coeff(x,e1)<25 and coeff(x,e2)<25;
borders:=external_border(sing);
iib:=rcurry(is_in_borders,borders);
calculate_branching_coefficient(lambda0-25*delta,f1,-(e2-e1)/2-eps,sing,t,iib);

#plot_projection(projection(map(x->x+eps*t[x],get_indices(t)),A22_roots_injected[1..-2],delta),A22_roots_injected,delta);
ap_draw2d(map(x->x+eps*t[x],get_indices(t)));
