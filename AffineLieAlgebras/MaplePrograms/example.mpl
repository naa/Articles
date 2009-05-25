read("affine.mpl"):

# Some basics


# Roots for given algebra

al:=algebra_roots(A2);

# Weights

weights(A2);

# Dominant weights of given level

dw:=dominant_weights(A2,2);

# to get a vector in simple roots base use function
# to_simple_roots_base

to_simple_roots_base(dw,A2);

# At first we want to get a star (fan). There is no need to do it
# explicitely if we need just weight multiplicities

star(A2,9);

# Also we can fold the star, which spans from some weight

fs:=fold_star(dw[1],star(A2,5),A2,5);
to_simple_roots_base(fs,A2);

# or we can use equivalent function which is slightly faster

fs:=fold_star_short(dw[2],star(A2,5),dw);
to_simple_roots_base(fs,A2);

# To calculate weight multiplicities with straightforward algorithm
# (without folding the star) we can use function multiplicities, which
# returns results as Maple table

r:=multiplicities(dw[1],A2,5);

map(x->[x,r[x]],[indices(r,'nolist')]);

# To filter out string function coefficients

filter_string_function(dw[1],r);

# We can also calculate string function coefficients using more
# advanced method with folding

string_function(dw[2],dw[1],A2,5);

string_function(dw[1],dw[1],A2,5);

# And finally, we can get all the string function coefficients for A2
# level 2 (for example) at once
map(y->[y,
        map(x->[x,string_function(x,y,A2,5)],dw)],dw);

# To get branching rules for branching of module $L^{[1,0,0]}$ of
# algebra $A_2$ to regular sub-algebra $A_1$ we can use following commands:

read("branching.mpl"):
wg:=weights(A2);
tt1:=branching_rules(wg[-1],A2,A1,{e1=alpha[1]},15);
tt1[lambda0];
tt1[lambda0+e1/2];
