read("affine.mpl"):

A3pr:=positive_roots(A3,10);

A1pr:=positive_roots(A1,10);
subtable:=table();
for r in A1pr do
    fpart0:=coroot_coeffs(r,A1);
    grade:=coeff(r,delta);

    fpart:=[op(fpart0),0,0];
    if not assigned(subtable[[fpart,grade]]) then
        subtable[[fpart,grade]]:=coeff(r,eps);
    else
        subtable[[fpart,grade]]:=subtable[[fpart,grade]]+coeff(r,eps);
    end;
    fpart:=[0,0,op(fpart0)];
    if not assigned(subtable[[fpart,grade]]) then
        subtable[[fpart,grade]]:=coeff(r,eps);
    else
        subtable[[fpart,grade]]:=subtable[[fpart,grade]]+coeff(r,eps);
    end;
end;



pmatrix:=[[2,1,2,1],
          [0,1,0,1],
          [2,1,0,1],
          [0,1,2,1]];

pretable:=table();

for r in A3pr do
    fpart:=coroot_coeffs(r,A3);
    grade:=coeff(r,delta);

    proj:=linalg[multiply](pmatrix,fpart);
    mult:=coeff(r,eps);
    if assigned(subtable[[fpart,grade]]) then
        sub:=min(mult,subtable[[fpart,grade]]);
        mult:=mult-sub;
        subtable[[fpart,grade]]:=subtable[[fpart,grade]]-sub;
    end;
    if not assigned(pretable[[fpart,grade]]) then
        pretable[[fpart,grade]]:=mult;
    else
        pretable[[fpart,grade]]:=pretable[[fpart,grade]]+mult;
    end;
end;

funtable:=table();
funtable[[[0,0,0,0],0]]:=1;

for rg in get_indices(pretable) do
    for j from 1 to pretable[rg] do
        newtable:=table(funtable);
        for r in get_indices(funtable) do
            if (r-rg)[-1]>=-10 then
                if not assigned(newtable[r-rg]) then
                    newtable[r-rg]:=-funtable[r];
                else
                    newtable[r-rg]:=newtable[r-rg]-funtable[r];
                end;
            end;
        end;
        funtable:=newtable;
    end;
end;

