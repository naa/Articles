read("affine.mpl"):

A3pr:=positive_roots(A3,1);

A1pr:=positive_roots(A1,1);

prepare_subtable:=proc(A1pr)
        local subtable, r, fpart, fpart0, grade;

            subtable:=table();
            for r in A1pr do
                fpart0:=coroot_coeffs(r,A1);
                grade:=coeff(r,delta);

                fpart:=[op(fpart0[1..-2]),0];
                if not assigned(subtable[[fpart,grade]]) then
                    subtable[[fpart,grade]]:=coeff(r,eps);
                else
                    subtable[[fpart,grade]]:=subtable[[fpart,grade]]+coeff(r,eps);
                end;
                fpart:=[0,op(fpart0[1..-2])];
                if not assigned(subtable[[fpart,grade]]) then
                    subtable[[fpart,grade]]:=coeff(r,eps);
                else
                    subtable[[fpart,grade]]:=subtable[[fpart,grade]]+coeff(r,eps);
                end;
            end;
            return subtable;
        end;


pmatrix:=[[2,1,2,1],
          [0,1,0,1],
          [2,1,0,1],
          [0,1,2,1]];

pm2:=[[1,0,1],
      [1,2,1]];

prepare_pretable:=proc(subtable0,pm2)
        local pretable, r, proj, grade, fpart, mult,sub,
            subtable, ind;
            subtable:=table(subtable0);

            pretable:=table();

            for r in A3pr do
                fpart:=coroot_coeffs(r,A3)[1..-2];
                grade:=coeff(r,delta);

                proj:=convert(linalg[multiply](pm2,fpart),'list');
                mult:=coeff(r,eps);
                ind:=[proj,grade];
                if assigned(subtable[ind]) and subtable[ind]<>0 then
                    sub:=min(mult,subtable[ind]);
                    mult:=mult-sub;
                    subtable[ind]:=subtable[ind]-sub;
                end;
                if not assigned(pretable[ind]) then
                    pretable[ind]:=mult;
                else
                    pretable[ind]:=pretable[ind]+mult;
                    print(pretable[ind]);
                end;

            end;
            return pretable;
        end;

prepare_funtable:=proc(pretable)
        local funtable, newtable, rg, j, r, elems;

            funtable:=table();
            funtable[[[0,0],0]]:=1;

            elems:=get_indices(pretable);
            for rg in elems do
                for j from 1 to pretable[rg] do
                    newtable:=table(funtable);
                    for r in get_indices(funtable) do
                        if (r-rg)[-1]>=-1 then
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
            return funtable;
        end;
