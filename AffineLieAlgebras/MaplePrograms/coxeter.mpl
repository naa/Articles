#
# Coxeter and weyl packages: version 2.4, vanilla edition.
# This version/edition requires one of the following versions of Maple:
#   V R1, R2, R3, R4, R5, 6, 7, 8, or 9.
#
# This is *not* a Maple worksheet.
#
# After loading this file during a Maple session, each package function
# may be accessed using one of the calling sequences
#
#   coxeter[<functionname>](<arguments>) or
#   weyl[<functionname>](<arguments>).
#
# In order to use package functions in the abbreviated form
#
#   <functionname>(<arguments>),
#
# run the commands 'withcoxeter()' and 'withweyl()' (respectively) after
# loading this file. If there is a conflict between the names of one of
# these functions and another name in the same session, a warning is printed.
#
# In order to use the abbreviated form for a subset of the procedures in
# one of the packages, run the commands
#
#   withcoxeter(<name1>,<name2>,...)  or  withweyl(<name1>,<name2>,...).
#
# For an introduction, see
#   http://www.math.lsa.umich.edu/~jrs/software/coxeter.ps
# For documentation on the individual functions, see
#   http://www.math.lsa.umich.edu/~jrs/software/coxeterhelp.html
#
# Copyright (c) 2004 by John R. Stembridge
#
#########################################################################
#
# In floating-point computations, angles and distances that
# are less than the following are assumed to be zero.
#
`coxeter/default`['epsilon']:=0.001:
#
coxeter:=table():
weyl:=table():
#
`coxeter/read`:=proc(R,a,b) eval(cat(`coxeter/`,a,'`/`',R)) end:
#
# assign short names, printing warnings if conflicts occur.
#
withcoxeter:=proc() local install,f,lpr;
  if [op(I)]<>[1] then lpr:=lprint else # maple6 has a broken lprint
    lpr:=proc(x,y) printf(`%0.70s\n`,cat(x,`  `,y)) end fi;
  install:=proc(x,wrn)
    if not assigned(coxeter[x]) then
      ERROR(cat(x,` is not a top level function in coxeter`))
    elif eval(x)<>eval(coxeter[x]) then
      if x='index' then unprotect(x) fi;
      if assigned(x) then wrn(`Warning: new definition for`,x) fi;
      assign(x,coxeter[x])
    fi; x
  end;
  if nargs>0 then map(install,[args],lpr) else
    f:=proc() map(op,[args]) end; # hack the names w/o full evaluation!
    map(install,f(indices(coxeter)),lpr)
  fi
end:
#
withweyl:=proc() local install,f,lpr;
  if [op(I)]<>[1] then lpr:=lprint else # maple6 has a broken lprint
    lpr:=proc(x,y) printf(`%0.70s\n`,cat(x,`  `,y)) end fi;
  install:=proc(x,wrn)
    if not assigned(weyl[x]) then
      ERROR(cat(x,` is not a top level function in weyl`))
    elif eval(x)<>eval(weyl[x]) then
      if x='tensor' then unprotect(x) fi;
      if assigned(x) then wrn(`Warning: new definition for`,x) fi;
      assign(x,weyl[x])
    fi; x
  end;
  if nargs>0 then map(install,[args],lpr) else
    f:=proc() map(op,[args]) end; # hack the names w/o full evaluation!
    map(install,f(indices(weyl)),lpr)
  fi
end:
#
# base(R) returns a base of simple roots for the root system R.
# co_base(R) returns the simple co-roots for base(R).
# In non-xtal cases, floating-point coordinates are used.
#
`coxeter/base`:=proc(X) local sys,R,i,S,pi;
  if X=1 then []
    elif type(X,'list') then X
    elif type(X,'matrix') then
      R:=coxeter['name_of'](X,'pi');
      S:=`coxeter/base`(R);
      subsop(seq(pi[i]=S[i],i=1..nops(pi)),S)
    else i:=1;
      sys:=sort([op(indets(X))],`coxeter/cox_order`);
      sys:=[seq(R$degree(X,R),R=sys)];
      [seq(`coxeter/base/simple`(R,i,'i'),R=sys)]
  fi
end:
#
`coxeter/co_base`:=proc(R) local r;
  [seq(2*r/coxeter['iprod'](r,r),r=`coxeter/base`(R))]
end:
#
`coxeter/base/simple`:=proc(R) local r,i,j,L,S,a;
  r:=coxeter['rank'](R); i:=args[2];
  if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
#  if R='A1' then
#    S:=cat('e',i)
#  el
    if L='A' then
    S:=seq(cat('e',j)-cat('e',j+1),j=i..i+r-1); i:=i+1
  elif L='B' then
    S:=seq(cat('e',j)-cat('e',j+1),j=i..i+r-2),cat('e',i+r-1)
  elif L='C' then
    S:=seq(cat('e',j)-cat('e',j+1),j=i..i+r-2),2*cat('e',i+r-1)
  elif L='D' and r>1 then
    S:=seq(cat('e',j)-cat('e',j+1),j=i..i+r-2),cat('e',i+r-2)+cat('e',i+r-1)
  elif L='E' and r>2 and r<9 then
    S:=
      seq(cat('e',j)-cat('e',j+1),j=i..i+r-3),
      cat('e',i+r-2)+cat('e',i+r-3),
      (-cat('e',i)-cat('e',i+1)-cat('e',i+2)-cat('e',i+3)-cat('e',i+4)
       -cat('e',i+5)-cat('e',i+6)-cat('e',i+7))/2; i:=i+8-r
  elif R='F4' then
    S:=
      cat('e',i)-cat('e',i+1),
      cat('e',i+1)-cat('e',i+2),
      cat('e',i+2),
      (-cat('e',i)-cat('e',i+1)-cat('e',i+2)+cat('e',i+3))/2
  elif R='G2' then
    S:=-2*cat('e',i)+cat('e',i+1)+cat('e',i+2),cat('e',i)-cat('e',i+1)
  elif R='H3' or R='H4' then a:=evalf((1+sqrt(5))/2);
    S:=-a*cat('e',i)+(a-1)*cat('e',i+1)-cat('e',i+2);
    S:=2*cat('e',i),S,2*cat('e',i+2);
    if r=4 then S:=S,-a*cat('e',i+1)-cat('e',i+2)+(a-1)*cat('e',i+3) fi;
  elif L='I2' and type(op(R)-1,'posint') then a:=op(R);
    S:=op(evalf([cat('e',i),-cos(Pi/a)*cat('e',i)+sin(Pi/a)*cat('e',i+1)]))
  else
    ERROR(cat(R,` is not a recognized root system name`))
  fi;
  assign(args[3],i+r); S
end:
#
# cartan_matrix(R) - generate the Cartan matrix of an xtal root system
# index(R) - determine the index of connection of an xtal root system
#
`coxeter/cartan_matrix`:=proc(M) local S,r,s,coS;
  if type(M,'matrix') and M[1,1]=2 then S:=M else
    S:=coxeter['base'](M);
    if S=[] then RETURN(S) fi;
    coS:=coxeter['co_base'](S);
    S:=array([seq([seq(coxeter['iprod'](r,s),s=S)],r=coS)])
  fi;
  if type(S,'matrix'('integer')) then eval(S)
    else ERROR(`not crystallographic`) fi
end:
#
`coxeter/index`:=proc(R) local M;
  M:=`coxeter/cartan_matrix`(R);
  if type(M,'matrix') then linalg['det'](M) else 1 fi
end:
#
# char_poly(w,R,q) returns det(1-q*w), relative to the reflection action.
# Modified to run faster for numeric q.
#
`coxeter/char_poly`:=proc(w) local S,z,vars,A,v,n,i,Z;
  S:=coxeter['base'](args[2]); A:=NULL;
  vars:=[op(indets(S))]; n:=nops(vars);
  for v in vars do;
    z:=coxeter['reflect'](seq(S[i],i=w),v);
    A:=A,[seq(coeff(z,vars[i]),i=1..n)];
  od;
  if n=0 then RETURN(1) fi;
  A:=array(1..n,1..n,[A]);
  if nargs>2 then z:=args[3] else z:=q fi;
  if n=nops(S) then
    A:=`coxeter/char_poly/mat`(A,n,z);
    expand(linalg['det'](A))
  elif type(z,'numeric') and z-1<>0 then
    A:=`coxeter/char_poly/mat`(A,n,z);
    linalg['det'](A)/(1-z)^(n-nops(S))
  else
    A:=linalg['det'](`coxeter/char_poly/mat`(A,n,Z));
    subs(Z=z,expand(quo(A,(1-Z)^(n-nops(S)),Z)))
  fi;
end:
#
# this is faster than evalm()
#
if `+`(0)=0 then # we are using Maple V.4 or later
  `coxeter/char_poly/mat`:=proc(B,n,z)
    linalg['matadd'](linalg['band']([1],n),B,1,-z) end
else
  `coxeter/char_poly/mat`:=proc(B,n,z)
    linalg['add'](linalg['band']([1],n),B,1,-z) end
fi:
#
# class_rep(w,R) returns a canonical rep for the conjugacy class of w.
# class_rep(R) returns a list of canonical reps for all classes.
#
`coxeter/class_rep`:=proc() local n,R,N,res,new,sys,i,j,k,w,L,pi,w0;
  R:=coxeter['name_of'](args[nargs],'pi'); N:=0;
  sys:=sort([op(indets(R))],`coxeter/cox_order`);
  sys:=[seq(i$degree(R,i),i=sys)];
  if nargs=1 then res:=[];
    for R in sys do
      n:=coxeter['rank'](R);
      if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
      new:=`coxeter/class_rep/list`(L,n,R);
      res:=seq(seq([op(i),seq(k+N,k=j)],i=[res]),j=new);
      N:=N+n;
    od;
  else res:=NULL;
    w0:=subs({seq(pi[i]=i,i=1..nops(pi))},args[1]);
    for R in sys do
      n:=coxeter['rank'](R);
      w:=map(proc(x,y,z) if x>y and x<=z then x-y fi end,w0,N,N+n);
      if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
      new:=`coxeter/mytype`(w,L,n,R);
      if L='A' then new:=`coxeter/par2A`(new)
        elif L='B' or L='C' then new:=`coxeter/par2B`(new)
        elif L='D' then new:=`coxeter/par2D`(new)
        elif L='E' and n<6 then
          new:=`coxeter/class_rep`(w,coxeter['base'](R))
        elif member(L,['E','F','H']) then
          new:=`coxeter/read`(R,'class_rep','classes')[new]
        else new:=`coxeter/class_rep/list`(L,n,R)[new]
      fi;
      res:=res,seq(i+N,i=new); N:=N+n;
    od;
  fi;
  subs({seq(i=pi[i],i=1..N)},[res]);
end:
#
`coxeter/class_rep/list`:=proc(L,n,R) local i,j,mu,nu,m,sa;
  if L='A' then
    [seq(`coxeter/par2A`(mu),mu=combinat['partition'](n+1))]
  elif L='B' or L='C' then
    [seq(seq(seq(`coxeter/par2B`(mu,nu),mu=combinat['partition'](n-i)),
      nu=combinat['partition'](i)),i=0..n)]
  elif L='D' and n>1 then
    sa:=seq(seq(seq(`coxeter/par2D`(mu,nu),mu=combinat['partition'](n-i)),
      nu=combinat['partition'](i)),i=0..n);
    if modp(n,2)=0 then [sa,seq(`coxeter/par2D`(map(x->2*x,mu),[],1),
      mu=combinat['partition'](n/2))] else [sa] fi;
  elif L='E' and n<6 then
    `coxeter/class_rep`(coxeter['base'](R))
  elif member(L,['E','F','H']) then
    `coxeter/read`(R,'class_rep','classes')
  elif R='G2' then [[],[1],[2],[1,2],[1,2,1,2],[1,2,1,2,1,2]]
  elif L='I2' and type(op(R)-1,'posint') then m:=op(R);
    if modp(m,2)=0 then
      [[],[1],[2],seq([seq(1+modp(i,2),i=0..2*j-1)],j=1..m/2)]
    else
      [[],[1],seq([seq(1+modp(i,2),i=0..2*j-1)],j=1..(m-1)/2)]
    fi;
  else
    ERROR(cat(R,` is not a recognized root system name`))
  fi;
end:
#
# class_size(w,R) returns the size of the conjugacy class of w.
# class_size(R) returns the sizes of all conjugacy classes.
#
`coxeter/class_size`:=proc() local n,R,N,res,new,sys,i,j,w,L,w0,pi;
  R:=coxeter['name_of'](args[nargs],'pi');
  sys:=sort([op(indets(R))],`coxeter/cox_order`);
  sys:=[seq(i$degree(R,i),i=sys)];
  if nargs=1 then res:=[1];
    for R in sys do
      n:=coxeter['rank'](R);
      if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
      new:=`coxeter/class_size/list`(L,n,R);
      res:=[seq(seq(i*j,i=res),j=new)];
    od;
  else
    N:=0; res:=1;
    w0:=subs({seq(pi[i]=i,i=1..nops(pi))},args[1]);
    for R in sys do
      n:=coxeter['rank'](R);
      w:=map(proc(x,y,z) if x>y and x<=z then x-y fi end,w0,N,N+n);
      if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
      new:=`coxeter/mytype`(w,L,n,R);
      if L='A' then new:=`coxeter/par2A`(new,0)
        elif L='B' or L='C' then new:=`coxeter/par2B`(new,0)
        elif L='D' then new:=`coxeter/par2D`(new[1..2],0)
        elif L='E' and n<6 then
          new:=`coxeter/class_size`(w,coxeter['base'](R))
        elif member(L,['E','F','H']) then
          new:=`coxeter/read`(R,'class_size','sizes')[new]
        else new:=`coxeter/class_size/list`(L,n,R)[new]
      fi;
      res:=res*new; N:=N+n;
    od
  fi; res
end:
#
`coxeter/class_size/list`:=proc(L,n,R) local i,mu,nu,m,sa;
  if L='A' then
    [seq(`coxeter/par2A`(mu,0),mu=combinat['partition'](n+1))]
  elif L='B' or L='C' then
    [seq(seq(seq(`coxeter/par2B`(mu,nu,0),mu=combinat['partition'](n-i)),
      nu=combinat['partition'](i)),i=0..n)]
  elif L='D' and n>1 then
    sa:=seq(seq(seq(`coxeter/par2D`(mu,nu,0),mu=combinat['partition'](n-i)),
      nu=combinat['partition'](i)),i=0..n);
    if modp(n,2)=0 then [sa,seq(`coxeter/par2D`(map(x->2*x,mu),[],0),
      mu=combinat['partition'](n/2))] else [sa] fi
  elif L='E' and n<6 then
    `coxeter/class_size`(coxeter['base'](R))
  elif member(L,['E','F','H']) then
    `coxeter/read`(R,'class_size','sizes')
  elif R='G2' then [1,3,3,2,2,1]
  elif L='I2' and type(op(R)-1,'posint') then m:=op(R);
    if modp(m,2)=0 then [1,m/2,m/2,2$(m/2-1),1]
      else [1,m,2$((m-1)/2)] fi;
  else
    ERROR(cat(R,` is not a recognized root system name`))
  fi;
end:
#
# cox_matrix(R) returns the Coxeter matrix of R.
# The empty matrix is represented as [].
#
`coxeter/cox_matrix`:=proc(M) local S,n,i,j,EPS,L;
  if type(M,'matrix') and M[1,1]=1 then eval(M)
    elif type(M,'matrix') then
      n:=linalg['rowdim'](M);
      array('symmetric',1..n,1..n,[seq([seq(
        `coxeter/cox_matrix/chk`({M[i,j],M[j,i]}),j=1..i-1),1],i=1..n)])
    else
      S:=coxeter['base'](M); n:=nops(S);
      if S=[] then RETURN(S) else EPS:=`coxeter/default`['epsilon'] fi;
      L:=[seq(coxeter['iprod'](i,i),i=S)];
      S:=[seq(evalf(S[i]/sqrt(L[i])),i=1..n)];
      array('symmetric',1..n,1..n,[seq([seq(`coxeter/cox_matrix/entry`
        (S[i],S[j],L[i],L[j],EPS),j=1..i-1),1],i=1..n)])
  fi
end:
#
`coxeter/cox_matrix/chk`:=proc(x)
  if x={0} then 2 elif x={-1} then 3
    elif x={-1,-2} then 4 elif x={-1,-3} then 6
    elif type(x,'set'('negint')) then infinity
    else ERROR(`not a valid Cartan matrix`)
  fi
end:
#
`coxeter/cox_matrix/entry`:=proc(r,s,lr,ls,EPS) local m,rm;
  m:=evalf(Pi/arccos(-coxeter['iprod'](r,s))); rm:=round(m);
  if abs(rm-m)>EPS then ERROR(`not a valid base`)
    elif modp(rm,2)=1 and abs(lr-ls)>EPS then ERROR(`not a valid base`)
    else rm
  fi
end:
#
# cox_order(R1,R2) returns true if R1 precedes R2 in the canonical ordering
# of irreducible components of a root system
#
`coxeter/cox_order`:=proc(R1,R2) local L1,L2;
  if type(R1,'indexed') then L1:=op(0,R1) else L1:=substring(R1,1..1) fi;
  if type(R2,'indexed') then L2:=op(0,R2) else L2:=substring(R2,1..1) fi;
  if L1<>L2 then lexorder(L1,L2)
    elif L1='I2' then evalb(op(R1)<=op(R2))
    elif length(R1)=length(R2) then lexorder(R1,R2)
    else evalb(length(R1)<length(R2))
  fi;
end:
#
# cprod(chi,theta,R) computes the inner product of class functions chi
#  and theta. Each class function should be a list (or array) of function
#  values on conjugacy classes ordered in the same way as class_rep(R).
# Note: This is the *real* (symmetric, not Hermitian) inner product.
#
`coxeter/cprod`:=proc(f,g,R) local i,sz;
  sz:=coxeter['class_size'](R);
  sz:=[seq(f[i]*g[i]*sz[i],i=1..nops(sz))];
  convert(sz,`+`)/coxeter['size'](R)
end:
#
#  size(R)=size of the Coxeter group W(R)
#  exponents(R)=list of exponents
#  num_refl(R)=number of reflections in the Coxeter group
#  cox_number(R)=order of the Coxeter element
#  degrees(R)=degrees of the fundamental invariants
#
`coxeter/size`:=proc(R) convert(`coxeter/degrees`(R),`*`) end:
`coxeter/exponents`:=proc(R) map(x->x-1,`coxeter/degrees`(R)) end:
`coxeter/num_refl`:=proc(R) convert(`coxeter/exponents`(R),`+`) end:
`coxeter/cox_number`:=proc(X) local R,r;
  R:=coxeter['name_of'](X);
  ilcm(seq(max(`coxeter/degrees/irr`(r,r)),r=indets(R)))
end:
#
`coxeter/degrees`:=proc(X) local R;
  R:=coxeter['name_of'](X);
  sort(map(`coxeter/degrees/irr`,[op(indets(R))],R))
end:
#
`coxeter/degrees/irr`:=proc(R,X) local L,r,i,d;
  r:=coxeter['rank'](R);
  if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
  if L='A' or R='E4' then d:=[$2..(r+1)]
    elif L='B' or L='C' then d:=[seq(2*i,i=1..r)]
    elif L='D' or R='E5' then d:=[seq(2*i,i=1..r-1),r]
    elif R='E6' then d:=[2,5,6,8,9,12]
    elif R='E7' then d:=[2,6,8,10,12,14,18]
    elif R='E8' then d:=[2,8,12,14,18,20,24,30]
    elif R='F4' then d:=[2,6,8,12]
    elif R='G2' then d:=[2,6]
    elif R='H3' then d:=[2,6,10]
    elif R='H4' then d:=[2,12,20,30]
    elif L='I2' and type(op(R)-1,'posint') then d:=[2,op(R)]
    elif R='E3' then d:=[2,2,3]
    else ERROR(cat(R,` is not a recognized root system name`))
  fi;
  op(map(op,[d$degree(X,R)]));
end:
#
# diagram(X) prints a Coxeter or Dynkin diagram for X, with
# the nodes numbered in the order used by base(X).
#
`coxeter/diagram`:=proc() local sys,R,L,i,r,j,pi,c,lpr;
  if [op(I)]=[1] then lpr:=proc() printf(`%0.75s\n`,args) end
    else lpr:=proc() lprint(args) end fi;
  lpr(` `); j:=1;
  R:=coxeter['name_of'](args[1],'pi');
  sys:=sort([op(indets(R))],`coxeter/cox_order`);
  sys:=[seq(r$degree(R,r),r=sys)];
  for R in sys do;
    r:=coxeter['rank'](R);
    if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1)  fi;
    if L='A' or R='I2[3]' or R='B1' or R='C1' then
      lpr(cat(`   `,pi[j],seq(cat(`---`,pi[i]),i=j+1..j+r-1)))
    elif L='B' then
      lpr(cat(`   `,pi[j],`=<=`,pi[j+1],seq(cat(`---`,pi[i]),i=j+2..j+r-1)))
    elif L='C' then
      lpr(cat(`   `,pi[j],`=>=`,pi[j+1],seq(cat(`---`,pi[i]),i=j+2..j+r-1)))
    elif L='D' and r>1 then
      c:=`      `,` `$length(pi[j]); lpr(cat(c,pi[j+1]));
      if r>2 then lpr(cat(c,`|`)) else lpr(` `) fi;
      lpr(cat(`   `,pi[j],seq(cat(`---`,pi[i]),i=j+2..j+r-1)))
    elif L='E' and r>2 then
      c:=`         `,` `$(length(pi[j])+length(pi[j+2]));
      lpr(cat(c,pi[j+1]));
      if r>3 then lpr(cat(c,`|`)) else lpr(` `) fi;
      lpr(cat(`   `,pi[j],seq(cat(`---`,pi[i]),i=j+2..j+r-1)))
    elif L='F' and r>2 then
      lpr(cat(`   `,pi[j],`---`,pi[j+1],`=<=`,pi[j+2],
        seq(cat(`---`,pi[i]),i=j+3..j+r-1)))
    elif R='G2' then
      lpr(cat(`   `,pi[j],`=<<=`,pi[j+1]))
    elif L='H' and r>1 then
      lpr(cat(`   `,pi[j],`--(5)--`,pi[j+1],
        seq(cat(`---`,pi[i]),i=j+2..j+r-1)))
    elif L='I2' and type(op(R)-1,'posint') then
      if op(R)=2 then c:=`   ` else c:=`--(`,op(R),`)--` fi;
      lpr(cat(`   `,pi[j],c,pi[j+1]))
    else
      ERROR(cat(R,` is not a recognized root system name`))
    fi;
    j:=j+r; lpr(` `);
  od;
end:
#
# If R is reducible or non-xtal, the result is simply a long root in the
# dominant chamber. (Highest roots are not defined in these cases).
# We could break the loop as soon as a longer root is found, but in
# reducible cases this might be short.
#
`coxeter/highest_root`:=proc(R) local S,v,l,ll,r;
  S:=coxeter['base'](R); ll:=0; v:=0;
  for r in S do
    l:=coxeter['iprod'](r,r);
    if l>ll then v:=r; ll:=l fi;
  od;
  coxeter['vec2fc'](v,S);
end:
#
# induce a class function f0 from some reflection subgroup to W(R).
# J can be either a sublist of [1,...,n] indicating a parabolic subgroup,
# or a base for some (not necessarily parabolic) reflection subgroup.
#
`coxeter/induce`:=proc(f0,J,R) local S,S0,w,gens,cc,cc0,sz,sz0,res,v,i,j,f;
  S:=coxeter['base'](R); S0:=J;
  if type(S0,'list'('integer')) then S0:=[seq(S[i],i=S0)] fi;
  v:=coxeter['interior_pt'](S);
  for i to nops(S0) do
    coxeter['vec2fc'](coxeter['reflect'](S0[i],v),S,'w');
    gens[i]:=op(w)
  od;
  cc:=coxeter['class_rep'](R);
  sz:=coxeter['class_size'](R);
  cc0:=coxeter['class_rep'](S0);
  sz0:=coxeter['class_size'](S0);
  res:=table('sparse');
  if f0=1 then f:=[1$nops(cc0)] else f:=f0 fi;
  for i to nops(cc0) do
    w:=[seq(gens[j],j=cc0[i])];
    w:=coxeter['class_rep'](w,R);
    if not member(w,cc,'j') then ERROR(`this cannot happen`) fi;
    res[j]:=res[j]+sz0[i]*f[i];
  od;
  j:=convert(sz,`+`)/convert(sz0,`+`);
  [seq(res[i]*j/sz[i],i=1..nops(cc))];
end:
#
`coxeter/perm_char`:=proc(J,R) `coxeter/induce`(1,[op(J)],R) end:
#
# produce an interior point v for the positive orthant defined by a
# list of independent vectors (or the base vectors of root system R).
# The coordinates must be real.
# It is important for weyl/co_rho that <v,r>=1 for each simple root r.
#
`coxeter/interior_pt`:=proc(R) local i,v,c,S,eqns;
  S:=coxeter['base'](R);
  v:=convert([seq(c[i]*S[i],i=1..nops(S))],`+`);
  v:=collect(v,indets(S));
  eqns:={seq(coxeter['iprod'](i,v)=1,i=S)};
  subs(`coxeter/linsolve`(eqns,indets(eqns)),v)
end:
#
# irr_chars(R) returns the list of irreducible characters of W(R).
# Each character is a list of values on conjugacy classes.
# The characters for A.n are ordered as in combinat[character](n+1),
#  so that the character indexed by the partition lambda appears in
#  the reverse of the order used by combinat[partition].
#
`coxeter/irr_chars`:=proc() local i,j,R,L,sys,n,irr,new,f,g;
  R:=coxeter['name_of'](args); irr:=[[1]];
  sys:=sort([op(indets(R))],`coxeter/cox_order`);
  sys:=[seq(i$degree(R,i),i=sys)];
  for R in sys do
    n:=coxeter['rank'](R);
    if type(R,'indexed') then L:=op(0,R) else L:=substring(R,1..1) fi;
    if L='A' then
      new:=convert(combinat['character'](n+1),'listlist')
    elif L='B' or L='C' then
      new:=`coxeter/irr_chars/B`(n)
    elif L='D' and n>1 then
      new:=`coxeter/irr_chars/D`(n)
    elif L='E' and n<6 then
      new:=`coxeter/irr_chars`(coxeter['base'](R))
    elif member(L,['E','F','H']) then
      new:=`coxeter/read`(R,'irr_chars','chars')
    elif R='G2' then
      new:=`coxeter/irr_chars/I2`(R,6)
    elif L='I2' and type(op(R)-1,'posint') then
      new:=`coxeter/irr_chars/I2`(R,op(R))
    else
      ERROR(cat(R,` is not a recognized root system name`))
    fi;
    irr:=[seq(seq([seq(seq(i*j,i=f),j=g)],f=irr),g=new)]
  od; irr
end:
#
`coxeter/irr_chars/B`:=proc(n) local i,j,k,ct,sp,sz,irr,cc,f,g,ch,p,q;
  for k from 0 to n do
    sp:=combinat['partition'](k);
    sz:=map(`coxeter/par2A`,sp,0);
    cc:=[seq(sz[j]*convert([seq(p[i]/2,i=sp[j])],`*`)/k!,j=1..nops(sp))];
    irr:=convert(combinat['character'](k),'listlist');
    ct[k]:=[seq(convert(zip((x,y)->x*y,f,cc),`+`),f=irr)];
  od;
  cc:=[seq(seq(seq(convert([seq(p[i],i=f),seq(q[i],i=g)],`*`),
    f=combinat['partition'](n-k)),g=combinat['partition'](k)),k=0..n)];
  sz:=[seq(seq(seq(`coxeter/par2B`(f,g,0),
    f=combinat['partition'](n-k)),g=combinat['partition'](k)),k=0..n)];
  irr:=NULL;
  for k from 0 to n do
    for g in ct[k] do
      g:=subs({seq(p[i]=p[i]-q[i],i=1..k)},g);
      for f in ct[n-k] do
        ch:=table('sparse');
        f:=expand(subs({seq(p[i]=p[i]+q[i],i=1..n-k)},f)*g);
        f:=[coeffs(f,indets(f),'sp')]; sp:=[sp];
        for i to nops(sp) do
          member(sp[i],cc,'j'); ch[j]:=2^n*n!*f[i]/sz[j]
        od;
        irr:=irr,[seq(ch[j],j=1..nops(sz))]
      od
    od
  od; [irr]
end:
#
`coxeter/irr_chars/D`:=proc(n) local i,j,jj,k,ct,sp,sz,
    irr,cc,f,g,ch,dch,p,q,evn,g0,skip,N;
  for k from 0 to n do
    sp:=combinat['partition'](k);
    sz:=map(`coxeter/par2A`,sp,0);
    cc:=[seq(sz[j]*convert([seq(p[i]/2,i=sp[j])],`*`)/k!,j=1..nops(sp))];
    irr:=convert(combinat['character'](k),'listlist');
    ct[k]:=[seq(convert(zip((x,y)->x*y,f,cc),`+`),f=irr)];
  od;
  evn:=proc(k) map(proc(x) if modp(nops(x),2)=0 then x fi end,
    combinat['partition'](k)) end;
  cc:=[seq(seq(seq(convert([seq(p[i],i=f),seq(q[i],i=g)],`*`),
    f=combinat['partition'](n-k)),g=evn(k)),k=0..n)];
  sz:=[seq(seq(seq(`coxeter/par2B`(f,g,0),
    f=combinat['partition'](n-k)),g=evn(k)),k=0..n)];
  if modp(n,2)=1 then evn:=[]
    else evn:=[seq(convert([seq(p[2*i],i=f)],`*`),
    f=combinat['partition'](n/2))] fi;
  irr:=NULL; N:=nops(sz)+nops(evn);
  for k from 0 to n/2 do
    for g0 in ct[k] do
      g:=subs({seq(p[i]=p[i]-q[i],i=1..k)},g0); skip:=false;
      for f in ct[n-k] do
        if skip then next else ch:=table('sparse') fi;
        if k=n/2 and f=g0 then skip:=true fi;
        f:=expand(subs({seq(p[i]=p[i]+q[i],i=1..n-k)},f)*g);
        f:=[coeffs(f,indets(f),'sp')]; sp:=[sp];
        for i to nops(sp) do
          if member(sp[i],cc,'j') then ch[j]:=2^n*n!*f[i]/sz[j];
            if member(sp[i],evn,'jj') then ch[jj+nops(cc)]:=ch[j] fi
          fi
        od;
        if skip then dch:=table('sparse');
          g0:=subs({seq(p[i]=p[2*i],i=1..k)},g0);
          g0:=[coeffs(g0,indets(g0),'sp')]; sp:=[sp];
          for i to nops(sp) do
            member(sp[i],cc,'j'); dch[j]:=2^n*n!*g0[i]/sz[j];
            member(sp[i],evn,'jj'); dch[nops(cc)+jj]:=-dch[j];
          od;
          irr:=irr,[seq((ch[j]+dch[j])/2,j=1..N)],
            [seq((ch[j]-dch[j])/2,j=1..N)]
        else
          irr:=irr,[seq(ch[j],j=1..N)]
        fi
      od
    od
  od; [irr]
end:
#
`coxeter/irr_chars/I2`:=proc(R,m) local cc,w,irr,i,k,st;
  cc:=coxeter['class_rep'](R); st:=2,0;
  irr:=[1$nops(cc)],[seq((-1)^nops(w),w=cc)];
  if modp(m,2)=0 then st:=st,0;
    irr:=irr,[seq((-1)^nops(w),w=subs(2=NULL,cc))],
    [seq((-1)^nops(w),w=subs(1=NULL,cc))];
  fi;
  [irr,seq([st,seq(2*cos(2*Pi*k*i/m),i=1..m/2)],k=1..(m-1)/2)]
end:
#
# length_gf(R) returns the generating function for length in W(R).
#
`coxeter/length_gf`:=proc(R) local f,d,z;
  f:=[seq(normal((1-z^d)/(1-z)),d=coxeter['degrees'](R))];
  f:=expand(convert(f,`*`));
  if nargs>1 then subs(z=args[2],f) else subs(z=q,f) fi;
end:
#
# Solve a linear system. This is a work-around, thanks to Maple
# "deprecating" the use of `solve/linear` starting with Maple8.
#
`coxeter/linsolve`:=proc() local x;
  if type(x(1)[1],'name') then # we are using Maple 7 or earlier
    readlib(`solve/linear`)(args)
  else # Maple 8 or later
    SolveTools:-Linear(args)
  fi
end:
#
# Let S be a set of equations of the form {s1=perm1,s2=perm2,...},
#  where the perms are permutations in disjoint cycle notation.
# For any word w=[i1,i2,i3,...], multperm(w,S) computes the product
#  s.i1 * s.i2 * ... in disjoint cycle form.
# The degree can be passed as an optional third argument.
# Alternatively, the second argument can be a permgroup; i.e.,
#  multperm(w,pg), where pg = permgroup(m,S), and S is as above.
#
`coxeter/multperm`:=proc(w) local S,alive,res,i,j,pi,c,cyc,g;
  if nargs>2 then
    S:=args[2]; alive:={$1..args[3]};
  elif type(args[2],'function') then
    S:=op(2,args[2]); alive:={$1..op(1,args[2])}
  else S:=args[2];
    g:={seq(cat('s',i),i=[op({op(w)})])};
    alive:=map(op,map(op,subs(S,g)))
  fi;
  res:=NULL; pi:=proc(k,x) x end;
  pi:=subsop(4=NULL,eval(pi));    # why is this necessary?
  for i in {op(w)} do
    for c in subs(S,cat('s',i)) do
      pi(i,c[nops(c)]):=c[1];
      for j to nops(c)-1 do pi(i,c[j]):=c[j+1] od
    od
  od;
  while nops(alive)>0 do
    j:=min(op(alive)); cyc:=j; i:=j;
    do
      for g in w do i:=pi(g,i) od;
      if i=j then break else cyc:=cyc,i fi;
    od;
    alive:=alive minus {cyc};
    if nops([cyc])>1 then res:=res,[cyc] fi;
  od; [res]
end:
#
# generate the label for the conjugacy class of w.
#
`coxeter/mytype`:=proc(w,L,n,R) local i,mu,sig,w0,m,pg;
  if L='A' then
    pg:={seq(cat('s',i)=[[i,i+1]],i=1..n)};
    mu:=map(nops,coxeter['multperm'](w,pg,n+1));
    [1$(n+1-convert(mu,`+`)),op(sort(mu))]
  elif L='B' or L='C' then
    pg:={s1=[[-1]],seq(cat('s',i)=[[i-1,i]],i=2..n)};
    `coxeter/mytype/bshape`(w,pg,n,'sig')
  elif L='D' then
    pg:={s1=[[-1,-2]],seq(cat('s',i)=[[i-1,i]],i=2..n)};
    mu:=`coxeter/mytype/bshape`(w,pg,n,'sig');
    if nops(mu[2])>0 or convert(map(modp,mu[1],2),`+`)>0 then sig:=NULL fi;
    mu,sig
  elif L='E' and n>5 then
    pg:=`coxeter/read`(R,'bpermrep','repn');
    mu:=`coxeter/mytype/bshape`(w,pg);
    if n<8 then mu:=(-1)^nops(w)*mu fi;
    pg:=`coxeter/read`(R,'mytype','types');
    if member(mu,pg,'i') then i else ERROR() fi;
  elif R='F4' then
    pg:=`coxeter/read`('F4','bpermrep','repn');
    mu:=`coxeter/mytype/bshape`(w,pg);
    sig:=coxeter['multperm'](subs(3=NULL,4=NULL,w),{s1=[[1,2]],s2=[[2,3]]},3);
    mu:=(-1)^nops(sig)*mu;
    pg:=`coxeter/read`('F4','mytype','types');
    if member(mu,pg,'i') then i else ERROR() fi;
  elif L='H' then
    sig:=round(2*coxeter['char_poly'](w,R,3));
    pg:=`coxeter/read`(R,'mytype','types');
    if member(sig,pg,'i') then i else ERROR() fi;
  elif L='I2' or R='G2' then
    w0:=coxeter['reduce'](w,R);
    if R='G2' then m:=6 else m:=op(R) fi;
    if w0=[] then 1
      elif modp(nops(w0),2)=1 then
        if modp(m,2)=1 then 2 else 1+w0[(nops(w0)+1)/2] fi
      else nops(w0)/2+3-modp(m,2)
    fi
  fi;
end:
#
# Determine the cycle type of a signed permutation
#
`coxeter/mytype/bshape`:=proc(w,pg,n) local alive,res,pi,c,i,j,cyc,sig;
  pi:=proc(k,x) if x<0 then -procname(k,-x) else x fi end;
  pi:=subsop(4=NULL,eval(pi));    # yuck
  for i in {op(w)} do
    for c in subs(pg,cat('s',i)) do;
      if nops(c)=1 then pi(i,-c[1]):=c[1]
        else pi(i,abs(c[2])):=c[1]; pi(i,abs(c[1])):=c[2]
      fi
    od
  od;
  alive:={$1..n}; res:=1; sig:=1;
  while nops(alive)>0 do
    j:=alive[1]; i:=j; cyc:=i;
    do
      for c in w do i:=pi(c,i) od;
      if i=j then res:=res*cat(``,nops([cyc])); break
        elif i=-j then res:=res*cat(``,-nops([cyc])); break
        else cyc:=cyc,abs(i); sig:=sig*sign(i)
      fi
    od;
    alive:=alive minus {cyc};
  od;
  if nargs=3 then subs(`1`=1,res)
    else assign(args[4],sig);
    [seq(i$degree(res,cat(``,i)),i=1..n)],
    [seq(i$degree(res,cat(``,-i)),i=1..n)]
  fi
end:
#
# name_of(X) produces the name of the Coxeter group with simple roots X,
#  or Coxeter matrix X, or Cartan matrix X. Note that Coxeter matrices
#  cannot distinguish between B/C, or G2/I2[6] (etc).
# name_of(X,'pi') does the same but also assigns to 'pi' a permutation of
#  [1,...,n] (n=rank) indicating an ordering of the roots of X that would
#  agree with the canonical ordering used for the root system of that name.
# Note that the L vector can be bogus in the case of cartan data, but we
#  only compare adjacent nodes, so it doesn't matter.
#
`coxeter/name_of`:=proc(S) local M,L,n,r,i,j;
  if type(S,'list') and nops(S)>0 then
    M:=coxeter['cox_matrix'](S); n:=nops(S);
    L:=[seq(coxeter['iprod'](r,r),r=S)];
  elif type(S,'matrix') and S[1,1]=2 then
    n:=linalg['rowdim'](S);
    L:=[seq(max(seq(-S[i,j],j=1..n)),i=1..n)];
    M:=coxeter['cox_matrix'](S);
  elif type(S,'matrix') then
    M:=S; n:=linalg['rowdim'](M); L:=[1$n];
  else
    if nargs>1 then assign(args[2],[$1..coxeter['rank'](S)]) fi;
    if S=[] then RETURN(1) else RETURN(S) fi
  fi;
  `coxeter/name_of/main`(M,L,n,args[2..nargs]);
end:

# Return list of names for direct product of algebras
`coxeter/names_of`:=proc(S) local M,L,n,r,i,j;
  if type(S,'list') and nops(S)>0 then
    M:=coxeter['cox_matrix'](S); n:=nops(S);
    L:=[seq(coxeter['iprod'](r,r),r=S)];
  elif type(S,'matrix') and S[1,1]=2 then
    n:=linalg['rowdim'](S);
    L:=[seq(max(seq(-S[i,j],j=1..n)),i=1..n)];
    M:=coxeter['cox_matrix'](S);
  elif type(S,'matrix') then
    M:=S; n:=linalg['rowdim'](M); L:=[1$n];
  else
    if nargs>1 then assign(args[2],[$1..coxeter['rank'](S)]) fi;
    if S=[] then RETURN(1) else RETURN(S) fi
  fi;
  `coxeter/name_of/main1`(M,L,n,args[2..nargs]);
end:

`coxeter/name_of/main`:=proc(M,L,n)
    convert(`coxeter/name_of/main1`(M,L,n,args[2..nargs]),`*`)
end proc:

#
# separate the connected components and toss out cases that are not
# trees or have an infinite edge or multiple edges of weight>3.
#
`coxeter/name_of/main1`:=proc(M,L,n) local X,i,j,edge,sp,res,k,l,N,EPS;
  X:=[$1..n]; EPS:=`coxeter/default`['epsilon'];
  for i to n do
    edge[i]:=NULL; sp[i]:=NULL;
    for j to i-1 do
      if M[i,j]=2 then next fi;
      if M[i,j]=infinity or X[i]=X[j] then
        ERROR(`not a finite Coxeter group`) fi;
      k:=min(X[i],X[j]); l:=max(X[i],X[j]);
      edge[k]:=edge[k],edge[l],{i,j}; sp[k]:=sp[k],sp[l];
      if M[i,j]>3 then
        if L[i]+EPS<L[j] then sp[k]:=sp[k],M[i,j],L[j]/L[i],[i,j]
        else sp[k]:=sp[k],M[i,j],L[i]/L[j],[j,i] fi
      fi;
      X:=subs(l=k,X);
    od;
    if nops([sp[X[i]]])>3 then ERROR(`not a finite Coxeter group`) fi;
  od;
  X:=sort([op({op(X)})]);
  for i in X do if sp[i]=NULL then sp[i]:=3,1,[i] fi od;
  if nargs>3 then
    res:=[seq([`coxeter/name_of/cc`({edge[i]},sp[i],true)],i=X)];
    res:=sort(res,(x,y)->`coxeter/cox_order`(x[1],y[1]));
    assign(args[4],map(x->op(op(2,x)),res));
    N:=map(x->op(1,x),res)
  else
    N:=[seq(`coxeter/name_of/cc`({edge[i]},sp[i],false),i=X)]
  fi;
    N
end:
#
# identify the connected component (and preferred ordering, if needed).
#
`coxeter/name_of/cc`:=proc(G,m,r,sp,flag) local pi,X,N,i,leafs,fork,J,n,k;
  if nops(G)<2 then
    if flag then pi:=sp else pi:=NULL fi;
    if nops(G)=0 then RETURN('A1',pi)
      elif m=4 and type(r,'rational') then RETURN('B2',pi)
      elif m=6 and type(r,'rational') then RETURN('G2',pi)
      elif m>3 then RETURN(I2[m],pi)
    fi
  fi;
  X:=sort([op(map(op,G))]); N:=table([seq(i=NULL,i=X)]);
  for k in G do N[k[1]]:=N[k[1]],k[2]; N[k[2]]:=N[k[2]],k[1] od;
  leafs:=NULL; n:=nops(X); pi:=NULL;
  for i in X do k:=nops([N[i]]);
    if k=1 then leafs:=leafs,i elif k=3 then fork:=i fi
  od;
  if nops([leafs])=2 then
    if flag then pi:=`coxeter/name_of/path`(leafs,N) fi;
    if m=3 then RETURN(cat('A',n),pi) fi;
    J:={leafs} intersect {sp[1],sp[2]};
    if nops(J)=1 then
      if flag and not pi[1]=J[1] then pi:=[seq(pi[-i],i=-n..-1)] fi;
      if m=4 and r=2 and sp[2]=J[1] then RETURN(cat('C',n),pi)
        elif m=4 then RETURN(cat('B',n),pi)
        elif m=5 and n<5 then RETURN(cat('H',n),pi)
      fi
    elif n=4 and m=4 then
      if flag and r=2 and pi[2]=sp[2] then pi:=[seq(pi[-i],i=-4..-1)] fi;
      RETURN('F4',pi)
    fi
  elif nops([leafs])=3 and m=3 then
    J:={leafs} intersect {N[fork]};
    if nops(J)>1 then
      if flag then pi:=op({leafs} minus {J[1],J[2]});
        pi:=`coxeter/name_of/path`(fork,pi,N); pi:=[J[1],J[2],op(pi)]
      fi;
      RETURN(cat('D',n),pi)
    elif nops(J)=1 and n<9 then
      X:={N[fork]} intersect {seq(N[i],i=[leafs])};
      if nops(X)>0 then
        if flag then k:=op({leafs} intersect {N[X[1]]});
          pi:=`coxeter/name_of/path`(fork,op({leafs} minus {k,J[1]}),N);
          pi:=[k,J[1],X[1],op(pi)]
        fi;
        RETURN(cat('E',n),pi)
      fi
    fi
  fi;
  ERROR(`not a finite Coxeter group`);
end:
#
# Find the unique path from vertex i to leaf j
#
`coxeter/name_of/path`:=proc(i,j,N) local path,k,old,L;
  path:=j; k:=j; old:={};
  while k<>i do
    L:=op({N[k]} minus old); old:={k}; k:=L; path:=k,path
  od; [path]
end:
#
# orbit(v,R) will return the list of all distinct images of the vector v
#  under the action of the reflection group corresponding to R.
# orbit(L,R) where L is a list or set of vectors, will start the algorithm
#  at L and work towards the anti-dominant chamber, instead of starting at
#  a point in the dominant chamber.
#
`coxeter/orbit`:=proc(v0,R) local S,coS,EPS,is_mem,orb,old,new,u,i,v;
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
  while nops(new)>0 do
    old:=new; new:=[];
    for u in old do
      for i to nops(S) do
        v:=coxeter['iprod'](coS[i],u);
        if v>EPS then v:=u-v*S[i];
          if not is_mem(v,new,EPS) then new:=[op(new),v] fi
        fi
      od
    od;
    orb:=orb,op(new);
  od; [orb]
end:
#
`coxeter/orbit/f`:=proc(v,L,EPS) local u;
  for u in L do if
    convert(map(abs,[coeffs(u-v)]),`+`)<EPS then RETURN(true)
  fi od; false
end:
#
# orbit_size(v,R) will produce the size of the orbit of the vector v under
#  the action of the reflection group corresponding to R.
# orbit_size(v,R,-1) will do the same, but for the (possibly larger) group
#  obtained by adjoining -1.
#
`coxeter/orbit_size`:=proc(v0,R) local v,S,res,EPS,J;
  S:=coxeter['base'](R);
  v:=coxeter['vec2fc'](v0,S);
  if type([op(S),v],'list'('polynom'('rational'))) then EPS:=0
    else EPS:=`coxeter/default`['epsilon']
  fi;
  J:=map(proc(x,y,z) if coxeter['iprod'](x,y)<=z then x fi end,S,v,EPS);
  res:=coxeter['size'](R)/coxeter['size'](J);
  if nargs>2 then
    v:=v-coxeter['vec2fc'](-v0,S);
    if EPS>0 and coxeter['iprod'](v,v)>EPS then res:=2*res
      elif EPS=0 and v<>0 then res:=2*res
    fi
  fi; res
end:
#
# par2A(mu)   = shortest, lex-first rep of type mu
# par2A(mu,*) = size of class indexed by mu
#
`coxeter/par2A`:=proc(mu) local sz,rep,mul,j,i;
  if nargs=1 then
    j:=0; rep:=NULL;
    for i from nops(mu) by -1 to 1 do
      rep:=rep,$j+1..j+mu[i]-1; j:=j+mu[i]
    od; [rep]
  else
    j:=0; sz:=convert(mu,`+`)!; mul:=1;
    for i in mu do
      if i=j then mul:=mul+1 else mul:=1; j:=i fi;
      sz:=sz/mul/j;
    od; sz
  fi
end:
#
# par2B(mu,nu)   = shortest, lex-first rep of type (mu,nu)
# par2B(mu,nu,*) = size of class indexed by (mu,nu)
#
`coxeter/par2B`:=proc(mu,nu) local m,n,i,k,w,w0,pi;
  if nargs=2 then
    m:=nops(nu)+1; n:=convert(nu,`+`);
    pi:=[$2..n+1]; w:=NULL; w0:=NULL;
    for i to nops(nu) do
      w0:=w0,seq(-k,k=-i..-2),$1..i;
      pi:=subsop(i=m,pi); m:=m+nu[i]-1; pi:=subs(m=i,pi)
    od;
    while n>1 do
      member(n,pi,'i'); n:=n-1;
      w:=seq(-k+1,k=-n..-i),w; pi:=subsop(i=NULL,pi)
    od;
    [w0,w,seq(i+m,i=`coxeter/par2A`(mu))]
  else
    m:=convert(mu,`+`); n:=convert(nu,`+`);
    2^(m+n-nops(mu)-nops(nu))*binomial(m+n,m)
      *`coxeter/par2A`(mu,0)*`coxeter/par2A`(nu,0)
  fi
end:
#
# par2D(mu,nu)    = shortest, lex-first rep of type (mu,nu)
#   (NULL if nops(nu) odd)
# par2D(mu,[],+1) = shortest, lex-first rep of type (mu,+), if all mu even
# par2D(mu,[],-1) = shortest, lex-first rep of type (mu,-), if all mu even
# par2D(mu,nu,0)  = class size(s) indexed by (mu,nu) (NULL if nops(nu) odd)
#
`coxeter/par2D`:=proc(mu,nu) local j,n;
  if modp(nops(nu),2)=1 then NULL
  elif nargs=2 and nops(nu)>0 then
    subs({1=2,2=1},subsop(1=NULL,seq(j*(j-1)+2=NULL,j=2..nops(nu)),
     `coxeter/par2B`(mu,nu)))
  elif nargs=2 or args[3]=-1 then
    subs(2=1,`coxeter/par2B`(mu,nu))
  elif args[3]=1 then
    `coxeter/par2B`(mu,nu)
  else
    n:=`coxeter/par2B`(mu,nu,0);
    if nu=[] and convert([seq(modp(j,2),j=mu)],`+`)=0 then n/2 else n fi
  fi
end:
#
# express permutation tau as a word in the generators of pg. If a
# stabilizer chain has been precomputed, it can be supplied as a third
# argument. FAIL is returned if tau does not belong to the group.
#
`coxeter/perm2word`:=proc(tau,pg) local res,pi,n,pt,i,w,g,sc,c;
  res:=NULL; pi:=tau; n:=1+nops(op(2,pg));
  if nargs>2 then sc:=args[3] else sc:=coxeter['stab_chain'](pg) fi;
  for g in sc while nops(pi)>0 do pt:=g[1];
    for c in pi do
      if member(pt,c,'i') then
        if nops(c)=i then pt:=c[1] else pt:=c[i+1] fi; break
      fi
    od;
    if not member(pt,g[2],'i') then RETURN(FAIL) fi;
    w:=g[3][i]; res:=op(w),res;
    pi:=coxeter['multperm']([n,seq(w[-i],i=-nops(w)..-1)],
      {op(op(2,pg)),cat('s',n)=pi});
  od;
  if pi=[] then [res] else FAIL fi;
end:
#
# perm_rep(J,R) produces the permutation representation of W(R) acting
# on cosets of the parabolic subgroup indexed by J (a list or set of
# indices in the range 1..n). If J is omitted, then the default is to
# take J to be [1,2,...n-1] (faithful iff R is irreducible).
# The format follows group[permgroup], or can be exported to gap format.
#
`coxeter/perm_rep`:=proc() local S,i,j,v,orb,pg,pi,m,find,EPS,alive,c,eq;
  m:=nargs; if args[m]='gap' then m:=m-1 fi;
  S:=coxeter['base'](args[m]);
  if m>1 then m:=args[1] else m:=[$1..(nops(S)-1)] fi;
  m:=subsop(seq(i=0,i=m),[1$nops(S)]);
  v:=convert([seq(c[i]*S[i],i=1..nops(S))],`+`);
  v:=collect(v,indets(S));
  eq:={seq(coxeter['iprod'](S[i],v)=m[i],i=1..nops(S))};
  v:=subs(`coxeter/linsolve`(eq,indets(eq)),v);
  orb:=coxeter['orbit'](v,S); pg:=NULL;
  if nops(orb)>256 and [op(I)]=[1] and args[nargs]='gap' then
    ERROR(`degree of permutation group too large to export`) fi;
  if type(S,'list'('polynom'('rational'))) then
    find:=proc(u,L) local j; member(u,L,'j'); j end;
  else
    EPS:=`coxeter/default`['epsilon'];
    find:=proc(u,L,ep) local j;
      for j to nops(L) do
        if convert(map(abs,[coeffs(u-L[j])]),`+`)<ep then RETURN(j) fi
      od; ERROR(`floating-point inaccuracies`)
    end;
  fi;
  for i to nops(S) do
    alive:=orb; m:=[$1..nops(orb)]; pi:=NULL;
    while nops(m)>1 do
      v:=coxeter['reflect'](S[i],alive[1]);
      j:=find(v,alive,EPS);
      if j>1 then pi:=pi,[m[1],m[j]] fi;
      alive:=subsop(1=NULL,j=NULL,alive);
      m:=subsop(1=NULL,j=NULL,m);
    od;
    pg:=pg,cat('s',i)=[pi]
  od;
  if args[nargs]<>'gap' then permgroup(nops(orb),{pg}) else
    `coxeter/perm_rep/gap`(nops(S),{pg},args[nargs-1])
  fi
end:
#
`coxeter/perm_rep/gap`:=proc(n,pg) local pi,i,c,lpr,R;
  if [op(I)]=[1] then lpr:=proc() printf(`%0.1200s\n`,args) end
    else lpr:=proc() lprint(args) end fi;
  lpr(` `);
  for i to n do
    pi:=subs(pg,cat('s',i));
    if nops(pi)=0 then c:=[`()`] else
      c:=map(y->convert(y,'string'),pi);
      c:=map(y->cat(`(`,substring(y,2..length(y)-1),`)`),c);
    fi;
    lpr(cat('s',i,`:=`,op(c),`;`));
  od;
  if n=0 then lpr(`G:=Group(());`)
    else lpr(cat(`G:=Group(s1`,seq(op([`,s`,i]),i=2..n),`);`))
  fi;
  R:=coxeter['name_of'](args[3]);
  lpr(cat(`G.name:="`,convert(R,'string'),`";`)); lpr(` `);
end:
#
# pos_roots(R) returns the list of all positive roots.
#
`coxeter/pos_roots`:=proc(R) local S;
  S:=coxeter['base'](R);
  map(x->-x,coxeter['orbit'](map(x->-x,S),S))
end:
#
# produce a presentation of W(R) in terms of generators and relations,
# suitable for use with the group package, or for export to gap.
#
`coxeter/presentation`:=proc(R) local n,M,i,j,gens,rels,lpr;
  M:=coxeter['cox_matrix'](R);
  n:=coxeter['rank'](M);
  if [op(I)]=[1] then lpr:=proc() printf(`%0.512s\n`,args) end
    else lpr:=proc() lprint(args) end fi;
  if not member('gap',[args]) then
    gens:={seq(cat('s',i),i=1..n)};
    rels:={seq([cat('s',i),cat('s',i)],i=1..n),
      seq(seq(map(op,[[cat('s',j),cat('s',i)]$M[i,j]]),j=1..i-1),i=2..n)};
    grelgroup(gens,rels)
  elif n=0 then
    lpr(` `); lpr(`F:=FreeGroup(0);`);
    lpr(`G:=F/[];`); lpr(`G.name:="1";`); lpr(` `)
  else lpr(` `);
    lpr(cat(`F:=FreeGroup("s1`,seq(op([`","s`,i]),i=2..n),`");`));
    lpr(cat(`G:=F/[F.1^2`,seq(cat(`,F.`,i,`^2`),i=2..n),
      seq(seq(cat(`,(F.`,j,`*F.`,i,`)^`,M[i,j]),j=1..i-1),i=2..n),`];`));
    lpr(cat(`G.name:="`,convert(coxeter['name_of'](R),'string'),`";`));
    lpr(` `);
  fi
end:
#
# rank(R) returns the rank of root system R.
#
`coxeter/rank`:=proc(X) local dd;
  if type(X,'list') then nops(X)
    elif type(X,'matrix') then linalg['rowdim'](X)
    else dd:=[`0`,`1`,`2`,`3`,`4`,`5`,`6`,`7`,`8`,`9`];
      convert(map(`coxeter/rank/number`,[op(indets(X))],X,dd),`+`)
  fi
end:
#
`coxeter/rank/number`:=proc(R,X,dd) local b,n,i,j;
  n:=0; if type(R,'indexed') then b:=op(0,R) else b:=R fi;
  for i from 2 to length(b) do
    member(substring(b,i..i),dd,'j'); n:=10*n+j-1;
  od;
  if type(n,'posint') then n*degree(X,R)
    else ERROR(cat(R,` is not a recognized root system name`))
  fi
end:
#
# reduce([i_1,...,i_r],R) returns a reduced expression representing the
#  same element of the Coxeter group as [i_1,...,i_r].
# If a third argument is present, the procedure returns a reduced *subword*
#  using the exchange property (a quadratic algorithm).
#
`coxeter/reduce`:=proc(w,R) local S,w0,i,v;
  S:=coxeter['base'](R);
  if nargs<3 then
    v:=coxeter['interior_pt'](S);
    v:=coxeter['reflect'](seq(S[i],i=w),v);
    coxeter['vec2fc'](v,S,'w0'); w0
  elif type(S,'list'('polynom'('rational'))) then
    `coxeter/reduce/subwd`(w,S,0)
  else
    `coxeter/reduce/subwd`(w,S,`coxeter/default`['epsilon'])
  fi
end:
#
`coxeter/reduce/subwd`:=proc(w,S,EPS) local w0,r,i;
  if nops(w)>1 then r:=S[w[1]] else RETURN(w) fi;
  w0:=procname(subsop(1=NULL,w),S,EPS);
  for i to nops(w0) do
    if EPS=0 then
      if S[w0[i]]=r then RETURN(subsop(i=NULL,w0)) fi
    elif convert(map(abs,[coeffs(S[w0[i]]-r)]),`+`)<EPS then
      RETURN(subsop(i=NULL,w0))
    fi;
    r:=coxeter['reflect'](S[w0[i]],r)
  od;
  [w[1],op(w0)];
end:
#
# reflect(r,v) applies the reflection through r to vector v.
# reflect(r1,r2,...,v) iterates the application of the reflections
#  r1,r2,... to vector v.
#
# iprod(u,v) computes the standard inner product of two vectors,
# expressed as a linear combo of the standard basis e1,e2,e3,...
#
`coxeter/reflect`:=proc() local i,v,r;
  v:=args[nargs];
  for i from nargs-1 by -1 to 1 do
    r:=args[i];
    v:=v-(2*`coxeter/iprod`(r,v)/`coxeter/iprod`(r,r))*r;
  od; v
end:
#
`coxeter/iprod`:=proc(u,v) local x;
  convert([seq(coeff(u,x)*coeff(v,x),x=indets(u))],`+`)
end:
#
# Restrict a class function f on W(R) to some reflection subgroup.
# J can be either a sublist of [1,...,n] indicating a parabolic subgroup,
# or a base for some (not necessarily parabolic) reflection subgroup.
#
`coxeter/restrict`:=proc(f,J,R) local S,S0,w,gens,cc,cc0,res,v,i,j;
  S:=coxeter['base'](R); S0:=J;
  if type(S0,'list'('integer')) then S0:=[seq(S[i],i=S0)] fi;
  v:=coxeter['interior_pt'](S);
  for i to nops(S0) do
    coxeter['vec2fc'](coxeter['reflect'](S0[i],v),S,'w');
    gens[i]:=op(w)
  od;
  cc:=coxeter['class_rep'](R);
  cc0:=coxeter['class_rep'](S0);
  for i to nops(cc0) do
    w:=[seq(gens[j],j=cc0[i])];
    w:=coxeter['class_rep'](w,R);
    if not member(w,cc,'j') then ERROR(`this cannot happen`) fi;
    res[i]:=f[j]
  od;
  [seq(res[i],i=1..nops(cc0))]
end:
#
# root_coords(v,R) returns the coordinates of v with respect to base(R).
# If v is not in the linear span of base(R), either an error will occur
# or garbage will be produced.
#
`coxeter/root_coords`:=proc(v,R) local S,i,c,u,vars,eq;
  S:=coxeter['base'](R);
  vars:=seq(c[i],i=1..nops(S));
  u:=convert([seq(c[i]*S[i],i=1..nops(S))],`+`)-v;
  if type(u,'polynom'('rational')) or nops(S)=nops(indets(S)) then
    eq:={coeffs(expand(u),indets(S))}
  else
    u:=collect(u,indets(S));
    eq:={seq(coxeter['iprod'](i,u),i=S)}
  fi;
  subs(`coxeter/linsolve`(eq,{vars}),[vars])
end:
#
# Stabilizer chain for a permutation rep pg of a reflection group W.
# Works under the assumption that the representation is isomorphic to
#  the action of W on some proper parabolic subgroup.
# For a chain W > W1 > W2 > ..., the output is [[b1,O1,X1],[b2,O2,X2],...],
#  where W2 is the stabilizer of b2 in W1, O2 is the list of points in
#  the W1-orbit of b2, X2 is a list of coset reps for W2/W1, ....
# In most cases, if the generators [1,...,n] of W are listed in the standard
#  order, one can take W.i to be generated by [1,...,n-i], so the algorithm
#  tries to do this first before trying any other possibilities.
#
`coxeter/stab_chain`:=proc(pg) local n,i,pi,alive,X,x,b,new,j,gens,reps,ct,k;
  n:=nops(op(2,pg)); alive:=[$1..n]; gens:=NULL;
  pi:=subs(op(2,pg),[seq(cat('s',i),i=1..n)]);
  pi:=[seq(map(proc(y) if nops(y)=2 then op(y) else
    ERROR(`improper permutation representation`) fi end,pi[i]),i=1..n)];
  while n>0 do;
    X:={seq(op(pi[alive[i]]),i=1..n-1)};
    b:={op(pi[alive[n]])} minus X;
    if nops(b)>0 then
      b:=min(op(b)); new:=subsop(n=NULL,alive);
    else
      new:=-1; ct:=table([seq(x=NULL,x=X)]);
      for i in alive do
        for x in X minus {op(pi[i])} do ct[x]:=ct[x],i od;
      od;
      for x in X do k:=nops({ct[x]});
        if k>new then new:=k; b:=x fi
      od;
      new:=[ct[b]]
    fi;
    X:=[b]; reps:=[[]];
    for j while j<=nops(X) do
      for i in alive do
        if member(X[j],pi[i],'k') then
          x:=pi[i][k-1+2*modp(k,2)];     # we know pi[i]^2=1
          if not member(x,X) then X:=[op(X),x];
            reps:=[op(reps),[op(reps[j]),i]]
          fi
        fi
      od
    od;
    gens:=gens,[b,X,reps]; alive:=new; n:=nops(alive);
    userinfo(2,'stab_chain',`Orbit size:`,nops(X),`Stabilizer:`,alive);
  od; [gens]
end:
#
# vec2fc(v,R) produces the unique vector v0 in the orbit of v that belongs
#  to the (closure of the) fundamental chamber of R.
# vec2fc(v,R,'w') does the same thing, but also assigns to w the shortest
#  element of W(R) with the property that w.v0 =v.
# longest_elt(R) produces the longest element of W(R).
#
`coxeter/vec2fc`:=proc() local S,EPS,v,w,i;
  S:=coxeter['base'](args[2]); v:=args[1]; w:=NULL;
  if type([op(S),v],'list'('polynom'('rational'))) then EPS:=0
    else EPS:=-`coxeter/default`['epsilon']
  fi;
  do
    for i to nops(S) do
      if coxeter['iprod'](S[i],v)<EPS then w:=w,i;
        v:=coxeter['reflect'](S[i],v); break
      fi
    od;
    if i>nops(S) then break fi;
  od;
  if nargs>2 then assign(args[3],[w]) fi; v
end:
#
`coxeter/longest_elt`:=proc(R) local w,S;
  S:=coxeter['base'](R);
  `coxeter/vec2fc`(-coxeter['interior_pt'](S),S,'w'); w;
end:
#
# branch(v,J,R) determines the irreducible decomposition when an irrep
#  of LieAlg(R) of highest weight v is restricted to the (reductive)
#  subalgebra generated by a Cartan subalgebra and the root spaces
#  supported on simple roots indexed by J.
# branch(v,J,R,'klimyk') uses the Brauer-Klimyk Rule (the default).
# branch(v,J,R,'mchain') uses the Minuscule Chain algorithm.
#  In this case, v can be a dominant weight or a (quasi-)minuscule chain.
#
`weyl/branch`:=proc(v,J,R) local n;
  if nargs<4 or args[4]='klimyk' then
    n:=coxeter['rank'](R);
    `weyl/klimyk`(v,0,{$1..n} minus {op(J)},R)
  elif args[4]='mchain' then
    `weyl/mchain`(v,0,{op(J)},R)
  else
    ERROR(cat(`unknown option `,args[4]))
  fi
end:
#
# Find a representation of dominant weight v as a sum of minuscule and
# quasi-minuscule weights with dominant partial sums (=minuscule chain).
# Try to use minuscules where possible, but not at the expense
# of creating a longer representation.
#
# The chain is built by concatenating minimum-length chains for the
# fundamental weights. The validity of the algorithm rests on the
# empirical observation that a minimum-length chain for a non-minuscule
# fund.wt. w[i] is obtained by appending a (quasi-)minuscule weight to
# a minimal chain for some previous w[j].
#
`weyl/findrep`:=proc(v,S) local n,w,qm,vc,J,L,new,old,rep,i,u;
  w:=weyl['weights'](S); n:=nops(w);
  vc:=weyl['weight_coords'](v,S);
  qm:=weyl['minuscule'](S,'quasi',w);
  new:=[0,0]; rep[0]:=NULL; L:={$1..n};
  J:=map(proc(x,y) if y[x]>0 then x fi end,L,vc);
  while nops(J)>0 do
    old:=[new]; new:=NULL;
    if old=[] then ERROR(`this should not happen`) fi;
    for i in L do
      for u in old do
        if member(coxeter['vec2fc'](w[i]-u[1],S),qm) then
          rep[i]:=rep[u[2]],w[i]-u[1];
          new:=new,[w[i],i]; J:=subs(i=NULL,J); break
        fi;
      od;
      if nops(J)=0 then break fi;
    od;
    L:=L minus map(x->op(2,x),{new});
  od;
  map(op,[seq([rep[i]]$vc[i],i=1..n)]);
end:
#
# Given a polynomial expression f in the variables M[..] and X[..],
# flatten it into a linear combo of M[..]'s and X[..]'s expeditiously.
# Heuristics: (1) Typical product expressions will have the smaller factors
# listed first, so parse from right to left. (2) Use the monomial*monomial
# algorithm only when the entire expression uses only M[..]-variables.
#
`weyl/flatten`:=proc(f,S,coS,wt,Nghbr) local v,g,res,i;
  v:=indets(f,'indexed');
  v:=map(proc(x) if op(0,x)='M' or op(0,x)='X' then x fi end,v);
  if not type(f,'polynom'('anything',v)) then
    ERROR(`invalid expression`) fi;
  if degree(f,v)<2 then res:=collect(f,v,'distributed')
    elif type(f,`+`) then
      g:=map(procname,args); v:=indets(g,'indexed');
      v:=map(proc(x) if op(0,x)='M' or op(0,x)='X' then x fi end,v);
      res:=collect(g,v,'distributed')
    elif type(f,`*`) then
      g:=map(procname,[op(f)],args[2..5]);
      g:=[seq(g[-i],i=-nops(g)..-1)];
      if map(x->op(0,x),v)={'M'} then res:=g[1];
        for i from 2 to nops(g) do
          res:=`weyl/flatten/multM`(g[i],res,S,coS,wt,Nghbr)
        od
      else
        for i to nops(g) do
          v:=map(x->op(0,x),indets(g[i],'indexed'));
          if member('X',v) and not member('M',v) then
            g:=subsop(1=g[i],i=g[1],g); break fi;
        od;
        if i>nops(g) then g:=[X[0$nops(S)],op(g)] fi;
        res:=g[1];
        for i from 2 to nops(g) do
          res:=`weyl/flatten/multX`(g[i],res,S,coS,wt,Nghbr)
        od
      fi
    elif type(f,`^`) and type(op(2,f),'posint') then
      g:=procname(op(1,f),args[2..5]); res:=g;
      if map(x->op(0,x),v)={'M'} then
        for i from 2 to op(2,f) do
          res:=`weyl/flatten/multM`(g,res,S,coS,wt,Nghbr)
        od
      else
        g:=weyl['toM'](g,S,coS,wt,'flat');
        if member('M',map(x->op(0,x),v)) then
          res:=`weyl/flatten/multX`(res,X[0$nops(S)],S,coS,wt,Nghbr)
        fi;
        for i from 2 to op(2,f) do
          res:=`weyl/flatten/multX`(g,res,S,coS,wt,Nghbr)
        od
      fi
    else ERROR(`invalid expression`)
  fi; res
end:
#
# The monomial*monomial algorithm.
# Given expressions f,g that are linear combos of M[..]'s (orbit sums),
# Convert f*g to a linear combo of the M[..]-basis. For max efficiency,
# the "smaller" of the two factors (measured by orbit sizes) should be f.
#
`weyl/flatten/multM`:=proc(f,g,S,coS,wt,Nghbr) local res,v,vf,vg,u,w;
  vf:=`weyl/flatten/vars`(f,wt,'M');
  vg:=`weyl/flatten/vars`(g,wt,'M');
  res:=0;
  for u in vf do
    for w in vg do
      res:=res+u[2]*w[2]*`weyl/orbit_sum`(u[3],w[3],S,coS,Nghbr,'symm');
      v:=indets(res,'indexed'); # collect at each step--avoid swell
      v:=map(proc(x) if op(0,x)='M' then x fi end,v);
      res:=collect(res,v,'distributed')
    od
  od; res
end:
#
# The monomial*character algorithm.
# Given an expression f0 that is a linear combo of M[..]'s and X[..]'s,
# and a g that is a linear combo of X[..]'s, convert f0 to the M[..]
# basis, and use Brauer-Klimyk to expand f0*g in the X[..] basis.
#
`weyl/flatten/multX`:=proc(f0,g,S,coS,wt,Nghbr) local res,f,r0,vf,vg,v,u,w;
  f:=weyl['toM'](f0,S,coS,wt,'flat');
  vf:=`weyl/flatten/vars`(f,wt,'M');
  vg:=`weyl/flatten/vars`(g,wt,'X');
  res:=0; r0:=convert(wt,`+`);
  for u in vf do
    for w in vg do
      res:=res+u[2]*w[2]*`weyl/orbit_sum`(u[3],w[3]+r0,S,coS,Nghbr,'skew');
      v:=indets(res,'indexed'); # collect at each step--avoid swell
      v:=map(proc(x) if op(0,x)='X' then x fi end,v);
      res:=collect(res,v,'distributed')
    od
  od; res
end:
#
# Assuming f0 is linear and collected with respect to a set of variables
# of the form L[..], where [..] represents fundamental weight coordinates,
# pull apart f0 into a list, with each operand having the form [V,c,v],
# where V=L[..] is the variable name, c is the coefficient of V, and v is
# the vector form of the dominant weight corresponding to V.
# Special attention is needed for the constant term.
#
`weyl/flatten/vars`:=proc(f0,wt,L) local v,f,u,ct;
  v:=indets(f0,'indexed');
  v:=map(proc(x,y) if op(0,x)=y then x fi end,v,L);
  ct:=subs({seq(u=0,u=v)},f0); f:=f0-ct+ct*L[0$nops(wt)];
  if ct<>0 then v:={op(v),L[0$nops(wt)]} fi;
  [seq([u,coeff(f,u),convert(zip((x,y)->x*y,[op(u)],wt),`+`)],u=v)]
end:
#
# Compute tensor product multiplicities and branching multiplicities via
#  the Brauer-Klimyk Rule. The smaller of the two irreps should be u.
# Calling sequence:  `weyl/klimyk`(u,v,J,R), where
#  J is the set of nodes to delete ( = {} for tensor products),
#  v is the second dominant weight ( = 0 for branching).
# Note that branching to arbitrary root subsystems is possible with
#  only slight modifications.
#
`weyl/klimyk`:=proc() local res,S,coS,v0,wt,wm,a,u,m,Nghbr;
  res:=0; S:=coxeter['base'](args[4]);
  a:=ilcm(op(map(denom,S))); S:=map((x,y)->y*x,S,a);
  coS:=coxeter['co_base'](S);
  wt:=weyl['weights'](S);
  v0:=a*args[2]+convert(wt,`+`);
  wm:=weyl['weight_mults'](a*args[1],S);
  Nghbr:=`weyl/neighbor`(S,args[3],0);
  if type(wm,`+`) then wm:=[op(wm)] else wm:=[wm] fi;
  for a in wm do
    u:=op(indets(a)); m:=coeff(a,u);
    u:=convert(zip((x,y)->x*y,[op(u)],wt),`+`);
    res:=res+m*`weyl/orbit_sum`(u,v0,S,coS,Nghbr,'skew');
  od; res
end:
#
# Determine the list of locally short dominant positive roots, including
#  the "exceptional" root in G2.
# A positive root r has this property if it is both short and dominant
#  relative to the root subsystem generated by the simple roots that
#  occur in the support of r. In G2, the sum of the two simple roots does
#  not have this property, but is part of the list by convention.
# If u and v are dominant weights such that u < v (i.e., v - u is a sum
#  of positive roots), then there must be a locally short dominant root r
#  such that v-r is dominant and u <= v-r.
#
# Reference: The partial order of dominant weights,
#  Adv. Math 136 (1998), 340--364.
#
`weyl/lsdpr`:=proc(R) local S,S0;
  S:=coxeter['base'](R);
  if nops(S)>0 then S0:=op(subsop(1=NULL,S));
    [`weyl/lsdpr/sub`(S[1],[S[1]],{S0}),op(procname([S0]))]
  else [] fi
end:
#
# Find all lsdpr's that include S1 in their support, where r is the unique
# lsdpr with support S1, and S2 is the complementary set of simple roots.
#
`weyl/lsdpr/sub`:=proc(r,S1,S2) local r0,c,s,S0;
  for s in S2 do
    c:=coxeter['iprod'](r,s);
    if c<0 then
      c:=2*c/coxeter['iprod'](s,s);
      if c<-1 then r0:=coxeter['vec2fc'](s,[op(S1),s])
        else r0:=coxeter['vec2fc'](r,[s,op(S1)]) fi;
      S0:=S2 minus {s};
      RETURN(procname(r,S1,S0),procname(r0,[s,op(S1)],S0))
    fi
  od;
  if nops(S1)=2 and r<>S1[1]+S1[2] then r,S1[1]+S1[2] else r fi
end:
#
# Compute tensor product multiplicities and branching multiplicities via
#  minuscule chain algorithm. The smaller of the two irreps should be u.
# Calling sequence:  `weyl/mchain`(u,v,J,R), where
#  u can either be a dominant weight or a chain
#  J is the set of nodes to keep ( = {1,...,n} for tensor products),
#  v is the second dominant weight ( = 0 for branching).
#
`weyl/mchain`:=proc(u,v,J,R) local S,coS,rep,chld,nc,c,i,rs,res,ct;
  S:=coxeter['base'](R); coS:=coxeter['co_base'](S);
  if type(u,'list') then rep:=u else rep:=`weyl/findrep`(u,S) fi;
  chld:=map(`weyl/mchain/downset`,rep,S,coS);
  nc:=map(nops,chld);
  rs[0]:=[map(coxeter['iprod'],coS,v),[0$nops(S)]];
  res:=0; c[1]:=0; i:=1; ct:=0;
  while i>0 do
    if i>nops(rep) then i:=i-1; res:=res+X[op(rs[i][1])]
    elif c[i]=nc[i] then i:=i-1
    else c[i]:=c[i]+1;
      rs[i]:=`weyl/mchain/insert`(c[i],chld[i],rs[i-1],J);
      if rs[i]=NULL then ct:=ct+1 else i:=i+1; c[i]:=0 fi;
    fi
  od;
  userinfo(2,'mchain',cat(`number of dead ends=`,ct)); res
end:
#
# Given a rise-sequence rs, generate the rise-sequence obtained
#  by inserting the j-th term of chld. Reject if it does not terminate
#  at the top element of chld, or if the total weight is not dominant.
# Possible inefficiency: this does not use seq/map to build the answer.
#
`weyl/mchain/insert`:=proc(j,chld,rs,J) local wc,u,new,r,i,eps;
  u:=chld[j]; wc:=zip((x,y)->x+y,u[1],rs[1]);
  if min(1,seq(wc[i],i=J))<0 then RETURN() fi;
  i:=nops(rs); new:=NULL;
  while i>1 do
    eps:=zip((x,y)->max(0,x+y),rs[i],u[1]);
    if u[3]>0 then eps:=subsop(u[3]=eps[u[3]]+1,eps) fi;
    new:=eps,new;
    for r in u[2] do
      if rs[i][r[2]]+u[1][r[2]]<0 then
        if r[1]=0 then RETURN() fi;
        u:=chld[r[1]]; i:=i+1; break
      fi
    od;
    i:=i-1;
  od;
  if u[3]<0 then if chld[j][3]>0 then wc:=rs[1] fi; [wc,new] fi;
end:
#
# Given a **nonzero** (quasi-)minuscule v0, generate a list of objects
# below v0 in the corresponding thin system. The list has the format
#   [weight-coords, edge-list, flag], where
# flag=-1 indicates that this is the top object, of weight v0<>0.
# flag=j>0 indicates that this is the object 0_j.
#   In this case, "weight-coords" = [0...-1...0].
# flag=0 for all other objects.
# Edge-lists are pairs [i,j], where i is the list position of the object
#  obtained by applying E_j. Set i=0 if applying E_j would yield an
#  object that cannot reach v0.
#
`weyl/mchain/downset`:=proc(v0,S,coS)
  local new,J,u,v,wt,i,j,k,old,mu,eps,nu,p,q;
  J:=[$1..nops(S)]; new:=[v0]; q:=0; eps[1]:=NULL;
  while nops(new)>0 do
    old:=new; new:=[]; p:=q; q:=q+nops(old);
    if old=[0] then k:=nops([eps[p+1]]);
      for u in [eps[p+1]] do
        j:=u[2]; nu:=-S[j]; new:=[op(new),nu];
        mu[q]:=[subsop(j=-1,[0$nops(S)]),[u],j];
        eps[q+k]:=[q,j]; q:=q+1;
      od; q:=q-1;
    else
      for i from p+1 to q do
        v:=old[i-p]; wt:=map(coxeter['iprod'],coS,v);
        nu:=subsop(seq(j[2]=0,j=[eps[i]]),wt);
        nu:=map(proc(x,y) if y[x]<0 then [0,x] fi end,J,nu);
        mu[i]:=[wt,[op(nu),eps[i]],0];
        for j to nops(S) do
          if wt[j]>0 then u:=v-S[j];
            if not member(u,new,'k') then
              new:=[op(new),u]; k:=nops(new); eps[q+k]:=NULL fi;
            eps[q+k]:=eps[q+k],[i,j];
          fi
        od
      od
    fi;
  od;
  [subsop(3=-1,mu[1]),seq(mu[i],i=2..q)];
end:
#
# minuscule(R) returns a list of all minuscule fundamental weights.
# minuscule(R,'quasi') returns a list of all minuscule fundamental and
#  quasi-minuscule weights, with the minuscules listed first.
# The fundamental weights can be passed along in args[3].
#
`weyl/minuscule`:=proc(R) local c,i,j,r,hcr,res,w,S,L;
  S:=coxeter['base'](R); c:=[seq(i,i=1..nops(S))];
  L:=[seq(coxeter['iprod'](r,r),r=S)];
  for i to nops(S) do
    for j to i-1 do
      r:=coxeter['iprod'](S[i],S[j]);
      if r>=0 then next
        elif L[c[i]]<L[c[j]] then c:=subs(c[j]=c[i],c)
        else c:=subs(c[i]=c[j],c)
      fi
    od
  od;
  if nargs>2 then w:=args[3] else w:=weyl['weights'](S) fi;
  hcr:=table([seq(i=coxeter['vec2fc'](2*S[i]/L[i],S),i={op(c)})]);
  res:=[seq([w[i],coxeter['iprod'](hcr[c[i]],w[i])],i=1..nops(S))];
  res:=map(proc(x) if x[2]=1 then x[1] fi end,res);
  if nargs=1 then RETURN(res) fi;
  hcr:=[seq(hcr[i],i=sort([op({op(c)})]))];
  [op(res),seq(2*r/coxeter['iprod'](r,r),r=hcr)]
end:
#
# Build the "neighbor table" for a root system R of rank n:
#  Entry i is a list of neighbors j of the i-th simple root with j>i.
#  For all such pairs, entry (i,j) is the list of roots obtained by
#  applying the j-th reflection to the simple roots i,i+1,...,j-1.
#  Entry n+1 is an empty list. Entry 0 is [n,n-1,...,2,1].
#  Entry -i is the list of neighbors and 1,...,i-1, sorted in reverse.
# Assuming J is a subset of {1,...,n}, then each member of J is purged
#  from the entries 0,-1,...,-n (used for moving into dominant chambers
#  of corresponding root subsystems).
# Note: neg = a small negative float (or 0, for the xtal case).
# Reference: Sec. 4 in MSJ Memoirs Vol. 11, Math Soc. Japan, 2001,
#  <http://www.math.lsa.umich.edu/~jrs/papers/carswc.ps.gz>.
#
`weyl/neighbor`:=proc(R,J,neg) local S,i,j,IP,sJ;
  S:=coxeter['base'](R);
  sJ:={seq(j=NULL,j=J)};
  for i to nops(S) do
    IP[i]:=NULL; IP[-i]:=seq(-j,j=1-i..-1);
    for j from i+1 to nops(S) do
      if coxeter['iprod'](S[i],S[j])<neg then
        IP[i]:=IP[i],j; IP[-i]:=j,IP[-i];
        IP[i,j]:=map((x,y)->coxeter['reflect'](y,x),[op(i..j-1,S)],S[j]);
      fi
    od;
    IP[i]:=[IP[i]]; IP[-i]:=subs(sJ,[IP[-i]]);
  od;
  IP[0]:=subs(sJ,[seq(-j,j=-nops(S)..-1)]);
  IP[nops(S)+1]:=[]; op(IP);
end:
#
# `weyl/orbit_sum`(u0,v0+rho,S,coS,Nghbr,'skew')
#  computes the irreducible character expansion of the product of the
#  orbit sum for u0 and the Weyl character indexed by v0 using the
#  Brauer-Klimyk Rule. The output is a linear combo of variables X[..],
#  where [..] are fund.wt.coords.
# Here, S,coS = simple roots and co-roots, Nghbr = the neighbor table.
# Similarly, `weyl/orbit_sum`(u0,v0,S,coS,Nghbr,'symm') computes the
#  expansion of the product of the orbit sums corresponding to u0 and v0
#  as a linear combination of orbit sums. The output is a linear combo
#  of variables M[..], where [..] are fund.wt.coords.
# Note that the orbit of u0 is efficiently traversed.
#
`weyl/orbit_sum`:=proc(u0,v0,S,coS,Nghbr)
  local u,res,i,j,f,g,r,rnk,Cstack,cc;
  rnk:=0; Cstack:=table();
  g:=cat(`weyl/orbit_sum/`,args[6]);  # args[6]='skew' or 'symm'
  i:=nops(S)+1; u:=u0; res:=0;
  do
    res:=res+g(u+v0,S,coS,Nghbr);
    cc:=[seq(`weyl/orbit_sum/ref`(u,j,S,coS),j=1..i-1)];
    for j in Nghbr[i] do f:=1;
      for r in Nghbr[i,j] do
        if coxeter['iprod'](r,u)<0 then f:=0; break fi
      od;
      if f=1 then cc:=[op(cc),`weyl/orbit_sum/ref`(u,j,S,coS)] fi;
    od;
    if nops(cc)>0 then rnk:=rnk+1 else
      do cc:=Cstack[rnk];
        if nops(cc)=0 then rnk:=rnk-1 else break fi;
      od;
      if rnk=0 then break fi
    fi;
    Cstack[rnk]:=subsop(1=NULL,cc);
    i:=cc[1][1]; u:=cc[1][2];
  od;
  if args[6]='symm' then f:=min(coeffs(res)); res/f else res fi
end:
#
# Given v0, return sgn(w)*X[..] or 0, where [..] are the fund.wt.coords
# of the unique dominant v such that w.(v+rho) = v0 (if v exists).
# This is faster than 'vec2fc', since we exploit neighborhood data.
#
`weyl/orbit_sum/skew`:=proc(v0,S,coS,N) local v,l,i,j,c,ndom;
  v:=v0; i:=0; ndom:=true;
  for l while ndom do
    ndom:=false;
    for j in N[-i] do
      c:=coxeter['iprod'](coS[j],v);
      if c<0 then v:=v-c*S[j]; ndom:=true; break
        elif c=0 then RETURN(0) fi;
    od;
    i:=j
  od;
  (-1)^l*X[seq(coxeter['iprod'](j,v)-1,j=coS)]
end:
#
# Given v0, return the expression c*M[..], where [..]=fund.wt.coords
# of the unique dominant weight in the W-orbit of v0, and c is the order
# of the W-stabilizer of v0.
#
`weyl/orbit_sum/symm`:=proc(v0,S,coS,N) local v,i,j,c,ndom,S0;
  v:=v0; i:=0; ndom:=true;
  while ndom do
    ndom:=false;
    for j in N[-i] do
      c:=coxeter['iprod'](coS[j],v);
      if c<0 then v:=v-c*S[j]; ndom:=true; break fi
    od;
    i:=j
  od;
  v:=[seq(coxeter['iprod'](j,v),j=coS)];
  S0:=map(proc(x,y,z) if y[x]=0 then z[x] fi end,[$1..nops(S)],v,S);
  coxeter['size'](S0)*M[op(v)]
end:
#
`weyl/orbit_sum/ref`:=proc(v,i,S,coS) local c;
  c:=coxeter['iprod'](coS[i],v);
  if c>0 then [i,v-c*S[i]] fi
end:
#
# Use the "qtensor" (quick and dirty) algorithm for decomposing the
#  tensor product of irreps with highest weights u and v.
# Note the secondary uses of X here as a global variable
# Reference: <http://www.math.lsa.umich.edu/~jrs/papers/carswc.ps.gz>,
#  MSJ Memoirs Vol. 11 pp. 1--38, Math. Soc. Japan, Tokyo, 2001.
#
`weyl/qtensor`:=proc(u,v,R)
  local coS,r,pr,copr,cow,r0,wts,i,c,vars,ht,wl,h0,wc,nm;
  coS:=coxeter['co_base'](R);
  copr:=coxeter['pos_roots'](coS);
  pr:=[seq(2*r/coxeter['iprod'](r,r),r=copr)];
  cow:=weyl['weights'](coS); r0:=convert(pr,`+`)/2;
  r:=coxeter['vec2fc'](u-coxeter['vec2fc'](-v,coS),coS);
  wts:=weyl['weight_sys'](u+v,R,'wc',r,cow);
  vars:=[seq(c[i],i=1..nops(wts))];
  wts:=map((x,y)->x+y,[u,v,op(wts)],r0);
  nm:=[seq(convert(map((x,y)->X[coxeter['iprod'](x,y)],copr,r),`*`),r=wts)];
  ht:=map(coxeter['iprod'],cow,r0);
  h0:=min(op(ht)); wl:=[[nops(coS),0,0]];
  for i while nops(indets(vars))>0 do
    if i>nops(wl) then wl:=`weyl/qtensor/morecow`(wl,cow,ht,h0) fi;
    userinfo(2,'qtensor',nops(indets(vars)),`using catalyst`,wl[i][3]);
    vars:=`weyl/qtensor/spec`(wts,wl[i][3],u+v+r0,nm,pr,cow,coS,vars);
  od;
  convert(zip((x,y)->x*X[op(y)],vars,wc),`+`);
end:
#
# Given a (structured) list of all co-weights up to some height ht,
# build a larger list of all co-weights up to height ht+h0.
#
`weyl/qtensor/morecow`:=proc(wl,cow,ht,h0) local new,m,u,h,i;
  new:=NULL; m:=wl[nops(wl)][2];
  for u in wl do
    for i to u[1] do h:=ht[i]+u[2];
      if h<=m+h0 and h>m then new:=new,[i,h,u[3]+cow[i]] fi
    od
  od;
  new:=sort([new],(x,y)->evalb(x[2]<=y[2]));
  [op(wl),op(new)];
end:
#
# Notice that specializations need not be symmetric with respect to
# replacing q->1/q, unless u0 is self-dual. However, if we only
# build equations up to the middle degree, the upper half will get
# built when we use the dual of u0.
#
`weyl/qtensor/spec`:=proc(wts,u0,v1,nm,pr,cow,coS)
  local u1,ev,i,j,q,N,k,dn,vars;
  N:=nops(wts); u1:=u0+convert(cow,`+`); vars:=args[8];
  dn:=convert([seq(X[coxeter['iprod'](j,u1)],j=pr)],`*`);
  dn:=map((x,y)->x/y,nm,dn);
  ev:=`weyl/qtensor/eval`(wts,u0,cow,coS,q);
  ev:=`weyl/qtensor/cyclosimp`(ev,dn,q);
  ev:=[[0,ev[1]*ev[2]],seq([coxeter['iprod'](v1-wts[i],u1),ev[i]],i=3..N)];
  k:=1+iquo(degree(ev[1][2],q),2);
  ev:=map((x,y,z)->[x[1],taylor(x[2],y,z-x[1])],ev,q,k);
  ev:=map((x,y,z)->taylor(y^x[1]*x[2],y,z),ev,q,k);
  for j from 0 to k-1 while nops(indets(vars))>0 do
    u1:=convert([seq(coeff(ev[i+1],q,j)*vars[i],i=1..N-2)],`+`)
      -coeff(ev[1],q,j);
    if u1=0 then next else i:=max(op(map(op,indets(u1))),0) fi;
    if i=0 then ERROR(`this cannot happen`) fi;
    u1:=-coeff(u1,vars[i],0)/coeff(u1,vars[i],1);
    vars:=subs(vars[i]=u1,vars);
  od; vars
end:
#
# Given a list of weights wts, and a dominant co-weight u0,
# return a list of evaluations chi(u0,-v), for v in wts,
# normalized so that the lowest order term is 1.
# We are assuming that the orbits are small, so traversal is overkill.
#
`weyl/qtensor/eval`:=proc(wts,u0,cow,coS,q) local wm,res,a,a0,m,orb,i;
  if u0=0 then RETURN([1$nops(wts)]) fi; # the most common case
  res:=[0$nops(wts)]; wm:=weyl['weight_mults'](u0,coS);
  if type(wm,`+`) then wm:=[op(wm)] else wm:=[wm] fi;
  for a in wm do
    a0:=op(indets(a,'indexed')); m:=coeff(a,a0);
    a0:=convert(zip((x,y)->x*y,[op(a0)],cow),`+`);
    orb:=map((x,y)->y-x,coxeter['orbit'](a0,coS),u0);
    res:=[seq(res[i]+m*convert(map((x,y,z)->
      z^coxeter['iprod'](x,y),orb,wts[i],q),`+`),i=1..nops(wts))];
  od; res
end:
#
# Note that we assume that numtheory[cyclotomic](m,q) always returns a
# polynomial with constant term 1. This is not true for m=1, but this
# will never occur in any of our expressions.
#
`weyl/qtensor/cyclosimp`:=proc(ev,dn,q) local vars,sa,i,j;
  vars:=map(op,indets(dn));
  sa:=subs({seq(X[i]=convert([seq(X[j],
    j=numtheory['divisors'](i))],`*`),i=vars)},dn);
  sa:=subs(X[1]=1-q,sa); vars:=map(op,indets(sa));
  sa:=zip((x,y)->[x/denom(y),numer(y)],ev,sa);
  sa:=subs({seq(X[i]=numtheory['cyclotomic'](i,q),i=vars)},sa);
  map(x->normal(x[1])*x[2],sa);
end:
#
# rho(R)    computes half the sum of the positive roots of R.
# co_rho(R) computes half the sum of the positive co-roots.
#
`weyl/rho`:=proc(R) local Rho,S,i,c,eq;
  S:=coxeter['base'](R);
  Rho:=[seq(c[i]*S[i],i=1..nops(S))];
  Rho:=collect(2*convert(Rho,`+`),indets(S));
  eq:={seq(coxeter['iprod'](i,Rho-i)=0,i=S)};
  subs(`coxeter/linsolve`(eq,indets(eq)),Rho/2)
end:
#
`weyl/co_rho`:=proc() coxeter['interior_pt'](args) end:
#
# decompose tensor products of Weyl characters.
# tensor(u,v,R) or tensor(u,v,R,'klimyk') computes the decomposition of
#  the tensor product of irreps of LieAlg(R) with highest weights u and v,
#  using an algorithm based on the Brauer-Klimyk Rule.
# tensor(u,v,R,'mchain') uses the minuscule chain algorithm.
# tensor(u,v,R,'qtensor') uses the quick-and-dirty algorithm based on
#  double specialization of Weyl characters.
#
`weyl/tensor`:=proc(u,v,R) local n;
  if nargs<4 or args[4]='klimyk' then
    `weyl/klimyk`(u,v,{},R)
  elif args[4]='mchain' then
    n:=coxeter['rank'](R);
    `weyl/mchain`(u,v,{$1..n},R)
  elif args[4]='qtensor' then
    `weyl/qtensor`(u,v,R)
  else
    ERROR(cat(`unknown option `,args[4]))
  fi
end:
#
# Given a polynomial expression f in the variables M[..] (orbit sums)
#  and X[..] (Weyl characters), where [..]=fund.wt. coordinates, convert
#  f into a linear combination of orbit sums M[..] using the Brauer-Klimyk
#  Rule (and related techniques).
# The output is "collected" with respect to the variables M[..].
#
`weyl/toM`:=proc(f,R) local S,coS,wt,r,N,g,v,pr,prwc,ct;
  if args[nargs]<>'flat' then
    if assigned(M) or assigned(X) then
      ERROR(`M and X cannot have assigned values`) fi;
    S:=coxeter['base'](R); r:=ilcm(op(map(denom,S)));
    S:=map((x,y)->y*x,S,r); wt:=weyl['weights'](S);
    coS:=coxeter['co_base'](S);
    N:=`weyl/neighbor`(S,{},0);
    g:=`weyl/flatten`(f,S,coS,wt,N);
    `weyl/toM`(g,S,coS,wt,'flat')
  else
    v:=`weyl/flatten/vars`(f,args[4],'X');
    v:=map(proc(x) if x[3]=0 then NULL else x fi end,v);
    g:=subs(X[0$nops(R)]=M[0$nops(R)],f);
    if nops(v)>0 then
      pr:=coxeter['pos_roots'](R);
      prwc:=[seq(map(coxeter['iprod'],args[3],r),r=pr)];
      g:=subs({seq(r[1]=weyl['weight_mults'](r[3],R,pr,prwc),r=v)},g)
    fi;
    v:=indets(g,'indexed');
    v:=map(proc(x) if op(0,x)='M' then x fi end,v);
    ct:=subs({seq(r=0,r=v)},g); g:=g-ct+ct*M[0$nops(R)];
    if ct<>0 then v:={op(v),M[0$nops(R)]} fi;
    collect(g,v,'distributed')
  fi
end:
#
# Given a polynomial expression f in the variables M[..] (orbit sums)
#  and X[..] (Weyl characters), where [..]=fund.wt. coordinates, convert
#  f into a linear combination of Weyl characters X[..] using the
#  Brauer-Klimyk Rule (and related techniques).
# The output is "collected" with respect to the variables X[..].
#
`weyl/toX`:=proc(f,R) local S,coS,wt,r,N,g,v,r0,ct;
  if assigned(M) or assigned(X) then
    ERROR(`M and X cannot have assigned values`) fi;
  S:=coxeter['base'](R); r:=ilcm(op(map(denom,S)));
  S:=map((x,y)->y*x,S,r); wt:=weyl['weights'](S);
  coS:=coxeter['co_base'](S); r0:=convert(wt,`+`);
  N:=`weyl/neighbor`(S,{},0);
  g:=`weyl/flatten`(f,S,coS,wt,N);
  v:=`weyl/flatten/vars`(g,wt,'M');
  g:=subs({seq(r[1]=`weyl/orbit_sum`(r[3],r0,S,coS,N,'skew'),r=v)},g);
  v:=indets(g,'indexed');
  v:=map(proc(x) if op(0,x)='X' then x fi end,v);
  ct:=subs({seq(r=0,r=v)},g); g:=g-ct+ct*X[0$nops(S)];
  if ct<>0 then v:={op(v),X[0$nops(S)]} fi;
  collect(g,v,'distributed');
end:
#
# weight_coords(v,R) returns the coordinates of v with respect to the
#  fundamental weights.
#
`weyl/weight_coords`:=proc(v,R) local S,s;
  S:=coxeter['base'](R);
  [seq(2*coxeter['iprod'](s,v)/coxeter['iprod'](s,s),s=S)]
end:
#
# weight_mults(v,R) produces the dimensions of the weight spaces of the
#  irreducible representation of LieAlg(R) with highest weight v.
#  The output is a sum  c1*M[w1]+c2*M[w2]+...,  where w1,w2,... are the
#  fund.wt.coords of the dominant weights in the weight system of v,
#  and c1,c2,... are their multiplicities.
# The positive roots and their fund.wt.coords can be passed along as
#  the third and fourth arguments.
# weight_mults(v,u,R) computes the multiplicity of u in the irrep of
#  highest weight v. The weight u need not be dominant (but v must be).
#
`weyl/weight_mults`:=proc(v0)
  local S,coS,pr,r0,wts,wl,prwc,m,n,i,j,k,J,SJ,u,v,mults,orb,u0;
  if nargs<>3 then S:=coxeter['base'](args[2]) else
    S:=coxeter['base'](args[3]);
    u:=coxeter['vec2fc'](args[2],S);
    v:=coxeter['root_coords'](v0-u,S);
    if not type(v,'list'('nonnegint')) then RETURN(0) fi
  fi;
  wts:=weyl['weight_sys'](v0,S,'wl');
  if nargs>3 then pr:=args[3]; prwc:=args[4] else
    pr:=coxeter['pos_roots'](S);
    coS:=coxeter['co_base'](S);
    prwc:=[seq(map(coxeter['iprod'],coS,i),i=pr)];
  fi;
  r0:=convert(pr,`+`)/2; mults[1]:=1;
  if nargs=3 then member(u,wts,'n') else n:=nops(wts) fi;
  for i from 2 to n do
    J:=map(proc(x,y) if y[x]=0 then x fi end,[$1..nops(S)],wl[i]);
    v:=wts[i]; SJ:=[seq(S[j],j=J)]; m:=0;
    for k to nops(pr) do
      if min(seq(prwc[k][j],j=J),1)<0 then next fi;
      u:=v+pr[k]; u0:=coxeter['vec2fc'](u,S);
      if not member(u0,wts,'j') then next fi;
      orb:=coxeter['orbit_size'](pr[k],SJ,-1);
      do
        m:=m+mults[j]*orb*coxeter['iprod'](u,pr[k]);
        u:=u+pr[k]; u0:=coxeter['vec2fc'](u,S);
        if not member(u0,wts,'j') then break fi;
      od
    od;
    k:=coxeter['iprod']((v0+r0)$2)-coxeter['iprod']((v+r0)$2);
    mults[i]:=m/k;
  od;
  if nargs=3 then mults[n] else
    convert([seq(mults[i]*M[op(wl[i])],i=1..n)],`+`) fi
end:
#
# weights(R) returns the list of fundamental weights for R.
#
`weyl/weights`:=proc(R) local S,n,A,i,j,B;
  S:=coxeter['base'](R); n:=nops(S);
  if n=0 then RETURN([]) fi;
  A:=[seq([seq(coxeter['iprod'](S[i],S[j]),j=1..i)],i=1..n)];
  A:=array('symmetric',1..n,1..n,A);
  B:=linalg['inverse'](A);
  [seq(A[i,i]*convert([seq(B[i,j]*S[j],j=1..n)],`+`)/2,i=1..n)];
end:
#
# weight_sys(v,R) produces the list of all dominant weights u <= v,
#  where v1 < v2 means v2-v1 is a sum of positive roots.
# weight_sys(v,R,'wl') will do the same, but also assigns to 'wl' the
#  coordinates of the weight vectors w.r.t. the fundamental weights.
# weight_sys(v0,R,'wl',v1) will do the analogous thing for the list of
#  dominant weights v such that v1 <= v <= v0.
#  It does not waste time by checking first to see if v1 <= v0.
# The output is sorted by decreasing height (ht(u):=iprod(u,co_rho(R))).
#  The sorting is an unnecessary expense, but is used in 'qtensor' and
#  'weight_mults' so that equations are solved in a sensible order.
#
`weyl/weight_sys`:=proc(v0,R)
  local coS,pr,prwc,v,w,v1,cr0,peaks,lb,i,j,res,wres,sat;
  coS:=coxeter['base'](R); pr:=`weyl/lsdpr`(coS);
  coS:=coxeter['co_base'](coS); i:=[$1..nops(coS)];
  prwc:=[seq(map(coxeter['iprod'],coS,v),v=pr)];
  peaks:=[seq(map(proc(x,y) if y[x]>0 then x fi end,i,v),v=prwc)];
  if nargs>3 then
    if nargs>4 then w:=args[5] else w:=weyl['weights'](coS) fi;
    v1:=args[4]; cr0:=convert(w,`+`);
    lb:=proc(u,L) local r; for r in L do
      if coxeter['iprod'](u,r)<0 then RETURN(false) fi
    od; true end;
  else
    v1:=0; cr0:=weyl['rho'](coS); lb:=proc() true end;
  fi;
  res:=[v0]; wres:=[map(coxeter['iprod'],coS,v0)];
  for sat while sat<=nops(res) do
    for i to nops(pr) do v:=0;
      for j in peaks[i] while v>=0 do v:=wres[sat][j]-prwc[i][j] od;
      if v>=0 then v:=res[sat]-pr[i];
        if not member(v,res) and lb(v-v1,w) then res:=[op(res),v];
          wres:=[op(wres),zip((x,y)->x-y,wres[sat],prwc[i])]
        fi
      fi
    od
  od;
  res:=[seq([res[i],wres[i],coxeter['iprod'](res[i],cr0)],i=1..nops(res))];
  res:=sort(res,(x,y)->evalb(x[3]>=y[3]));
  if nargs>2 then assign(args[3],map(x->op(2,x),res)) fi;
  map(x->op(1,x),res);
end:
#
# Compute the dimension of the irrep of LieAlg(R) of highest weight v.
# With a third argument q, compute the principal specialization.
# A fourth argument can be used to pass the positive co-roots.
# (For computing dimensions roots or co_roots can be used).
#
`weyl/weyl_dim`:=proc(v,R) local r,pr,r0,q,vr0,f;
  if nargs<3 or args[3]-1=0 then
    if nargs>3 then pr:=args[4]; r0:=weyl['rho'](R)
      else pr:=coxeter['pos_roots'](R); r0:=convert(pr,`+`)/2
    fi;
    f:=[seq(1+coxeter['iprod'](r,v)/coxeter['iprod'](r,r0),r=pr)];
    convert(f,`*`)
  else q:=args[3];
    if nargs>3 then pr:=args[4] else
      pr:=coxeter['pos_roots'](coxeter['co_base'](R)) fi;
    r0:=weyl['rho'](R); vr0:=v+r0;
    f:=[seq((1-q^coxeter['iprod'](r,vr0))/
      (1-q^coxeter['iprod'](r,r0)),r=pr)];
    q^coxeter['iprod'](-v,convert(pr,`+`)/2)*convert(f,`*`)
  fi
end:
#
`coxeter/class_size/E6` := [1,36,240,270,1620,540,1440,1440,5184,4320
,720,540,3240,6480,480,4320,5760,2160,4320,5184,45,540,1440,1440
,80]:
`coxeter/class_rep/E6` := [[],[1],[1,3],[1,2],[1,3,4],[2,3,4,2,5,
4],[1,2,3],[2,3,4,5],[1,2,3,4],[1,2,3,4,2,5,4],[1,2,3,1,
4,2,3,4,5,4,6,5],[1,2,5],[1,2,4,5],[1,2,3,4,5],[1,3,5,6]
,[1,3,4,5,6],[1,2,3,4,2,5,4,6],[1,2,3,5],[1,2,3,4,5,6],
[1,2,3,4,6],[2,3,4,2,3,4,5,4,2,3,4,5],[1,2,3,4,2,3,4,5
,4,2,3,4,5],[1,2,3,5,6],[1,2,3,4,2,3,4,5,4,2,3,4,5,6],
[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5,6]]:
`coxeter/mytype/E6` := [1,-`2`^10,`3`^10,`2`^14,-`4`^7*`2`^3,`2`^6*`4`^6,
-`3`^4*`6`^3*`2`,`3`^4*`6`^3*`2`^3,`5`^7,-`3`^2*`2`*`12`^2*`4`,`6`^4*`3`^4
,-`2`^16,`4`^7*`2`^3,-`4`^3*`8`^3,`3`^11,-`6`^5*`2`*`3`,`9`^4,`3`^2*`2`^
2*`6`^4,`12`^2*`6`^2,-`5`^3*`10`^2,`2`^12,-`2`*`4`^7,-`3`^5*`6`^3*`2`,`6`
^4*`3`^3,`3`^12]:
`coxeter/bpermrep/E6` := {s1 = [[1,35],[9,10],[13,14],[16,18],[17,19]
,[20,22],[21,23],[24,27],[25,28],[29,32]],s2 = [[5,34],[12,15],[
17,20],[19,22],[21,24],[23,27],[25,29],[26,30],[28,32],[31,33]]
,s4 = [[1,32],[3,33],[6,34],[7,15],[13,19],[14,17],[16,21],[18,
23],[29,35],[30,36]],s3 = [[2,35],[5,32],[7,9],[8,14],[11,16],[
19,26],[22,30],[23,31],[27,33],[28,34]],s5 = [[4,36],[5,33],[7,8
],[9,14],[10,13],[21,25],[23,28],[24,29],[27,32],[31,34]],s6 = [
[3,36],[8,11],[13,18],[14,16],[17,21],[19,23],[20,24],[22,27],[
26,31],[30,33]]},36:
`coxeter/irr_chars/E6` := [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1],[1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,-1
,1,1,1,-1,1,-1,-1,1,1],[6,-4,3,2,-2,2,-1,1,1,-1,1,0,0,0
,0,0,0,-1,-1,1,-2,2,2,-2,-3],[6,4,3,2,2,2,1,1,1,1,1,0,0
,0,0,0,0,-1,-1,-1,-2,-2,-2,-2,-3],[10,0,-2,2,0,2,0,0,0,0
,-3,0,-2,0,4,0,1,2,-1,0,-6,0,0,0,1],[15,5,0,3,-1,-1,2,-2
,0,0,1,1,1,-1,3,1,0,0,-1,0,7,3,-1,1,-3],[15,-5,0,3,1,-1
,-2,-2,0,0,1,-1,1,1,3,-1,0,0,-1,0,7,-3,1,1,-3],[15,5,3,-
1,1,3,-1,-1,0,1,2,-3,-1,-1,0,0,0,-1,0,0,-1,1,2,2,6],[15,
-5,3,-1,-1,3,1,-1,0,-1,2,3,-1,1,0,0,0,-1,0,0,-1,-1,-2,2,
6],[20,10,5,4,2,0,1,1,0,-1,-2,2,0,0,-1,-1,-1,1,0,0,4,2,1
,1,2],[20,-10,5,4,-2,0,-1,1,0,1,-2,-2,0,0,-1,1,-1,1,0,0,
4,-2,-1,1,2],[20,0,2,-4,0,4,0,-2,0,0,1,0,0,0,2,0,-1,2,1
,0,4,0,0,-2,-7],[24,4,0,0,0,0,-2,2,-1,0,2,4,0,0,3,1,0,0
,0,-1,8,0,1,-1,6],[24,-4,0,0,0,0,2,2,-1,0,2,-4,0,0,3,-1
,0,0,0,1,8,0,-1,-1,6],[30,10,3,2,0,-2,1,-1,0,-1,-1,-2,0,
0,3,1,0,-1,1,0,-10,-4,1,-1,3],[30,-10,3,2,0,-2,-1,-1,0,1,
-1,2,0,0,3,-1,0,-1,1,0,-10,4,-1,-1,3],[60,10,-3,4,-2,0,1,
-1,0,1,2,2,0,0,-3,-1,0,1,0,0,-4,-2,1,-1,6],[60,-10,-3,4,2
,0,-1,-1,0,-1,2,-2,0,0,-3,1,0,1,0,0,-4,2,-1,-1,6],[60,0,
-6,4,0,4,0,0,0,0,-3,0,0,0,0,0,0,-2,1,0,12,0,0,0,-3],[64
,-16,4,0,0,0,2,0,-1,0,0,0,0,0,-2,0,1,0,0,-1,0,0,2,0,-8]
,[64,16,4,0,0,0,-2,0,-1,0,0,0,0,0,-2,0,1,0,0,1,0,0,-2,0
,-8],[80,0,-4,0,0,0,0,2,0,0,2,0,0,0,2,0,-1,0,0,0,-16,0,
0,2,-10],[81,-9,0,-3,1,-3,0,0,1,0,0,3,-1,-1,0,0,0,0,0,1,
9,-3,0,0,0],[81,9,0,-3,-1,-3,0,0,1,0,0,-3,-1,1,0,0,0,0,0
,-1,9,3,0,0,0],[90,0,0,-6,0,2,0,0,0,0,-3,0,2,0,0,0,0,0
,-1,0,-6,0,0,0,9]]:
#
`coxeter/class_size/E7` := [1,63,672,945,7560,3780,10080,10080,48384,
60480,20160,2240,3780,315,45360,7560,90720,13440,120960,40320,161280
,40320,11340,90720,96768,30240,30240,120960,145152,145152,120960,
60480,207360,207360,60480,3780,315,45360,7560,90720,40320,161280,
120960,40320,13440,96768,11340,90720,10080,10080,48384,60480,2240,
20160,945,7560,3780,672,63,1]:
`coxeter/class_rep/E7` := [[],[1],[1,3],[1,2],[1,3,4],[2,3,4,2,5,
4],[1,2,3],[2,3,4,5],[1,2,3,4],[1,2,3,4,2,5,4],[1,2,3,1,
4,2,3,4,5,4,6,5],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,4,6,5,4
,7,6,5],[1,2,5],[2,5,7],[1,2,4,5],[2,4,5,7],[1,2,3,4,5],
[1,3,5,6],[1,3,4,5,6],[2,4,5,6,7],[1,2,3,4,2,5,4,6],[2,3
,4,2,5,4,6,5,4,7],[2,3,4,2,5,4,7],[2,3,4,2,5,4,6,7],[1
,2,3,1,4,2,3,4,5,4,6,5,7],[1,2,3,5],[2,3,4,5,7],[1,2,3
,4,5,6],[1,2,3,4,6],[2,3,4,5,6,7],[1,2,3,4,2,5,4,6,5,4
,7],[1,3,4,6,7],[1,3,4,5,6,7],[1,2,3,4,2,5,4,6,7],[1,2,
3,4,2,5,4,7],[1,2,5,7],[2,3,4,2,3,4,5,4,2,3,4,5],[1,2,4
,5,7],[1,2,3,4,2,3,4,5,4,2,3,4,5],[1,2,3,4,5,7],[1,2,3
,5,6],[1,2,3,4,5,6,7],[1,2,4,5,6,7],[1,2,3,4,2,3,4,5,4
,2,3,4,5,6],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,5,4,6,
5,7,6],[1,2,3,4,6,7],[2,3,4,2,3,4,5,4,2,3,4,5,6,5,7,6]
,[1,2,3,4,2,3,4,5,4,2,3,4,5,6,5,7,6],[1,2,3,5,7],[2,3,
4,2,3,4,5,4,2,3,4,5,6,7],[1,2,3,4,2,3,4,5,4,2,3,4,5,6
,7],[1,2,3,5,6,7],[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,5,6,5
,4,2,3,4,5,6],[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,5,6,5,4,
2,3,4,5,6,7],[2,3,4,2,3,4,5,4,2,3,4,5,7],[1,2,3,4,2,3,
4,5,4,2,3,4,5,7],[1,2,3,1,4,2,3,4,5,4,2,3,1,4,5,6,5,4
,2,3,4,5,6,7,6,5,4,2,3,4,5,6,7],[1,2,3,4,2,3,4,5,4,2,
3,4,5,6,5,4,2,3,4,5,6,7,6,5,4,2,3,4,5,6,7],[2,3,4,2,3
,4,5,4,2,3,4,5,6,5,4,2,3,4,5,6,7,6,5,4,2,3,4,5,6,7],[
1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,5,4,2,6,5,4,2,3,1,4
,3,5,4,2,6,5,4,3,1,7,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1
,7,6,5,4,2,3,4,5,6,7]]:
`coxeter/mytype/E7` := [1,-`2`^10,`3`^10,`2`^14,-`4`^7*`2`^3,`2`^6*`4`^6,
-`3`^4*`6`^3*`2`,`3`^4*`6`^3*`2`^3,`5`^7,-`3`^2*`2`*`12`^2*`4`,`6`^4*`3`^4
,-`-3`^12,-`2`^16,-`2`^12*`-1`^12,`4`^7*`2`^3,`4`^7*`-1`^6*`2`,-`4`^3*`8`
^3,`3`^11,-`6`^5*`2`*`3`,-`-3`^3*`-1`^3*`6`^4,`9`^4,`6`^3*`-3`^5*`2`*`-1`
,-`4`^6*`-1`^4*`2`^4,`8`^3*`4`^2*`2`*`-1`^2,-`-15`^2*`-5`*`-1`,`3`^2*`2`^2*
`6`^4,-`6`^4*`2`^2*`-1`^2*`-3`^2,`12`^2*`6`^2,-`5`^3*`10`^2,`10`^2*`-5`^3*
`-1`,-`12`^2*`6`^2,-`6`*`12`^2*`4`,`7`^5,-`-7`^5*`-1`,`6`*`12`^2*`4`*`-1`^
2,`2`^16*`-1`^4,`2`^12,-`4`^7*`2`^3*`-1`^2,-`2`*`4`^7,`4`^3*`8`^3,-`3`^5*
`6`^3*`2`,-`-9`^4,`6`^5*`2`*`-3`*`-1`,`6`^4*`3`^3,-`-3`^11*`-1`^3,`15`^2*
`5`,`4`^6*`2`^4,-`8`^3*`4`^2*`2`,-`6`^3*`-3`^4*`2`^3,`6`^3*`-3`^4*`2`*`-1`^
4,-`-5`^7*`-1`,`12`^2*`-3`^2*`2`*`4`,`3`^12,-`6`^4*`-3`^4,-`2`^14*`-1`^8,
`4`^7*`2`^3*`-1`^2,-`2`^6*`4`^6,-`-3`^10*`-1`^6,`2`^10*`-1`^16,-`-1`^36]:
`coxeter/bpermrep/E7` := {s4 = [[5,9],[2,19],[6,24],[15,25],[14,26],[
17,28],[3,29],[22,32],[-31,-34],[33,35]],s3 = [[7,9],[11,14],[3,
21],[4,22],[16,24],[1,25],[8,28],[-19,-30],[-20,-31],[-18,-35]],
s2 = [[5,13],[-18,-19],[3,20],[1,22],[4,25],[12,26],[17,27],[-21
,-31],[-30,-35],[-6,-36]],s1 = [[7,10],[6,19],[12,20],[2,24],[17
,25],[3,26],[4,27],[15,28],[14,29],[18,36]],s6 = [[8,11],[4,20]
,[1,21],[3,25],[17,26],[12,27],[14,28],[15,29],[-22,-31],[32,34
]],s7 = [[11,16],[3,19],[-18,-20],[14,24],[6,26],[2,29],[-21,-30]
,[33,34],[-31,-35],[-12,-36]],s5 = [[8,9],[-6,-12],[10,15],[3,18]
,[-19,-20],[7,28],[-30,-31],[23,32],[-21,-35],[26,36]]},36:
`coxeter/irr_chars/E7` := [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,-1,1,-1,1,
1,-1,1,-1,-1,-1,1,1,-1,1,-1,-1,1,1,-1,1,-1,1,-1,1,-1,1,-1
,-1,1,-1,1,1,1,-1,-1,1,-1,-1,1,1,-1,1,1,-1,-1,1,-1,1,1,-
1,-1,1,-1,-1,1,-1],[7,-5,4,3,-3,3,-2,2,2,-2,2,-2,-1,-1,1,
1,-1,1,-1,-1,1,1,-1,1,-1,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,
-1,1,1,-1,-1,1,-1,-1,1,2,-2,2,-2,-2,2,3,-3,3,4,-5,7],[7,5
,4,3,3,3,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0
,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,
-2,-2,-2,-2,-2,-3,-3,-3,-4,-5,-7],[15,-5,0,3,1,-1,-2,-2,0,0
,1,-3,-1,7,1,-3,1,3,-1,1,0,1,3,-1,0,0,0,-1,0,0,-1,-2,1,
1,-2,-1,7,1,-3,1,1,0,-1,1,3,0,3,-1,-2,-2,0,0,-3,1,3,1,-1
,0,-5,15],[15,5,0,3,-1,-1,2,-2,0,0,1,3,1,-7,1,-3,-1,3,1,
-1,0,1,-3,-1,0,0,0,-1,0,0,1,2,1,-1,-2,-1,7,-1,3,1,-1,0,-
1,1,-3,0,3,1,2,-2,0,0,-3,-1,-3,1,1,0,-5,-15],[21,-11,6,5,
-3,1,-2,2,1,0,-1,3,-3,5,1,-3,-1,0,0,2,0,-2,1,-1,1,2,2,1
,-1,-1,1,0,0,0,0,-3,5,1,-3,-1,-2,0,0,2,0,1,1,-1,2,-2,1,
0,3,-1,5,-3,1,6,-11,21],[21,11,6,5,3,1,2,2,1,0,-1,-3,3,-5
,1,-3,1,0,0,-2,0,-2,-1,-1,-1,2,-2,1,1,-1,-1,0,0,0,0,-3,5
,-1,3,-1,2,0,0,2,0,1,1,1,-2,-2,-1,0,3,1,-5,-3,-1,-6,-11,
-21],[21,9,6,1,3,5,0,0,1,2,3,3,-3,-3,-1,-1,-1,0,0,0,0,0,
1,1,1,-2,-2,-1,-1,-1,-1,0,0,0,0,-3,-3,-1,-1,-1,0,0,0,0,0
,1,1,1,0,0,1,2,3,3,1,3,5,6,9,21],[21,-9,6,1,-3,5,0,0,1
,-2,3,-3,3,3,-1,-1,1,0,0,0,0,0,-1,1,-1,-2,2,-1,1,-1,1,0
,0,0,0,-3,-3,1,1,-1,0,0,0,0,0,1,1,-1,0,0,-1,2,3,-3,-1,3
,-5,-6,9,-21],[27,15,9,7,5,3,3,3,2,1,0,0,3,3,1,1,1,0,0,
0,0,0,-1,-1,-1,1,1,0,0,0,0,-1,-1,-1,-1,3,3,1,1,1,0,0,0,
0,0,-1,-1,-1,3,3,2,1,0,0,7,5,3,9,15,27],[27,-15,9,7,-5,3
,-3,3,2,-1,0,0,-3,-3,1,1,-1,0,0,0,0,0,1,-1,1,1,-1,0,0,0
,0,1,-1,1,-1,3,3,-1,-1,1,0,0,0,0,0,-1,-1,1,-3,3,-2,1,0,
0,-7,5,-3,-9,15,-27],[35,15,5,7,1,-1,3,-1,0,-1,-1,-1,3,11,
1,5,-1,2,0,2,-1,0,3,1,0,1,1,-1,0,0,-1,1,0,0,1,3,11,1,5
,-1,0,-1,0,2,2,0,3,1,-1,3,0,-1,-1,-1,7,1,-1,5,15,35],[35
,-15,5,7,-1,-1,-3,-1,0,1,-1,1,-3,-11,1,5,1,2,0,-2,-1,0,-3
,1,0,1,-1,-1,0,0,1,-1,0,0,1,3,11,-1,-5,-1,0,1,0,2,-2,0,
3,-1,1,3,0,-1,-1,1,-7,1,1,-5,15,-35],[35,-5,5,-5,-1,7,1,-3
,0,-1,3,-1,3,3,-1,-1,1,2,0,0,-1,-2,-1,1,0,1,1,1,0,0,1,-
1,0,0,-1,3,3,-1,-1,1,-2,-1,0,0,2,0,-1,1,-3,1,0,-1,-1,3,-
5,-1,7,5,-5,35],[35,5,5,-5,1,7,-1,-3,0,1,3,1,-3,-3,-1,-1,
-1,2,0,0,-1,-2,1,1,0,1,-1,1,0,0,-1,1,0,0,-1,3,3,1,1,1,2
,1,0,0,-2,0,-1,-1,3,1,0,-1,-1,-3,5,-1,-7,-5,-5,-35],[56,-
24,11,8,-4,0,-3,1,1,1,-2,2,0,-8,0,4,0,2,0,-2,-1,0,0,0,1
,-1,-1,0,1,1,0,-1,0,0,-1,0,-8,0,4,0,0,-1,0,-2,2,1,0,0,1
,-3,1,1,2,-2,8,-4,0,11,-24,56],[56,24,11,8,4,0,3,1,1,-1,-
2,-2,0,8,0,4,0,2,0,2,-1,0,0,0,-1,-1,1,0,-1,1,0,1,0,0,-1
,0,-8,0,-4,0,0,1,0,-2,-2,1,0,0,-1,-3,-1,1,2,2,-8,-4,0,-
11,-24,-56],[70,-10,-5,6,2,2,-1,-1,0,-1,-1,7,-2,-10,-2,2,0,
1,1,-1,1,-1,2,0,0,3,3,-1,0,0,-1,-1,0,0,-1,-2,-10,-2,2,0,
-1,1,1,-1,1,0,2,0,-1,-1,0,-1,7,-1,6,2,2,-5,-10,70],[70,10
,-5,6,-2,2,1,-1,0,1,-1,-7,2,10,-2,2,0,1,-1,1,1,-1,-2,0,0
,3,-3,-1,0,0,1,1,0,0,-1,-2,-10,2,-2,0,1,-1,1,-1,-1,0,2,0
,1,-1,0,-1,7,1,-6,2,-2,5,-10,-70],[84,4,-6,4,0,4,-2,2,-1,
0,-1,3,4,20,0,0,0,3,1,-1,0,1,4,0,-1,-2,-2,1,-1,-1,1,0,0
,0,0,4,20,0,0,0,1,0,1,-1,3,-1,4,0,2,-2,-1,0,3,-1,4,0,4
,-6,4,84],[84,-4,-6,4,0,4,2,2,-1,0,-1,-3,-4,-20,0,0,0,3,-
1,1,0,1,-4,0,1,-2,2,1,1,-1,-1,0,0,0,0,4,20,0,0,0,-1,0,1
,-1,-3,-1,4,0,-2,-2,1,0,3,1,-4,0,-4,6,4,-84],[105,-35,15,5
,-5,5,1,1,0,-1,1,-3,1,1,-1,-1,1,-3,1,1,0,1,1,-1,0,-1,-1
,-1,0,0,-1,1,0,0,1,1,1,-1,-1,1,1,0,1,1,-3,0,1,-1,1,1,0
,-1,-3,1,5,-5,5,15,-35,105],[105,35,15,5,5,5,-1,1,0,1,1,3
,-1,-1,-1,-1,-1,-3,-1,-1,0,1,-1,-1,0,-1,1,-1,0,0,1,-1,0,0
,1,1,1,1,1,1,-1,0,1,1,3,0,1,1,-1,1,0,-1,-3,-1,-5,-5,-5,
-15,-35,-105],[105,5,0,-3,-1,-3,2,2,0,0,2,6,-7,17,-1,3,1,3
,-1,-1,0,-1,1,-1,0,0,0,0,0,0,0,2,0,0,2,-7,17,-1,3,1,-1,
0,-1,-1,3,0,1,-1,2,2,0,0,6,2,-3,-1,-3,0,5,105],[105,-5,0,
-3,1,-3,-2,2,0,0,2,-6,7,-17,-1,3,-1,3,1,1,0,-1,-1,-1,0,0
,0,0,0,0,0,-2,0,0,2,-7,17,1,-3,1,1,0,-1,-1,-3,0,1,1,-2,
2,0,0,6,-2,3,-1,3,0,5,-105],[105,25,0,9,-3,-3,4,-4,0,0,2,
6,1,-7,1,-3,-1,3,1,-1,0,1,-3,-1,0,0,0,0,0,0,0,0,0,0,0,1
,-7,1,-3,-1,1,0,1,-1,3,0,-3,-1,-4,4,0,0,6,2,9,-3,-3,0,25
,105],[105,-25,0,9,3,-3,-4,-4,0,0,2,-6,-1,7,1,-3,1,3,-1,1
,0,1,3,-1,0,0,0,0,0,0,0,0,0,0,0,1,-7,-1,3,-1,-1,0,1,-1
,-3,0,-3,1,4,4,0,0,6,-2,-9,-3,3,0,25,-105],[120,40,15,8,4
,0,1,1,0,-1,-2,-6,0,-8,0,-4,0,0,0,-2,0,-2,0,0,0,-1,-1,0
,0,0,0,1,1,1,1,0,-8,0,-4,0,-2,0,0,-2,0,0,0,0,1,1,0,-1,
-6,-2,8,4,0,15,40,120],[120,-40,15,8,-4,0,-1,1,0,1,-2,6,0,
8,0,-4,0,0,0,2,0,-2,0,0,0,-1,1,0,0,0,0,-1,1,-1,1,0,-8,0
,4,0,2,0,0,-2,0,0,0,0,-1,1,0,-1,-6,2,-8,4,0,-15,40,-120]
,[168,40,6,8,0,0,-2,2,-2,0,2,6,8,8,0,0,0,-3,-1,-1,0,1,0
,0,1,2,2,0,0,0,0,0,0,0,0,8,8,0,0,0,1,0,-1,-1,-3,1,0,0
,2,-2,-2,0,6,2,8,0,0,6,40,168],[168,-40,6,8,0,0,2,2,-2,0
,2,-6,-8,-8,0,0,0,-3,1,1,0,1,0,0,-1,2,-2,0,0,0,0,0,0,0
,0,8,8,0,0,0,-1,0,-1,-1,3,1,0,0,-2,-2,2,0,6,-2,-8,0,0,-
6,40,-168],[189,-39,9,1,-1,-3,3,3,-1,1,0,0,-3,21,-1,-5,-1,0
,0,0,0,0,1,1,-1,1,1,0,1,1,0,-1,0,0,-1,-3,21,-1,-5,-1,0,
0,0,0,0,-1,1,1,3,3,-1,1,0,0,1,-1,-3,9,-39,189],[189,39,9,
1,1,-3,-3,3,-1,-1,0,0,3,-21,-1,-5,1,0,0,0,0,0,-1,1,1,1,-
1,0,-1,1,0,1,0,0,-1,-3,21,1,5,-1,0,0,0,0,0,-1,1,-1,-3,3
,1,1,0,0,-1,-1,3,-9,-39,-189],[189,-51,9,13,1,-3,-3,-3,-1,1
,0,0,-3,-3,1,1,1,0,0,0,0,0,-3,1,-1,1,1,0,-1,-1,0,1,0,0
,1,-3,-3,1,1,1,0,0,0,0,0,-1,-3,1,-3,-3,-1,1,0,0,13,1,-3
,9,-51,189],[189,51,9,13,-1,-3,3,-3,-1,-1,0,0,3,3,1,1,-1,0
,0,0,0,0,3,1,1,1,-1,0,1,-1,0,-1,0,0,1,-3,-3,-1,-1,1,0,0
,0,0,0,-1,-3,-1,3,-3,1,1,0,0,-13,1,3,-9,-51,-189],[189,21,
9,-11,1,9,-3,-3,-1,1,0,0,-3,-3,1,1,-1,0,0,0,0,0,1,-1,-1,
1,1,0,1,1,0,1,0,0,1,-3,-3,1,1,-1,0,0,0,0,0,-1,1,-1,-3,-
3,-1,1,0,0,-11,1,9,9,21,189],[189,-21,9,-11,-1,9,3,-3,-1,-1
,0,0,3,3,1,1,1,0,0,0,0,0,-1,-1,1,1,-1,0,-1,1,0,-1,0,0,
1,-3,-3,-1,-1,-1,0,0,0,0,0,-1,1,1,3,-3,1,1,0,0,11,1,-9,-
9,21,-189],[210,10,-15,10,-2,6,1,1,0,1,-2,-6,2,-14,-2,-2,0,
3,-1,1,0,1,-2,0,0,1,1,0,0,0,0,1,0,0,1,2,-14,-2,-2,0,1,0
,-1,1,3,0,-2,0,1,1,0,1,-6,-2,10,-2,6,-15,10,210],[210,-10,
-15,10,2,6,-1,1,0,-1,-2,6,-2,14,-2,-2,0,3,1,-1,0,1,2,0,0
,1,-1,0,0,0,0,-1,0,0,1,2,-14,2,2,0,-1,0,-1,1,-3,0,-2,0,
-1,1,0,1,-6,2,-10,-2,-6,15,10,-210],[210,50,15,2,2,-2,-1,-1
,0,-1,-1,3,-6,2,-2,2,0,0,0,2,0,2,-2,0,0,-1,-1,1,0,0,1,-
1,0,0,-1,-6,2,-2,2,0,2,0,0,2,0,0,-2,0,-1,-1,0,-1,3,-1,2
,2,-2,15,50,210],[210,-50,15,2,-2,-2,1,-1,0,1,-1,-3,6,-2,-2
,2,0,0,0,-2,0,2,2,0,0,-1,1,1,0,0,-1,1,0,0,-1,-6,2,2,-2
,0,-2,0,0,2,0,0,-2,0,1,-1,0,-1,3,1,-2,2,2,-15,50,-210],[
216,-24,-9,8,4,0,-3,-3,1,-1,0,0,0,24,0,-4,0,0,0,0,0,0,0,
0,1,-1,-1,0,1,1,0,1,-1,-1,1,0,24,0,-4,0,0,0,0,0,0,1,0,0
,-3,-3,1,-1,0,0,8,4,0,-9,-24,216],[216,24,-9,8,-4,0,3,-3,1
,1,0,0,0,-24,0,-4,0,0,0,0,0,0,0,0,-1,-1,1,0,-1,1,0,-1,-
1,1,1,0,24,0,4,0,0,0,0,0,0,1,0,0,3,-3,-1,-1,0,0,-8,4,0
,9,-24,-216],[280,40,-5,8,-4,0,1,-3,0,1,0,-8,0,24,0,4,0,-2
,0,0,1,-2,0,0,0,-1,-1,0,0,0,0,-1,0,0,-1,0,24,0,4,0,-2,1
,0,0,-2,0,0,0,-3,1,0,1,-8,0,8,-4,0,-5,40,280],[280,-40,-5
,8,4,0,-1,-3,0,-1,0,8,0,-24,0,4,0,-2,0,0,1,-2,0,0,0,-1,
1,0,0,0,0,1,0,0,-1,0,24,0,-4,0,2,-1,0,0,2,0,0,0,3,1,0,
1,-8,0,-8,-4,0,5,40,-280],[280,-40,10,-8,0,0,2,-2,0,0,-2,10
,8,-8,0,0,0,1,-1,1,1,-1,0,0,0,-2,-2,0,0,0,0,0,0,0,0,8,
-8,0,0,0,-1,1,-1,1,1,0,0,0,-2,2,0,0,10,-2,-8,0,0,10,-40,
280],[280,40,10,-8,0,0,-2,-2,0,0,-2,-10,-8,8,0,0,0,1,1,-1,
1,-1,0,0,0,-2,2,0,0,0,0,0,0,0,0,8,-8,0,0,0,1,-1,-1,1,-1
,0,0,0,2,2,0,0,10,2,8,0,0,-10,-40,-280],[315,-45,0,3,3,-5
,0,0,0,0,3,-9,3,-21,-1,3,-1,0,0,0,0,0,3,-1,0,0,0,1,0,0
,1,0,0,0,0,3,-21,-1,3,-1,0,0,0,0,0,0,3,-1,0,0,0,0,-9,3
,3,3,-5,0,-45,315],[315,45,0,3,-3,-5,0,0,0,0,3,9,-3,21,-1
,3,1,0,0,0,0,0,-3,-1,0,0,0,1,0,0,-1,0,0,0,0,3,-21,1,-3
,-1,0,0,0,0,0,0,3,1,0,0,0,0,-9,-3,-3,3,5,0,-45,-315],[336
,-16,6,-16,0,0,2,-2,1,0,-2,-6,0,16,0,0,0,0,0,-2,0,2,0,0
,1,2,2,0,-1,-1,0,0,0,0,0,0,16,0,0,0,2,0,0,-2,0,1,0,0,-
2,2,1,0,-6,-2,-16,0,0,6,-16,336],[336,16,6,-16,0,0,-2,-2,1
,0,-2,6,0,-16,0,0,0,0,0,2,0,2,0,0,-1,2,-2,0,1,-1,0,0,0
,0,0,0,16,0,0,0,-2,0,0,-2,0,1,0,0,2,2,-1,0,-6,2,16,0,0
,-6,-16,-336],[378,-30,-9,2,2,6,3,3,-2,-1,0,0,-6,-6,2,2,0,
0,0,0,0,0,-2,0,1,-1,-1,0,0,0,0,-1,0,0,-1,-6,-6,2,2,0,0,
0,0,0,0,1,-2,0,3,3,-2,-1,0,0,2,2,6,-9,-30,378],[378,30,-9
,2,-2,6,-3,3,-2,1,0,0,6,6,2,2,0,0,0,0,0,0,2,0,-1,-1,1,
0,0,0,0,1,0,0,-1,-6,-6,-2,-2,0,0,0,0,0,0,1,-2,0,-3,3,2,
-1,0,0,-2,2,-6,9,-30,-378],[405,45,0,-3,-3,-3,0,0,0,0,0,0,
-3,-27,1,-3,1,0,0,0,0,0,5,1,0,0,0,0,0,0,0,0,-1,-1,0,-3,
-27,1,-3,1,0,0,0,0,0,0,5,1,0,0,0,0,0,0,-3,-3,-3,0,45,405
],[405,-45,0,-3,3,-3,0,0,0,0,0,0,3,27,1,-3,-1,0,0,0,0,0,
-5,1,0,0,0,0,0,0,0,0,-1,1,0,-3,-27,-1,3,1,0,0,0,0,0,0,5
,-1,0,0,0,0,0,0,3,-3,3,0,45,-405],[420,20,0,-12,0,-4,-4,4
,0,0,1,-3,4,4,0,0,0,3,1,1,0,-1,-4,0,0,0,0,-1,0,0,-1,0,
0,0,0,4,4,0,0,0,-1,0,1,1,3,0,-4,0,4,-4,0,0,-3,1,-12,0,-
4,0,20,420],[420,-20,0,-12,0,-4,4,4,0,0,1,3,-4,-4,0,0,0,3
,-1,-1,0,-1,4,0,0,0,0,-1,0,0,1,0,0,0,0,4,4,0,0,0,1,0,1
,1,-3,0,-4,0,-4,-4,0,0,-3,-1,12,0,4,0,20,-420],[512,0,-16,
0,0,0,0,0,2,0,0,8,0,0,0,0,0,-4,0,0,-1,0,0,0,-1,0,0,0,0
,0,0,0,1,1,0,0,0,0,0,0,0,-1,0,0,-4,-1,0,0,0,0,2,0,8,0
,0,0,0,-16,0,512],[512,0,-16,0,0,0,0,0,2,0,0,-8,0,0,0,0,
0,-4,0,0,-1,0,0,0,1,0,0,0,0,0,0,0,1,-1,0,0,0,0,0,0,0,1
,0,0,4,-1,0,0,0,0,-2,0,8,0,0,0,0,16,0,-512]]:
#
`coxeter/class_size/E8` := [1,120,2240,3780,45360,37800,4480,80640,
100800,580608,1209600,806400,268800,37800,907200,1814400,268800,
4838400,6451200,2419200,453600,5443200,11612160,2419200,1161216,604800
,1209600,4838400,5806080,8709120,14515200,3628800,24883200,24883200,
7257600,23224320,14515200,17418240,12902400,2419200,89600,1209600,
113400,3150,2721600,151200,10886400,2419200,1612800,19353600,14515200,
3225600,1612800,29030400,34836480,19353600,1612800,11612160,5443200,
680400,29030400,43545600,10886400,3628800,11612160,9676800,9676800,
907200,403200,15120,1209600,604800,14515200,8709120,5806080,4838400,
7257600,23224320,14515200,24883200,24883200,3628800,89600,2419200,
17418240,12902400,1209600,37800,907200,1814400,2419200,6451200,4838400
,268800,2419200,11612160,453600,5443200,1161216,100800,80640,580608,
1209600,268800,806400,3780,45360,37800,4480,2240,120,1]:
`coxeter/class_rep/E8` := [[],[1],[1,3],[1,2],[1,3,4],[2,3,4,2,5,
4],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,5,6,5,4,2,3,1,4,5
,6,7,6,5,4,2,3,4,5,6,7,8,7,6,5],[1,2,3],[2,3,4,5],[1,2
,3,4],[1,2,3,4,2,5,4],[1,2,3,1,4,2,3,4,5,4,6,5],[1,2,3
,1,4,2,3,1,4,3,5,4,2,3,4,6,5,4,7,6,5],[1,2,5],[1,2,4,5
],[1,2,3,4,5],[1,3,5,6],[1,3,4,5,6],[1,2,3,4,2,5,4,6],[2
,3,4,2,5,4,6,5,4,7],[2,3,4,2,5,4,7],[2,3,4,2,5,4,6,7],
[1,2,3,1,4,2,3,4,5,4,6,5,7],[1,2,3,1,4,2,3,1,4,3,5,4,2
,3,4,6,5,4,7,6,5,8],[1,2,3,1,4,2,3,1,4,5,4,2,3,4,5,6,
5,4,7,6,5,8,7,6],[1,2,3,5],[2,3,4,5,7],[1,2,3,4,5,6],[1
,2,3,4,6],[2,3,4,5,6,7],[1,2,3,4,2,5,4,6,5,4,7],[1,3,4
,6,7],[1,3,4,5,6,7],[1,2,3,4,2,5,4,6,7],[1,2,3,4,2,5,4
,7],[1,2,3,1,4,2,3,4,5,4,6,5,7,6,5,8],[2,3,4,2,5,4,6,5
,4,7,8],[2,3,4,2,5,4,6,7,8],[1,2,3,1,4,2,3,4,5,4,6,5,7
,8],[1,2,3,1,4,2,3,4,5,4,6,5,8],[1,2,3,1,4,2,3,1,4,3,5
,4,2,3,1,4,3,5,4,2,6,5,4,2,3,1,4,5,6,7,6,5,4,2,3,4,5
,6,7,8,7,6],[2,3,4,2,5,4,7,8],[1,2,5,7],[2,3,4,2,3,4,5
,4,2,3,4,5],[1,2,4,5,7],[1,2,3,4,2,3,4,5,4,2,3,4,5],[1
,2,3,4,5,7],[1,2,3,1,4,2,3,4,5,4,2,3,6,5,4,7,6,5,4,8]
,[1,2,3,5,6],[1,2,3,4,5,6,7],[1,2,4,5,6,7],[2,3,4,5,7,8
],[1,2,3,4,2,3,4,5,4,2,3,4,5,6],[2,3,4,5,6,7,8],[1,2,3
,4,2,5,4,6,5,4,7,8],[1,2,3,4,2,5,4,6,8],[1,2,3,1,4,2,3
,1,4,3,5,4,2,3,1,4,3,5,4,6,5,7,6],[1,2,3,4,6,7],[1,3,4
,6,7,8],[2,3,4,2,3,4,5,4,2,3,4,5,6,5,7,6],[1,2,3,4,2,5
,4,6,7,8],[1,3,4,5,6,7,8],[1,2,3,4,2,3,4,5,4,2,3,4,5,6
,5,7,6],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,6,5,4,2,3,4,
7,6,5,4,8,7,6,5],[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,3,5,6
,5,4,7,6,5,8,7,6],[1,2,3,4,2,5,4,7,8],[1,2,3,4,2,3,4,5
,4,2,3,4,5,6,5,4,7,6,5,8,7,6],[2,3,4,2,3,4,5,4,2,3,4,
5,6,5,4,7,6,5,8,7,6],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4
,3,5,4,6,5,4,2,3,4,5,6,7,6,5,4,3,1,8,7,6,5,4,2,3,4,5
,6,7],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,5,4,2,6,5,4,
2,3,4,5,7,6,5,4,2,3,1,4,3,5,4,8,7,6,5,4,2,3,1,4,3,5,4
,2,6,5,4,7,6,5,8,7,6],[1,2,3,5,7],[2,3,4,2,3,4,5,4,2,3
,4,5,6,7],[1,2,3,4,5,6,8],[1,2,3,4,6,8],[1,2,3,4,2,3,4
,5,4,2,3,4,5,6,7],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,
5,4,6,5,7,6,8],[1,2,3,5,6,7],[1,2,3,4,5,6,7,8],[1,2,3,4
,5,7,8],[1,2,4,5,6,7,8],[1,2,3,4,2,3,4,5,4,2,3,4,5,6,5
,7,6,8],[2,3,4,2,3,4,5,4,2,3,4,5,6,5,7,6,8],[1,2,3,1,4
,2,3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5,6],[1,2,3,1,4,2,
3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5,6,7],[1,2,3,4,6,7,8]
,[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5,6,7,6
,8,7],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,5,6,5,4,2,3,1,
4,3,5,6,7,6,5,4,2,3,1,4,3,5,4,6,5,7,6,8,7,6,5],[2,3,4
,2,3,4,5,4,2,3,4,5,7],[1,2,3,4,2,3,4,5,4,2,3,4,5,7],[2
,3,4,2,3,4,5,4,2,3,4,5,6,7,8],[1,2,3,5,6,8],[1,2,3,4,2
,3,4,5,4,2,3,4,5,6,7,8],[1,2,3,4,2,3,4,5,4,2,3,4,5,6,
8],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,4,5,6,5,4,2,3,1,4,5,6
,7,6,5,4,2,3,1,4,3,5,4,6,5,7,6,8,7,6,5],[1,2,3,1,4,2,
3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5,6,7,8],[1,2,3,5,6,7
,8],[1,2,3,1,4,2,3,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5,6,7,
6,5,4,2,3,4,5,6,7],[1,2,3,1,4,2,3,4,5,4,2,3,1,4,5,6,5
,4,2,3,4,5,6,7,6,5,4,2,3,4,5,6,7,8],[1,2,3,1,4,2,3,1,
4,5,4,2,3,1,4,5,6,5,4,2,3,1,4,5,6,7,6,5,4,2,3,1,4,5,6
,7,8,7,6,5,4,2,3,4,5,6,7,8],[2,3,4,2,3,4,5,4,2,3,4,5,
7,8],[1,2,3,4,2,3,4,5,4,2,3,4,5,6,5,4,2,3,4,5,6,7,6,5
,4,2,3,4,5,6,7],[1,2,3,4,2,3,4,5,4,2,3,4,5,6,5,4,2,3,
4,5,6,7,6,5,4,2,3,4,5,6,7,8],[1,2,3,4,2,3,4,5,4,2,3,4
,5,7,8],[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,
5,6,8],[1,2,3,1,4,2,3,1,4,5,4,2,3,1,4,5,6,5,4,2,3,4,5
,6,7,6,5,4,2,3,4,5,6,7,8,7,6,5,4,2,3,4,5,6,7,8],[2,3,
4,2,3,4,5,4,2,3,4,5,6,5,4,2,3,4,5,6,7,6,5,4,2,3,4,5,6
,7],[2,3,4,2,3,4,5,4,2,3,4,5,6,5,4,2,3,4,5,6,7,6,5,4,
2,3,4,5,6,7,8],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,5,4
,2,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7,6,5,4,2,3,1,4,3
,5,4,2,6,5,4,3,1,7,8,7,6,5,4,2,3,4,5,6,7,8],[1,2,3,1,
4,2,3,1,4,3,5,4,2,3,1,4,3,5,6,5,4,2,3,1,4,3,5,4,2,6,7
,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7,8,7,6,5,4,2,3,1,4
,3,5,4,2,6,5,4,3,1,7,6,5,4,2,3,4,5,6,7,8,7,6,5],[1,2,
3,1,4,2,3,1,4,3,5,4,2,3,1,4,3,5,4,2,6,5,4,2,3,1,4,3,5
,4,2,6,5,4,3,1,7,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7,6
,5,4,2,3,4,5,6,7,8],[1,2,3,1,4,2,3,1,4,3,5,4,2,3,1,4,
3,5,4,2,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7,6,5,4,2,3,1
,4,3,5,4,2,6,5,4,3,1,7,6,5,4,2,3,4,5,6,7],[1,2,3,1,4,
2,3,1,4,3,5,4,2,3,1,4,3,5,4,2,6,5,4,2,3,1,4,3,5,4,2,6
,5,4,3,1,7,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7,6,5,4,2
,3,4,5,6,7,8,7,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7,6,5
,4,2,3,4,5,6,7,8,7,6,5,4,2,3,1,4,3,5,4,2,6,5,4,3,1,7
,6,5,4,2,3,4,5,6,7,8]]:
`coxeter/mytype/E8` := [1,`2`^28*`-1`,`3`^28,`2`^44*`-1`^2,`4`^22*`2`^5*
`-2`,`4`^24*`-2`^6,`-3`^40,`6`^6*`3`^16*`2`^10*`-1`,`6`^12*`2`^12*`-3`^4,
`5`^22,`12`^5*`4`^7*`6`^3*`-3`^2*`-2`,`6`^16*`3`^3*`-3`^4,`6`^9*`-3`^21*`2`
,`2`^52*`-1`^3,`4`^22*`2`^11*`-2`*`-1`,`8`^12*`2`^3*`-4`^3,`3`^38,`6`^17*
`3`^3*`2`*`-3`,`9`^13,`6`^14*`2`^2*`-3`^10,`2`^4*`4`^24*`-2`^6*`-1`,`8`^12*
`4`^2*`-4`^3*`-2`,`10`^2*`30`*`6`*`-5`^3*`-15`^3*`-3`,`-6`^2*`12`^8*`-3`^4,
`-5`^24,`2`^14*`6`^10*`3`^8*`-1`^2,`2`^16*`6`^12*`-3`^4*`-1`,`12`^8*`3`^3*
`-6`^2,`2`^3*`5`^12*`10`^5*`-1`,`10`^8*`2`^4*`-5`^6,`12`^8*`6`*`-6`^2*`2`*
`-3`,`3`^6*`4`^7*`2`^2*`12`^5*`6`*`-2`,`7`^17,`14`^4*`-7`^9,`2`^2*`12`^5*
`4`^7*`6`^3*`-3`^2*`-1`*`-2`,`15`^8,`8`^3*`6`*`24`^3*`-4`^3*`-3`^2,`10`*`4`^
2*`20`^4*`-5`^4*`-2`,`-9`^13*`-3`,`6`^17*`2`*`3`*`-3`^4*`-1`,`-3`^39*`-1`^3
,`3`^4*`4`^6*`12`^6*`-2`^6,`2`^56*`-1`^4,`2`^48*`-1`^12,`2`^13*`4`^22*`-2`*
`-1`^2,`4`^22*`2`^9*`-2`*`-1`^6,`2`^5*`8`^12*`-4`^3*`-1`,`-6`^20,`3`^20*`6`
^9*`2`*`-1`,`18`^3*`-9`^7*`2`,`6`^18*`2`^2*`3`*`-3`*`-1`,`3`^4*`2`^3*`6`^15*
`-3`^4,`6`^16*`3`^3*`-3`^3*`-1`^3,`12`^7*`2`*`4`*`-6`^5,`-10`^12,`2`*`9`^7*
`18`^3*`-1`,`6`^9*`-3`^20*`2`*`-1`^3,`3`^3*`5`^7*`15`^5,`2`^2*`4`^28*`-2`^2
,`4`^24*`2`^4*`-2`^6*`-1`^2,`-12`^10,`8`^14*`4`*`-4`,`8`^12*`-4`^3*`4`^2*
`-2`*`-1`,`-4`^30,`-15`^5*`-5`^7*`-3`^3*`-1`,`3`^2*`12`^7*`4`*`6`^3*`-3`^2*
`-2`,`-6`^19*`-2`^3,`2`*`4`^22*`-2`^15,`6`^16*`3`^4*`-3`^4,`-2`^60,`2`^16*
`6`^12*`3`^4*`-1`^3,`2`^14*`6`^10*`-3`^8*`-1`^6,`2`*`12`^8*`3`*`6`*`-6`^2*
`-1`,`2`^4*`5`^6*`10`^8*`-1`^2,`10`^5*`2`^3*`-5`^12*`-1`^3,`-6`^2*`12`^8*
`-3`^3*`-1`^3,`4`^7*`6`^3*`12`^5*`3`^2*`2`^2*`-2`*`-1`,`-15`^8,`3`^2*`8`^3*
`24`^3*`6`*`-4`^3,`7`^9*`14`^4*`-1`,`-7`^17*`-1`,`6`*`4`^7*`2`^2*`12`^5*`-3`
^6*`-2`*`-1`^2,`3`^39,`6`^17*`2`*`-3`*`3`^4,`4`^2*`5`^4*`20`^4*`10`*`-2`,
`9`^13*`3`,`-3`^4*`12`^6*`-2`^6*`4`^6,`2`^52*`-1`^13,`2`^11*`4`^22*`-2`*`-1`
^7,`8`^12*`2`^3*`-4`^3*`-1`^6,`3`^10*`6`^14*`2`^2*`-1`^2,`-9`^13*`-1`^3,`6`
^17*`2`*`3`*`-3`^3*`-1`^4,`-3`^38*`-1`^6,`12`^8*`3`^4*`-6`^2,`5`^3*`6`*`15`^
3*`10`^2*`30`*`3`*`-1`,`4`^24*`-2`^6*`2`^4*`-1`^3,`8`^12*`-4`^3*`4`^2*`-2`*
`-1`^2,`5`^24,`3`^4*`2`^12*`6`^12*`-1`^12,`6`^6*`2`^10*`-3`^16*`-1`^15,`-5`
^22*`-1`^10,`3`^2*`4`^7*`12`^5*`6`^3*`-2`*`-1`^6,`3`^21*`2`*`6`^9*`-1`,`-3`^
3*`6`^16*`3`^4*`-1`^3,`2`^44*`-1`^30,`4`^22*`2`^5*`-2`*`-1`^20,`4`^24*`-2`^6
*`-1`^12,`3`^40,`-3`^28*`-1`^36,`2`^28*`-1`^63,`-1`^120]:
`coxeter/bpermrep/E8` := {s6 = [[3,4],[14,17],[16,19],[18,22],[21,26]
,[25,29],[35,40],[38,44],[41,46],[43,49],[45,51],[48,54],[50,57
],[53,59],[55,62],[60,66],[70,77],[76,84],[81,89],[83,91],[88,
96],[90,97],[94,102],[95,103],[99,106],[101,109],[107,114],[108,
116],[-115]],s2 = [[6,7],[8,9],[10,11],[21,25],[26,29],[30,33],[
31,34],[35,38],[36,39],[40,44],[41,45],[46,51],[52,58],[67,74],
[71,78],[75,82],[79,85],[80,86],[83,90],[87,93],[88,94],[91,97]
,[95,101],[96,102],[98,105],[103,109],[104,110],[111,117],[-119]]
,s1 = [[8,10],[9,11],[12,13],[14,16],[17,19],[20,23],[24,28],[42
,47],[48,53],[54,59],[55,60],[61,65],[62,66],[67,71],[68,72],[
69,73],[74,78],[75,80],[76,81],[82,86],[83,88],[84,89],[90,94],
[91,96],[97,102],[98,104],[105,110],[112,118],[-120]],s8 = [[1,2],
[20,24],[23,28],[27,32],[31,36],[34,39],[35,41],[38,45],[40,46]
,[43,50],[44,51],[48,55],[49,57],[53,60],[54,62],[56,63],[59,66
],[61,68],[65,72],[67,75],[71,80],[74,82],[78,86],[79,87],[85,
93],[92,100],[99,107],[106,114],[-113]],s5 = [[4,5],[12,14],[13,16
],[15,18],[26,30],[29,33],[31,35],[34,38],[36,41],[39,45],[49,
56],[54,61],[57,63],[59,65],[62,68],[64,70],[66,72],[69,76],[73
,81],[91,98],[92,99],[96,104],[97,105],[100,107],[102,110],[103,
111],[108,115],[109,117],[-116]],s4 = [[5,6],[9,12],[11,13],[18,21
],[22,26],[27,31],[32,36],[33,37],[38,43],[44,49],[45,50],[51,
57],[58,64],[61,67],[65,71],[68,75],[72,80],[76,83],[81,88],[84
,91],[85,92],[89,96],[93,100],[101,108],[105,112],[109,116],[110
,118],[111,119],[-117]],s7 = [[2,3],[17,20],[19,23],[22,27],[26,
31],[29,34],[30,35],[33,38],[37,43],[42,48],[46,52],[47,53],[51
,58],[57,64],[62,69],[63,70],[66,73],[68,76],[72,81],[75,83],[
80,88],[82,90],[86,94],[87,95],[93,101],[100,108],[106,113],[107
,115],[-114]],s3 = [[6,8],[7,9],[13,15],[16,18],[19,22],[23,27],
[28,32],[37,42],[43,48],[49,54],[50,55],[56,61],[57,62],[63,68]
,[64,69],[70,76],[71,79],[77,84],[78,85],[80,87],[86,93],[88,95
],[94,101],[96,103],[102,109],[104,111],[110,117],[112,120],[-118]
]},120:
`coxeter/irr_chars/E8` := [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],[1,-1,1,1,-1,1,1,-1,
1,1,-1,1,-1,-1,1,-1,1,-1,1,1,-1,1,-1,1,1,1,-1,1,-1,1,-1,
-1,1,-1,1,1,-1,-1,1,-1,1,1,1,1,-1,-1,1,1,-1,-1,1,1,1,-1,
1,-1,-1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,1
,-1,-1,1,-1,1,-1,-1,1,1,-1,1,-1,1,1,-1,1,1,-1,-1,1,1,1,-
1,1,-1,-1,1,1,-1,1,1,1,-1,1],[8,-6,5,4,-4,4,4,-3,3,3,-3,
3,-3,-2,2,-2,2,-2,2,2,-2,2,-2,2,2,1,-1,1,-1,1,-1,-1,1,-1
,1,1,-1,-1,1,-1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,0,1,-1,1,-1,1,-1,-1,-1,1,1,-1,1,-1
,1,1,-1,-1,2,-2,2,-2,-2,2,-2,-2,2,2,-2,-2,-3,3,-3,3,3,-3
,-4,4,-4,-4,-5,6,-8],[8,6,5,4,4,4,4,3,3,3,3,3,3,2,2,2,2
,2,2,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2
,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-3,-3,-3,-3,-3,-3,-4,-4,-4,-
4,-5,-6,-8],[28,14,10,4,6,8,10,2,2,3,4,5,5,-2,0,0,1,1,1,
1,2,2,2,2,3,-2,-2,-1,-1,-1,-1,0,0,0,0,0,0,1,1,1,1,2,-4,
-4,-2,-2,-2,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,0,0,0,0
,1,1,2,2,4,-2,-2,-1,-1,-1,-1,0,0,0,0,0,0,1,1,1,1,2,-2,0
,0,1,1,1,1,2,2,2,2,3,2,2,3,4,5,5,4,6,8,10,10,14,28],[28
,-14,10,4,-6,8,10,-2,2,3,-4,5,-5,2,0,0,1,-1,1,1,-2,2,-2,
2,3,-2,2,-1,1,-1,1,0,0,0,0,0,0,-1,1,-1,1,2,-4,-4,2,2,-2
,-2,1,1,-1,-1,-1,1,-1,1,1,0,0,0,0,0,0,0,0,-1,1,-2,2,4,2
,-2,1,-1,1,-1,0,0,0,0,0,0,1,-1,-1,1,2,2,0,0,1,1,-1,1,2
,-2,-2,2,3,2,-2,3,-4,-5,5,4,-6,8,10,10,-14,28],[35,-21,14,
11,-9,7,5,-6,6,5,-4,3,-3,-5,3,-3,2,-2,2,2,-1,1,-1,1,0,2,
-2,1,-1,1,-1,0,0,0,0,0,0,1,-1,1,-1,-2,3,3,-1,-1,1,1,0,0
,0,0,0,0,0,0,0,-1,-1,-1,-1,1,1,-1,-1,2,-2,3,-3,-5,-2,2,-
1,1,-1,1,0,0,0,0,0,0,-1,1,1,-1,-2,-5,3,-3,2,2,-2,2,1,-1
,-1,1,0,6,-6,5,-4,-3,3,11,-9,7,5,14,-21,35],[35,21,14,11,9
,7,5,6,6,5,4,3,3,5,3,3,2,2,2,2,1,1,1,1,0,2,2,1,1,1,1
,0,0,0,0,0,0,-1,-1,-1,-1,-2,3,3,1,1,1,1,0,0,0,0,0,0,0,
0,0,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2,-3,-3,-5,2,2,1,1,1,1,0,
0,0,0,0,0,-1,-1,-1,-1,-2,5,3,3,2,2,2,2,1,1,1,1,0,6,6,5
,4,3,3,11,9,7,5,14,21,35],[50,-20,5,10,0,-2,5,-5,-3,0,1,0
,-2,-4,2,2,5,-1,-1,1,0,0,0,1,0,1,-1,-2,0,0,0,-3,1,1,-1,
0,-1,0,-1,2,-4,1,2,18,0,-8,0,1,1,1,-1,3,3,-1,0,1,1,0,2,
6,-1,0,-2,2,0,1,1,-4,-3,10,-1,1,0,0,0,-2,-1,0,-1,1,1,-3,
-4,2,0,-1,1,-4,2,2,1,-1,-1,5,1,0,0,0,0,-3,-5,0,1,-2,0,10
,0,-2,5,5,-20,50],[50,20,5,10,0,-2,5,5,-3,0,-1,0,2,4,2,-2
,5,1,-1,1,0,0,0,1,0,1,1,-2,0,0,0,3,1,-1,-1,0,1,0,-1,-2
,-4,1,2,18,0,8,0,1,-1,-1,-1,3,3,1,0,-1,-1,0,2,6,-1,0,2,
2,0,-1,1,4,-3,10,1,1,0,0,0,-2,-1,0,1,-1,1,3,-4,-2,0,-1,1
,4,2,-2,1,-1,1,5,1,0,0,0,0,-3,5,0,-1,2,0,10,0,-2,5,5,20
,50],[56,-14,11,-4,-4,12,16,1,-3,1,-3,6,-4,6,-2,2,2,0,-1,-
2,-2,2,-1,0,4,-1,3,0,1,-1,2,-1,0,0,-1,-1,1,-1,1,0,-2,3,0
,0,0,0,0,0,-2,-1,0,0,0,0,0,1,2,1,0,0,0,0,0,0,-1,0,0,0
,0,0,-3,1,-2,1,-1,0,1,1,-1,0,0,1,2,0,1,-1,-3,-6,2,-2,2,
1,0,-2,0,1,2,-2,-4,3,-1,-1,3,4,-6,4,4,-12,-16,-11,14,-56],
[56,14,11,-4,4,12,16,-1,-3,1,3,6,4,-6,-2,-2,2,0,-1,-2,2,2
,1,0,4,-1,-3,0,-1,-1,-2,1,0,0,-1,-1,-1,1,1,0,-2,3,0,0,0
,0,0,0,2,1,0,0,0,0,0,-1,-2,1,0,0,0,0,0,0,-1,0,0,0,0,0
,3,1,2,1,1,0,1,1,1,0,0,-1,2,0,-1,-1,-3,6,2,2,2,1,0,-2,
0,-1,-2,-2,-4,3,1,-1,-3,-4,-6,4,-4,-12,-16,-11,-14,-56],[70,
0,10,-10,0,14,19,0,-6,0,0,6,0,0,-2,0,4,0,-2,-4,0,2,0,-1,
5,2,0,2,0,0,0,0,0,0,-2,-1,0,0,1,0,-2,2,6,6,0,0,2,3,0,0
,0,0,0,0,1,0,0,0,-2,-2,-1,0,0,2,0,0,0,0,3,6,0,2,0,0,0
,2,-2,-1,0,0,0,0,-2,0,0,1,2,0,-2,0,-4,-2,0,4,-1,0,0,2,5
,-6,0,0,0,0,6,-10,0,14,19,10,0,70],[84,42,21,20,10,4,-6,9
,5,4,1,-1,-3,10,4,2,3,1,0,-1,2,0,-1,-2,-1,5,1,1,2,0,-1,
1,0,0,1,-1,-1,0,0,1,3,1,4,20,2,10,0,-2,3,0,1,-1,5,-1,-1
,0,3,1,0,4,0,0,2,0,1,1,1,2,2,4,1,5,-1,0,2,1,1,-1,-1,0
,0,1,3,1,0,0,1,10,4,2,-1,0,1,3,-2,-1,2,0,-1,5,9,4,1,-3
,-1,20,10,4,-6,21,42,84],[84,-42,21,20,-10,4,-6,-9,5,4,-1,-
1,3,-10,4,-2,3,-1,0,-1,-2,0,1,-2,-1,5,-1,1,-2,0,1,-1,0,0
,1,-1,1,0,0,-1,3,1,4,20,-2,-10,0,-2,-3,0,1,-1,5,1,-1,0,-
3,1,0,4,0,0,-2,0,1,-1,1,-2,2,4,-1,5,1,0,-2,1,1,-1,1,0,0
,-1,3,-1,0,0,1,-10,4,-2,-1,0,-1,3,-2,1,-2,0,-1,5,-9,4,-1
,3,-1,20,-10,4,-6,21,-42,84],[112,-56,31,24,-16,8,-4,-11,9,7
,-3,0,2,-8,4,-4,4,-2,1,0,0,0,1,-2,-2,3,-1,2,-1,1,0,-1,0
,0,-1,-1,1,1,-1,2,-4,-1,0,0,0,0,0,0,-2,-1,0,0,0,0,0,1,2
,1,0,0,0,0,0,0,-1,0,0,0,0,0,1,-3,0,-1,1,-2,1,1,-1,0,0,
1,4,-2,-1,1,1,8,-4,4,0,-1,2,-4,2,-1,0,0,2,-9,11,-7,3,-2,
0,-24,16,-8,4,-31,56,-112],[112,56,31,24,16,8,-4,11,9,7,3,0
,-2,8,4,4,4,2,1,0,0,0,-1,-2,-2,3,1,2,1,1,0,1,0,0,-1,-1
,-1,-1,-1,-2,-4,-1,0,0,0,0,0,0,2,1,0,0,0,0,0,-1,-2,1,0,
0,0,0,0,0,-1,0,0,0,0,0,-1,-3,0,-1,-1,-2,1,1,1,0,0,-1,4,
2,1,1,1,-8,-4,-4,0,-1,-2,-4,2,1,0,0,2,-9,-11,-7,-3,2,0,-
24,-16,-8,4,-31,-56,-112],[160,64,34,16,16,16,20,4,6,5,6,6,8
,0,0,0,-2,0,1,2,0,0,1,2,0,-2,0,-2,-1,-1,0,-2,-1,-1,0,0,
0,-1,-1,0,2,-2,0,0,0,0,0,0,-2,-1,0,0,0,0,0,1,2,-1,0,0,0
,0,0,0,1,0,0,0,0,0,0,2,0,1,1,2,0,0,0,1,1,2,-2,0,1,1,2
,0,0,0,-2,-1,0,2,-2,-1,0,0,0,-6,-4,-5,-6,-8,-6,-16,-16,-16
,-20,-34,-64,-160],[160,-64,34,16,-16,16,20,-4,6,5,-6,6,-8,0
,0,0,-2,0,1,2,0,0,-1,2,0,-2,0,-2,1,-1,0,2,-1,1,0,0,0,1
,-1,0,2,-2,0,0,0,0,0,0,2,1,0,0,0,0,0,-1,-2,-1,0,0,0,0,
0,0,1,0,0,0,0,0,0,2,0,1,-1,2,0,0,0,-1,1,-2,-2,0,-1,1,2
,0,0,0,-2,-1,0,2,-2,1,0,0,0,-6,4,-5,6,8,-6,-16,16,-16,-20
,-34,64,-160],[168,0,-12,8,0,8,15,0,4,-2,0,-2,0,0,0,0,6,0
,0,2,0,0,0,-1,3,-4,0,2,0,-2,0,0,0,0,0,0,0,0,0,0,6,-4,8
,40,0,0,0,3,0,0,2,4,-2,0,-1,0,0,-2,0,8,1,0,0,4,-2,0,0,
0,7,24,0,-4,0,-2,0,2,0,0,0,0,0,0,6,0,0,0,-4,0,0,0,2,0,
0,6,-1,0,0,0,3,4,0,-2,0,0,-2,8,0,8,15,-12,0,168],[175,-35
,-5,15,5,-1,-5,-5,-5,0,-1,1,1,-3,-1,1,4,0,1,0,5,-1,0,-1,
0,3,3,-1,0,0,-1,-1,0,0,-1,0,1,0,1,-3,13,-1,-1,-17,-3,5,-1
,3,-2,1,2,-2,-2,0,0,1,-2,0,3,-1,-1,-1,1,-1,0,2,0,-3,-5,
15,3,3,-1,0,0,-1,-1,0,1,0,0,-1,13,-3,0,1,-1,-3,-1,1,0,1,
0,4,-1,0,5,-1,0,-5,-5,0,-1,1,1,15,5,-1,-5,-5,-35,175],[175
,35,-5,15,-5,-1,-5,5,-5,0,1,1,-1,3,-1,-1,4,0,1,0,-5,-1,0
,-1,0,3,-3,-1,0,0,1,1,0,0,-1,0,-1,0,1,3,13,-1,-1,-17,3,-
5,-1,3,2,-1,2,-2,-2,0,0,-1,2,0,3,-1,-1,1,-1,-1,0,-2,0,3,
-5,15,-3,3,1,0,0,-1,-1,0,-1,0,0,1,13,3,0,1,-1,3,-1,-1,0,
1,0,4,-1,0,-5,-1,0,-5,5,0,1,-1,1,15,-5,-1,-5,-5,35,175],[
210,-84,39,26,-16,6,-15,-9,7,5,-1,-2,6,-4,2,-2,3,-1,0,-1,0
,0,1,-3,0,-1,-1,0,1,1,0,-1,0,0,-1,0,1,-1,0,2,-6,3,2,-14
,0,8,0,1,3,0,-1,1,-5,-1,0,0,3,-1,2,-2,1,0,2,-2,-1,-1,1,
-4,1,10,-1,-1,0,1,1,0,-1,0,1,0,0,-1,-6,2,-1,0,3,-4,2,-2,
-1,0,-1,3,-3,1,0,0,0,7,-9,5,-1,6,-2,26,-16,6,-15,39,-84,
210],[210,84,39,26,16,6,-15,9,7,5,1,-2,-6,4,2,2,3,1,0,-1,0
,0,-1,-3,0,-1,1,0,-1,1,0,1,0,0,-1,0,-1,1,0,-2,-6,3,2,-14
,0,-8,0,1,-3,0,-1,1,-5,1,0,0,-3,-1,2,-2,1,0,-2,-2,-1,1,1
,4,1,10,1,-1,0,1,-1,0,-1,0,-1,0,0,1,-6,-2,1,0,3,4,2,2,-
1,0,1,3,-3,-1,0,0,0,7,9,5,1,-6,-2,26,16,6,-15,39,84,210],
[300,-90,30,20,-10,8,30,0,6,0,-2,3,-9,-10,0,0,-6,2,0,2,2,-
2,0,2,0,2,-4,-1,0,0,-1,2,-1,1,0,0,0,0,0,-1,3,2,12,12,-2
,-2,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,2,-6,6,20
,-4,2,-1,0,0,-1,0,0,0,1,-1,2,3,-1,0,0,2,-10,0,0,2,0,2,-
6,2,0,2,-2,0,6,0,0,-2,-9,3,20,-10,8,30,30,-90,300],[300,90
,30,20,10,8,30,0,6,0,2,3,9,10,0,0,-6,-2,0,2,-2,-2,0,2,0
,2,4,-1,0,0,1,-2,-1,-1,0,0,0,0,0,1,3,2,12,12,2,2,2,2,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2,6,6,20,4,2,1,0,0
,-1,0,0,0,-1,-1,-2,3,1,0,0,2,10,0,0,2,0,-2,-6,2,0,-2,-2
,0,6,0,0,2,9,3,20,10,8,30,30,90,300],[350,-70,35,-10,-10,26
,35,5,-5,0,-5,7,-7,10,-2,4,-1,1,-1,-1,-2,0,0,-1,0,-1,1,-1
,0,0,1,-1,0,0,1,0,1,0,-1,1,-1,-1,-2,-2,-2,-2,0,-1,-1,-1,
1,1,1,-1,0,-1,-1,0,2,2,1,0,0,-2,0,1,-1,2,-5,-10,1,-1,1,0
,0,-1,1,0,1,0,0,-1,-1,1,0,-1,-1,10,-2,4,-1,-1,1,-1,-1,0,
-2,0,0,-5,5,0,-5,-7,7,-10,-10,26,35,35,-70,350],[350,70,35,-
10,10,26,35,-5,-5,0,5,7,7,-10,-2,-4,-1,-1,-1,-1,2,0,0,-1,0
,-1,-1,-1,0,0,-1,1,0,0,1,0,-1,0,-1,-1,-1,-1,-2,-2,2,2,0,
-1,1,1,1,1,1,1,0,1,1,0,2,2,1,0,0,-2,0,-1,-1,-2,-5,-10,-1
,-1,-1,0,0,-1,1,0,-1,0,0,1,-1,-1,0,-1,-1,-10,-2,-4,-1,-1,
-1,-1,-1,0,2,0,0,-5,-5,0,5,7,7,-10,10,26,35,35,70,350],[400
,120,25,40,0,-8,20,15,-9,0,-3,0,6,8,4,-4,10,2,-2,2,0,0,0
,2,0,1,1,-2,0,0,0,3,1,-1,-1,0,1,0,-1,-2,-4,1,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,-1,
0,0,0,2,1,0,-1,1,-1,-3,4,2,0,1,-1,-8,-4,4,-2,2,-2,-10,-2
,0,0,0,0,9,-15,0,3,-6,0,-40,0,8,-20,-25,-120,-400],[400,-120
,25,40,0,-8,20,-15,-9,0,3,0,-6,-8,4,4,10,-2,-2,2,0,0,0,2
,0,1,-1,-2,0,0,0,-3,1,1,-1,0,-1,0,-1,2,-4,1,0,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0,0
,0,2,1,0,1,-1,-1,3,4,-2,0,1,-1,8,-4,-4,-2,2,2,-10,-2,0,0
,0,0,9,15,0,-3,6,0,-40,0,8,-20,-25,120,-400],[420,0,-30,20,
0,12,24,0,2,0,0,-4,0,0,-4,0,6,0,0,2,0,0,0,0,5,2,0,0,0,
0,0,0,0,0,2,-1,0,0,0,0,-12,-6,4,-28,0,0,0,0,0,0,-2,-4,2
,0,1,0,0,0,4,-4,0,0,0,0,0,0,0,0,8,36,0,2,0,0,0,0,2,-1
,0,0,0,0,-12,0,0,0,-6,0,-4,0,2,0,0,6,0,0,0,0,5,2,0,0,0
,0,-4,20,0,12,24,-30,0,420],[448,0,28,-32,0,32,44,0,-12,-2,
0,6,0,0,0,0,4,0,-2,-4,0,0,0,-2,2,4,0,2,0,2,0,0,0,0,0,1
,0,0,-1,0,2,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-2,0,0,0
,0,0,0,2,0,0,0,0,0,0,-4,0,-2,0,-2,0,-1,0,0,0,0,-2,0,0,
1,4,0,0,0,4,2,0,-4,2,0,0,0,-2,12,0,2,0,0,-6,32,0,-32,-44
,-28,0,-448],[448,-112,16,32,0,0,-16,-4,0,-2,0,0,4,-16,0,0,
-2,2,1,-2,0,0,-1,0,2,8,4,0,-2,-2,0,0,0,0,0,1,0,0,-1,4,-
16,0,0,0,0,0,0,0,-4,1,0,0,0,0,0,-1,4,1,0,0,0,0,0,0,-1,
0,0,0,0,0,-4,-8,0,2,2,0,0,-1,0,0,0,0,16,-4,0,1,0,16,0,0
,2,-1,-2,2,0,1,0,0,-2,0,4,2,0,-4,0,-32,0,0,16,-16,112,-
448],[448,112,16,32,0,0,-16,4,0,-2,0,0,-4,16,0,0,-2,-2,1,-2
,0,0,1,0,2,8,-4,0,2,-2,0,0,0,0,0,1,0,0,-1,-4,-16,0,0,0
,0,0,0,0,4,-1,0,0,0,0,0,1,-4,1,0,0,0,0,0,0,-1,0,0,0,0
,0,4,-8,0,2,-2,0,0,-1,0,0,0,0,16,4,0,1,0,-16,0,0,2,-1,2
,2,0,-1,0,0,-2,0,-4,2,0,4,0,-32,0,0,16,-16,-112,-448],[525
,-105,30,5,-5,-7,30,0,6,0,2,0,-6,7,-3,-3,3,1,0,-1,3,-1,0
,2,0,2,4,2,0,0,0,-2,0,0,0,0,0,0,0,-2,12,2,-19,45,3,-13,
-1,2,-3,0,-1,3,3,-1,0,0,-3,0,1,1,0,1,1,-3,0,-1,-1,-1,6,5
,4,2,0,0,0,2,0,0,0,0,0,-2,12,-2,0,0,2,7,-3,-3,-1,0,1,3
,2,0,3,-1,0,6,0,0,2,-6,0,5,-5,-7,30,30,-105,525],[525,105,
30,5,5,-7,30,0,6,0,-2,0,6,-7,-3,3,3,-1,0,-1,-3,-1,0,2,0,
2,-4,2,0,0,0,2,0,0,0,0,0,0,0,2,12,2,-19,45,-3,13,-1,2,3
,0,-1,3,3,1,0,0,3,0,1,1,0,-1,-1,-3,0,1,-1,1,6,5,-4,2,0
,0,0,2,0,0,0,0,0,2,12,2,0,0,2,-7,-3,3,-1,0,-1,3,2,0,-3
,-1,0,6,0,0,-2,6,0,5,5,-7,30,30,105,525],[560,196,74,56,24
,8,-20,16,6,5,0,-3,-7,12,4,0,2,0,-1,-2,4,0,-1,-2,0,2,0,-
1,1,-1,-1,0,0,0,2,0,0,1,1,3,7,2,0,0,0,0,0,0,-2,-1,0,0,
0,0,0,1,2,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,-2,1,1,-1,1,-2,
0,0,0,0,0,-7,-3,-1,-1,-2,-12,-4,0,2,1,0,-2,2,1,-4,0,0,-6
,-16,-5,0,7,3,-56,-24,-8,20,-74,-196,-560],[560,-196,74,56,-24
,8,-20,-16,6,5,0,-3,7,-12,4,0,2,0,-1,-2,-4,0,1,-2,0,2,0,
-1,-1,-1,1,0,0,0,2,0,0,-1,1,-3,7,2,0,0,0,0,0,0,2,1,0,0
,0,0,0,-1,-2,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,-2,-1,1,1,1,
-2,0,0,0,0,0,-7,3,1,-1,-2,12,-4,0,2,1,0,-2,2,-1,4,0,0,-6
,16,-5,0,-7,3,-56,24,-8,20,-74,196,-560],[567,189,81,39,29,15
,0,9,9,7,3,0,0,-3,-1,1,0,0,0,0,-3,-1,-1,0,-3,-3,-3,0,-1
,-1,0,-1,0,0,-1,0,1,-1,0,0,0,-3,-9,-9,-3,-3,-1,0,0,0,0,0
,0,0,1,0,0,1,-1,-1,0,1,1,3,1,0,0,-3,0,-9,-3,-3,0,-1,-1,
0,-1,0,1,0,0,-1,0,0,-1,0,-3,-3,-1,1,0,0,0,0,0,-1,-3,-1,-
3,9,9,7,3,0,0,39,29,15,0,81,189,567],[567,-189,81,39,-29,15
,0,-9,9,7,-3,0,0,3,-1,-1,0,0,0,0,3,-1,1,0,-3,-3,3,0,1,-
1,0,1,0,0,-1,0,-1,1,0,0,0,-3,-9,-9,3,3,-1,0,0,0,0,0,0,0
,1,0,0,1,-1,-1,0,-1,-1,3,1,0,0,3,0,-9,3,-3,0,-1,1,0,-1,
0,-1,0,0,1,0,0,1,0,-3,3,-1,-1,0,0,0,0,0,1,3,-1,-3,9,-9,
7,-3,0,0,39,-29,15,0,81,-189,567],[700,-210,55,60,-10,-4,10,-
15,-1,0,3,-1,-3,-18,4,2,4,0,-2,0,-2,0,0,2,0,3,-3,-1,0,0,
1,-1,0,0,1,0,-1,0,1,-3,7,-1,12,-4,-2,6,0,-2,0,0,0,2,-4,0
,0,0,0,0,0,-4,0,0,2,0,0,0,-2,6,2,-20,-3,3,1,0,0,-1,1,0
,-1,0,0,-1,7,-3,0,1,-1,-18,4,2,0,-2,0,4,2,0,-2,0,0,-1,-
15,0,3,-3,-1,60,-10,-4,10,55,-210,700],[700,210,55,60,10,-4,
10,15,-1,0,-3,-1,3,18,4,-2,4,0,-2,0,2,0,0,2,0,3,3,-1,0,0
,-1,1,0,0,1,0,1,0,1,3,7,-1,12,-4,2,-6,0,-2,0,0,0,2,-4,0
,0,0,0,0,0,-4,0,0,-2,0,0,0,-2,-6,2,-20,3,3,-1,0,0,-1,1,
0,1,0,0,1,7,3,0,1,-1,18,4,-2,0,-2,0,4,2,0,2,0,0,-1,15,0
,-3,3,-1,60,10,-4,10,55,210,700],[700,-70,-20,20,10,0,-20,-10
,-4,0,-2,2,2,10,0,0,7,1,1,-1,6,-2,0,0,0,-4,-2,0,0,0,0,-
2,0,0,0,0,0,0,1,-2,-2,0,-4,92,2,-14,2,-4,5,-1,-1,-1,-1,1
,0,-1,5,0,0,8,0,0,0,0,0,1,-1,-2,-4,20,-2,-4,0,0,0,0,0,0
,0,0,0,-2,-2,-2,0,1,0,10,0,0,-1,1,1,7,0,0,6,-2,0,-4,-10
,0,-2,2,2,20,10,0,-20,-20,-70,700],[700,70,-20,20,-10,0,-20,
10,-4,0,2,2,-2,-10,0,0,7,-1,1,-1,-6,-2,0,0,0,-4,2,0,0,0,
0,2,0,0,0,0,0,0,1,2,-2,0,-4,92,-2,14,2,-4,-5,1,-1,-1,-1,
-1,0,1,-5,0,0,8,0,0,0,0,0,-1,-1,2,-4,20,2,-4,0,0,0,0,0,
0,0,0,0,2,-2,2,0,1,0,-10,0,0,-1,1,-1,7,0,0,-6,-2,0,-4,10
,0,2,-2,2,20,-10,0,-20,-20,70,700],[840,126,21,4,4,20,60,-9
,3,-5,3,3,9,10,2,-2,-6,-2,0,2,2,-2,-1,-2,0,1,5,-1,1,1,1
,1,0,0,1,0,-1,1,0,-1,-3,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,-5,-1,-1,-1,-1,1,-1,0,1,0
,0,-1,3,1,-1,0,-5,-10,-2,2,-2,0,2,6,2,1,-2,2,0,-3,9,5,-3
,-9,-3,-4,-4,-20,-60,-21,-126,-840],[840,-126,21,4,-4,20,60,9
,3,-5,-3,3,-9,-10,2,2,-6,2,0,2,-2,-2,1,-2,0,1,-5,-1,-1,1
,-1,-1,0,0,1,0,1,-1,0,1,-3,5,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,5,-1,1,-1,1,1,-1,0,-1,
0,0,1,3,-1,1,0,-5,10,-2,-2,-2,0,-2,6,2,-1,2,2,0,-3,-9,5,
3,9,-3,-4,4,-20,-60,-21,126,-840],[840,-84,-24,24,4,16,30,6,8
,-5,-2,-1,-3,-20,0,0,3,1,0,3,-4,0,1,-2,0,0,-2,1,1,-1,-1,
-2,0,0,0,0,0,-1,0,1,3,4,8,8,4,4,0,2,-3,0,-1,-1,-1,1,0,0
,-3,1,0,0,0,0,0,0,1,1,-1,4,-10,-40,-2,0,-1,-1,1,1,0,0,0
,0,0,-2,3,1,-1,0,4,-20,0,0,3,0,1,3,-2,1,-4,0,0,8,6,-5,-
2,-3,-1,24,4,16,30,-24,-84,840],[840,84,-24,24,-4,16,30,-6,8
,-5,2,-1,3,20,0,0,3,-1,0,3,4,0,-1,-2,0,0,2,1,-1,-1,1,2,
0,0,0,0,0,1,0,-1,3,4,8,8,-4,-4,0,2,3,0,-1,-1,-1,-1,0,0,
3,1,0,0,0,0,0,0,1,-1,-1,-4,-10,-40,2,0,1,-1,-1,1,0,0,0,0
,0,2,3,-1,1,0,4,20,0,0,3,0,-1,3,-2,-1,4,0,0,8,-6,-5,2,3
,-1,24,-4,16,30,-24,84,840],[972,162,0,36,-6,0,0,0,0,-3,0,0
,0,18,0,0,0,0,0,0,6,2,0,0,-3,0,0,0,-3,1,0,0,-1,1,0,0,0
,-1,0,0,0,0,12,108,2,18,-2,0,0,0,0,0,0,0,1,0,0,0,0,8,0
,0,0,0,0,0,0,6,0,36,0,0,0,1,-3,0,0,0,0,1,-1,0,0,0,-1,0
,0,18,0,0,0,0,0,0,0,0,6,2,-3,0,0,-3,0,0,0,36,-6,0,0,0,
162,972],[972,-162,0,36,6,0,0,0,0,-3,0,0,0,-18,0,0,0,0,0,0
,-6,2,0,0,-3,0,0,0,3,1,0,0,-1,-1,0,0,0,1,0,0,0,0,12,108
,-2,-18,-2,0,0,0,0,0,0,0,1,0,0,0,0,8,0,0,0,0,0,0,0,-6,
0,36,0,0,0,1,3,0,0,0,0,-1,-1,0,0,0,1,0,0,-18,0,0,0,0,0
,0,0,0,-6,2,-3,0,0,-3,0,0,0,36,6,0,0,0,-162,972],[1008,-
252,90,24,-24,8,-36,0,6,3,0,-3,9,12,-4,0,0,0,0,0,4,0,0,-2
,2,-6,0,-1,3,1,-1,0,0,0,-2,1,0,-1,0,3,-9,2,0,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,1,-1
,-3,1,2,-1,0,0,0,0,9,-3,1,0,-2,-12,4,0,0,0,0,0,2,0,-4,0
,-2,-6,0,-3,0,-9,3,-24,24,-8,36,-90,252,-1008],[1008,252,90,
24,24,8,-36,0,6,3,0,-3,-9,-12,-4,0,0,0,0,0,-4,0,0,-2,2,-6
,0,-1,-3,1,1,0,0,0,-2,1,0,1,0,-3,-9,2,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,-1,-1,3,1
,2,-1,0,0,0,0,9,3,-1,0,-2,12,4,0,0,0,0,0,2,0,4,0,-2,-6
,0,-3,0,9,3,-24,-24,-8,36,-90,-252,-1008],[1050,-210,15,50,10
,-10,15,-15,-17,0,1,1,-3,-2,2,4,6,-2,0,2,2,0,0,-1,0,-1,1
,-1,0,0,-1,1,0,0,-1,0,1,0,0,1,-3,-1,-6,58,2,-14,0,3,0,0
,0,-2,4,0,0,0,0,0,-2,-2,1,0,0,-2,0,-2,0,6,7,-30,1,-1,-1
,0,0,-1,-1,0,1,0,0,1,-3,1,0,0,-1,-2,2,4,2,0,-2,6,-1,0,2
,0,0,-17,-15,0,1,-3,1,50,10,-10,15,15,-210,1050],[1050,210,15
,50,-10,-10,15,15,-17,0,-1,1,3,2,2,-4,6,2,0,2,-2,0,0,-1,0
,-1,-1,-1,0,0,1,-1,0,0,-1,0,-1,0,0,-1,-3,-1,-6,58,-2,14,0
,3,0,0,0,-2,4,0,0,0,0,0,-2,-2,1,0,0,-2,0,2,0,-6,7,-30,-
1,-1,1,0,0,-1,-1,0,-1,0,0,-1,-3,-1,0,0,-1,2,2,-4,2,0,2,6
,-1,0,-2,0,0,-17,15,0,-1,3,1,50,-10,-10,15,15,210,1050],[
1134,0,0,-18,0,30,81,0,0,-6,0,0,0,0,6,0,0,0,0,0,0,-2,0,-
3,4,0,0,0,0,2,0,0,0,0,0,1,0,0,0,0,0,0,-18,-18,0,0,-2,-3
,0,0,0,0,0,0,0,0,0,0,-2,-2,-1,0,0,2,0,0,0,0,9,30,0,0,0
,2,0,0,0,1,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,-3,0,0,-2,
4,0,0,-6,0,0,0,-18,0,30,81,0,0,1134],[1296,-216,81,-24,-16,
24,0,9,-9,1,-3,0,0,24,-4,4,0,0,0,0,0,0,1,0,-6,-3,3,0,-1
,-1,0,-1,1,-1,1,0,-1,1,0,0,0,-3,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,-3,3,0,1,1,0,-1,0,1
,1,-1,1,0,0,-1,0,3,-24,4,-4,0,0,0,0,0,-1,0,0,6,9,-9,-1,
3,0,0,24,16,-24,0,-81,216,-1296],[1296,216,81,-24,16,24,0,-9,
-9,1,3,0,0,-24,-4,-4,0,0,0,0,0,0,-1,0,-6,-3,-3,0,1,-1,0,
1,1,1,1,0,1,-1,0,0,0,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,1,0,0,0,0,0,0,-1,0,0,0,0,0,3,3,0,1,-1,0,-1,0,-1,-1,-1
,-1,0,0,1,0,3,24,4,4,0,0,0,0,0,1,0,0,6,9,9,-1,-3,0,0,
24,-16,-24,0,-81,-216,-1296],[1344,336,84,64,16,0,-24,6,4,-1,-
2,-2,-6,16,0,0,-6,-2,0,-2,0,0,1,0,4,4,-2,0,1,-1,0,-2,0,0
,0,1,0,1,0,-2,-6,0,0,64,0,16,0,0,0,0,0,-2,4,0,0,0,0,-1
,0,0,0,0,0,0,-1,-2,0,0,-8,0,-2,4,0,-1,1,0,0,1,0,0,0,-2
,-6,-2,1,0,0,16,0,0,-2,0,-2,-6,0,1,0,0,4,4,6,-1,-2,-6,-2
,64,16,0,-24,84,336,1344],[1344,-336,84,64,-16,0,-24,-6,4,-1,
2,-2,6,-16,0,0,-6,2,0,-2,0,0,-1,0,4,4,2,0,-1,-1,0,2,0,0
,0,1,0,-1,0,2,-6,0,0,64,0,-16,0,0,0,0,0,-2,4,0,0,0,0,-1
,0,0,0,0,0,0,-1,2,0,0,-8,0,2,4,0,-1,-1,0,0,1,0,0,0,2,-
6,2,-1,0,0,-16,0,0,-2,0,2,-6,0,-1,0,0,4,4,-6,-1,2,6,-2,
64,-16,0,-24,84,-336,1344],[1344,0,-60,32,0,32,60,0,12,-6,0,-
6,0,0,0,0,12,0,0,4,0,0,0,-2,6,-4,0,2,0,-2,0,0,0,0,0,0,
0,0,0,0,6,-4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,4,0,2,0,-2,0,0,0,0,0,0,-6,0,0,0,4,0
,0,0,-4,0,0,-12,2,0,0,0,-6,-12,0,6,0,0,6,-32,0,-32,-60,60
,0,-1344],[1400,210,-25,60,-20,-4,-20,15,-15,0,3,3,-3,6,-2,-2
,8,0,2,0,-10,-2,0,-2,0,3,-3,-1,0,0,1,1,0,0,-1,0,-1,0,1,
3,13,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,0,0,0,3,-3,-1,0,0,1,1,0,1,0,0,-1,-13,-3,0,-1,1,-6,
2,2,0,-2,0,-8,2,0,10,2,0,15,-15,0,-3,3,-3,-60,20,4,20,25,
-210,-1400],[1400,-210,-25,60,20,-4,-20,-15,-15,0,-3,3,3,-6,-2
,2,8,0,2,0,10,-2,0,-2,0,3,3,-1,0,0,-1,-1,0,0,-1,0,1,0,1
,-3,13,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,0,0,0,0,-3,-3,1,0,0,1,1,0,-1,0,0,1,-13,3,0,-1,1,6
,2,-2,0,-2,0,-8,2,0,-10,2,0,15,15,0,3,-3,-3,-60,-20,4,20,
25,210,-1400],[1400,-280,50,40,0,-16,50,-10,-6,0,4,0,-10,8,0,
0,5,-1,-1,1,0,0,0,2,0,-2,2,2,0,0,0,0,0,0,0,0,0,0,-1,2,
-4,2,-8,-72,0,16,0,-2,-1,-1,1,-3,-3,1,0,-1,-1,0,0,0,0,0,0
,0,0,1,1,-8,-6,40,2,-2,0,0,0,2,0,0,0,0,0,0,-4,2,0,-1,2
,8,0,0,1,-1,-1,5,2,0,0,0,0,-6,-10,0,4,-10,0,40,0,-16,50,
50,-280,1400],[1400,280,50,40,0,-16,50,10,-6,0,-4,0,10,-8,0,0
,5,1,-1,1,0,0,0,2,0,-2,-2,2,0,0,0,0,0,0,0,0,0,0,-1,-2,
-4,2,-8,-72,0,-16,0,-2,1,1,1,-3,-3,-1,0,1,1,0,0,0,0,0,0,
0,0,-1,1,8,-6,40,-2,-2,0,0,0,2,0,0,0,0,0,0,-4,-2,0,-1,2
,-8,0,0,1,-1,1,5,2,0,0,0,0,-6,10,0,-4,10,0,40,0,-16,50,
50,280,1400],[1400,-350,95,60,-20,-4,40,-5,9,0,3,0,-10,-10,-2
,-2,-4,2,-1,0,6,-2,0,4,0,3,1,2,0,0,0,1,0,0,-1,0,-1,0,1
,-2,4,-1,0,0,0,0,0,0,-2,-1,0,0,0,0,0,1,2,0,0,0,0,0,0,0
,0,0,0,0,0,0,-1,-3,0,0,0,-2,1,0,1,0,0,-1,-4,2,0,-1,1,10
,2,2,0,1,-2,4,-4,0,-6,2,0,-9,5,0,-3,10,0,-60,20,4,-40,-95
,350,-1400],[1400,350,95,60,20,-4,40,5,9,0,-3,0,10,10,-2,2,-
4,-2,-1,0,-6,-2,0,4,0,3,-1,2,0,0,0,-1,0,0,-1,0,1,0,1,2,
4,-1,0,0,0,0,0,0,2,1,0,0,0,0,0,-1,-2,0,0,0,0,0,0,0,0,0
,0,0,0,0,1,-3,0,0,0,-2,1,0,-1,0,0,1,-4,-2,0,-1,1,-10,2,
-2,0,1,2,4,-4,0,6,2,0,-9,-5,0,3,-10,0,-60,-20,4,-40,-95,-
350,-1400],[1400,0,20,-40,0,-8,65,0,4,0,0,-2,0,0,0,0,8,0,2
,-4,0,0,0,1,0,-4,0,-2,0,0,0,0,0,0,0,0,0,0,-1,0,14,4,24
,-8,0,0,0,1,0,0,0,-2,4,0,0,0,0,0,0,-8,1,0,0,4,0,0,-2,0
,1,40,0,-4,0,0,0,-2,0,0,0,0,0,0,14,0,0,-1,4,0,0,0,-4,2
,0,8,1,0,0,0,0,4,0,0,0,0,-2,-40,0,-8,65,20,0,1400],[1575,
-315,90,15,-15,11,-45,0,-6,0,0,-3,9,21,-1,3,0,0,0,0,-7,1,0
,-1,0,-6,0,-1,0,0,-1,0,0,0,2,0,0,0,0,-3,9,2,-9,-57,1,9,
1,3,0,0,0,0,0,0,0,0,0,0,-1,3,-1,1,-1,-1,0,0,0,-3,3,15,0
,-6,-1,0,0,-1,2,0,0,0,0,0,9,-3,0,0,2,21,-1,3,0,0,0,0,-1
,0,-7,1,0,-6,0,0,0,9,-3,15,-15,11,-45,90,-315,1575],[1575,
315,90,15,15,11,-45,0,-6,0,0,-3,-9,-21,-1,-3,0,0,0,0,7,1,0
,-1,0,-6,0,-1,0,0,1,0,0,0,2,0,0,0,0,3,9,2,-9,-57,-1,-9,
1,3,0,0,0,0,0,0,0,0,0,0,-1,3,-1,-1,1,-1,0,0,0,3,3,15,0
,-6,1,0,0,-1,2,0,0,0,0,0,9,3,0,0,2,-21,-1,-3,0,0,0,0,-1
,0,7,1,0,-6,0,0,0,-9,-3,15,15,11,-45,90,315,1575],[1680,0,
60,-80,0,32,6,0,-20,0,0,-2,0,0,0,0,6,0,0,-2,0,0,0,2,-5,4
,0,2,0,0,0,0,0,0,0,1,0,0,0,0,6,-4,16,16,0,0,0,2,0,0,-2
,-2,-2,0,-1,0,0,0,0,0,0,0,0,0,0,0,2,0,-2,-16,0,4,0,0,0
,2,0,1,0,0,0,0,6,0,0,0,-4,0,0,0,-2,0,0,6,2,0,0,0,-5,-
20,0,0,0,0,-2,-80,0,32,6,60,0,1680],[2016,0,-90,48,0,16,36,
0,-6,6,0,-6,0,0,-8,0,0,0,0,0,0,0,0,2,4,6,0,-2,0,2,0,0,
0,0,2,-1,0,0,0,0,-18,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0,-6,0,-2,0,2,-2,1,0,0,0,0,
18,0,0,0,2,0,8,0,0,0,0,0,-2,0,0,0,-4,6,0,-6,0,0,6,-48,0
,-16,-36,90,0,-2016],[2100,-210,75,-60,-10,12,-60,15,-5,0,-1,-
2,6,14,-4,2,3,-1,0,3,6,0,0,0,0,3,-1,0,0,0,0,-1,0,0,-1,0
,-1,0,0,2,-6,3,4,52,-2,-10,0,-4,-3,0,1,1,1,1,0,0,-3,0,0
,-4,0,0,2,0,0,-1,-1,-2,4,20,-1,3,0,0,0,0,-1,0,-1,0,0,-1
,-6,2,0,0,3,14,-4,2,3,0,-1,3,0,0,6,0,0,-5,15,0,-1,6,-2,
-60,-10,12,-60,75,-210,2100],[2100,210,75,-60,10,12,-60,-15,-5,
0,1,-2,-6,-14,-4,-2,3,1,0,3,-6,0,0,0,0,3,1,0,0,0,0,1,0,
0,-1,0,1,0,0,-2,-6,3,4,52,2,10,0,-4,3,0,1,1,1,-1,0,0,3,
0,0,-4,0,0,-2,0,0,1,-1,2,4,20,1,3,0,0,0,0,-1,0,1,0,0,1
,-6,-2,0,0,3,-14,-4,-2,3,0,1,3,0,0,-6,0,0,-5,-15,0,1,-6,
-2,-60,10,12,-60,75,210,2100],[2100,0,30,-60,0,-20,30,0,14,0,
0,2,0,0,-4,0,12,0,0,0,0,0,0,-2,0,6,0,-2,0,0,0,0,0,0,2,
0,0,0,0,0,-6,-2,-12,116,0,0,0,2,0,0,0,2,-4,0,0,0,0,0,-4
,-4,0,0,0,0,0,0,2,0,-10,20,0,6,0,0,0,-2,2,0,0,0,0,0,-6
,0,0,0,-2,0,-4,0,0,0,0,12,-2,0,0,0,0,14,0,0,0,0,2,-60,0
,-20,30,30,0,2100],[2240,-336,-4,64,16,0,-40,-6,-4,-5,-2,2,6
,-16,0,0,2,2,2,-2,0,0,-1,0,0,4,2,0,-1,-1,0,-2,0,0,0,0,0
,1,-1,2,-10,0,0,-64,0,16,0,0,0,0,0,2,-4,0,0,0,0,1,0,0,0
,0,0,0,1,-2,0,0,8,0,2,4,0,-1,-1,0,0,0,0,0,0,-2,-10,2,1
,-1,0,-16,0,0,-2,2,2,2,0,-1,0,0,0,-4,-6,-5,-2,6,2,64,16,
0,-40,-4,-336,2240],[2240,336,-4,64,-16,0,-40,6,-4,-5,2,2,-6,
16,0,0,2,-2,2,-2,0,0,1,0,0,4,-2,0,1,-1,0,2,0,0,0,0,0,-1
,-1,-2,-10,0,0,-64,0,-16,0,0,0,0,0,2,-4,0,0,0,0,1,0,0,0
,0,0,0,1,2,0,0,8,0,-2,4,0,-1,1,0,0,0,0,0,0,2,-10,-2,-1
,-1,0,16,0,0,-2,2,-2,2,0,1,0,0,0,-4,6,-5,2,-6,2,64,-16,0
,-40,-4,336,2240],[2268,378,81,12,10,-12,0,-9,9,-2,-3,0,0,-6
,-4,2,0,0,0,0,-6,0,1,0,3,-3,3,0,-2,2,0,1,0,0,-1,0,-1,0
,0,0,0,-3,12,-36,2,-6,0,0,0,0,0,0,0,0,-1,0,0,1,0,4,0,0
,2,0,1,0,0,-6,0,-36,3,-3,0,2,-2,0,-1,0,-1,0,0,1,0,0,0,0
,-3,-6,-4,2,0,0,0,0,0,1,-6,0,3,9,-9,-2,-3,0,0,12,10,-12,
0,81,378,2268],[2268,-378,81,12,-10,-12,0,9,9,-2,3,0,0,6,-4,
-2,0,0,0,0,6,0,-1,0,3,-3,-3,0,2,2,0,-1,0,0,-1,0,1,0,0,0
,0,-3,12,-36,-2,6,0,0,0,0,0,0,0,0,-1,0,0,1,0,4,0,0,-2,0
,1,0,0,6,0,-36,-3,-3,0,2,2,0,-1,0,1,0,0,-1,0,0,0,0,-3,6
,-4,-2,0,0,0,0,0,-1,6,0,3,9,9,-2,3,0,0,12,-10,-12,0,81,-
378,2268],[2400,120,60,-80,0,16,-60,0,-12,0,0,-3,-3,-24,0,0,6
,0,0,2,8,0,0,2,0,4,0,1,0,0,1,0,-1,-1,0,0,0,0,0,3,3,4,
0,0,0,0,0,0,6,0,0,0,0,0,0,0,-6,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,-4,-1,0,0,-1,0,0,0,1,1,0,-3,-3,0,0,-4,24,0,0,-2,0
,0,-6,-2,0,-8,0,0,12,0,0,0,3,3,80,0,-16,60,-60,-120,-2400]
,[2400,-120,60,-80,0,16,-60,0,-12,0,0,-3,3,24,0,0,6,0,0,2,
-8,0,0,2,0,4,0,1,0,0,-1,0,-1,1,0,0,0,0,0,-3,3,4,0,0,0,
0,0,0,-6,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
-4,1,0,0,-1,0,0,0,-1,1,0,-3,3,0,0,-4,-24,0,0,-2,0,0,-6,-
2,0,8,0,0,12,0,0,0,-3,3,80,0,-16,60,-60,120,-2400],[2688,0,
-48,0,0,0,60,0,-16,8,0,-4,0,0,0,0,-12,0,0,0,0,0,0,0,3,0
,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-12,0,0,128,0,0,0,4,0,0,0
,2,-4,0,-1,0,0,2,0,0,0,0,0,0,2,0,-2,0,-4,64,0,0,0,0,0,
0,0,0,0,0,0,0,-12,0,0,0,0,0,0,0,0,0,0,-12,0,0,0,0,3,-16
,0,8,0,0,-4,0,0,0,60,-48,0,2688],[2800,280,55,-40,0,-24,80,
-5,9,0,-3,0,8,-24,-4,4,10,0,1,-2,0,0,0,0,0,-1,-3,0,0,0,0
,3,0,0,1,0,-1,0,-1,0,8,3,0,0,0,0,0,0,-2,-1,0,0,0,0,0,1
,2,0,0,0,0,0,0,0,0,0,0,0,0,0,3,1,0,0,0,0,-1,0,1,0,0,-
3,-8,0,0,1,-3,24,4,-4,2,-1,0,-10,0,0,0,0,0,-9,5,0,3,-8,0
,40,0,24,-80,-55,-280,-2800],[2800,-280,55,-40,0,-24,80,5,9,0
,3,0,-8,24,-4,-4,10,0,1,-2,0,0,0,0,0,-1,3,0,0,0,0,-3,0,
0,1,0,1,0,-1,0,8,3,0,0,0,0,0,0,2,1,0,0,0,0,0,-1,-2,0,0
,0,0,0,0,0,0,0,0,0,0,0,-3,1,0,0,0,0,-1,0,-1,0,0,3,-8,0
,0,1,-3,-24,4,4,2,-1,0,-10,0,0,0,0,0,-9,-5,0,-3,8,0,40,0
,24,-80,-55,280,-2800],[2835,-189,-81,51,19,3,0,-9,-9,5,-3,0,
0,3,-5,-1,0,0,0,0,3,-1,1,0,0,3,3,0,1,1,0,1,0,0,1,0,-1,
-1,0,0,0,3,3,-45,-5,3,-1,0,0,0,0,0,0,0,0,0,0,-1,-1,3,0,
1,-1,3,-1,0,0,3,0,-45,3,3,0,1,1,0,1,0,-1,0,0,1,0,0,-1,0
,3,3,-5,-1,0,0,0,0,0,1,3,-1,0,-9,-9,5,-3,0,0,51,19,3,0,
-81,-189,2835],[2835,189,-81,51,-19,3,0,9,-9,5,3,0,0,-3,-5,1
,0,0,0,0,-3,-1,-1,0,0,3,-3,0,-1,1,0,-1,0,0,1,0,1,1,0,0
,0,3,3,-45,5,-3,-1,0,0,0,0,0,0,0,0,0,0,-1,-1,3,0,-1,1,3
,-1,0,0,-3,0,-45,-3,3,0,1,-1,0,1,0,1,0,0,-1,0,0,1,0,3,-
3,-5,1,0,0,0,0,0,-1,-3,-1,0,-9,9,5,3,0,0,51,-19,3,0,-81,
189,2835],[3150,0,-90,30,0,22,45,0,6,0,0,-6,0,0,-2,0,0,0,0
,0,0,2,0,1,0,6,0,-2,0,0,0,0,0,0,-2,0,0,0,0,0,18,-2,-18
,-114,0,0,2,-3,0,0,0,0,0,0,0,0,0,0,-2,6,1,0,0,-2,0,0,0
,0,-3,30,0,6,0,0,0,-2,-2,0,0,0,0,0,18,0,0,0,-2,0,-2,0,0
,0,0,0,1,0,0,2,0,6,0,0,0,0,-6,30,0,22,45,-90,0,3150],[
3200,-160,-40,0,0,0,-40,20,8,0,0,2,2,-32,0,0,-4,-2,-1,0,0,
0,0,0,0,0,4,0,0,0,0,0,1,1,0,0,0,0,-1,-2,14,0,0,128,0,0
,0,0,2,-1,0,-4,-4,0,0,-1,2,0,0,0,0,0,0,0,0,0,0,0,8,0,4
,0,0,0,0,0,0,0,0,1,1,0,14,-2,0,-1,0,-32,0,0,0,-1,-2,-4,
0,0,0,0,0,8,20,0,0,2,2,0,0,0,-40,-40,-160,3200],[3200,160,-
40,0,0,0,-40,-20,8,0,0,2,-2,32,0,0,-4,2,-1,0,0,0,0,0,0,0
,-4,0,0,0,0,0,1,-1,0,0,0,0,-1,2,14,0,0,128,0,0,0,0,-2,1
,0,-4,-4,0,0,1,-2,0,0,0,0,0,0,0,0,0,0,0,8,0,-4,0,0,0,0
,0,0,0,0,-1,1,0,14,2,0,-1,0,32,0,0,0,-1,2,-4,0,0,0,0,0
,8,-20,0,0,-2,2,0,0,0,-40,-40,160,3200],[3240,594,81,84,-4,-
12,0,9,-9,-5,-3,0,0,6,2,-2,0,0,0,0,6,2,1,0,0,-3,3,0,-1,
1,0,-1,-1,1,1,0,-1,-1,0,0,0,-3,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,-3,3,0,-1,1,0,-1,0,1
,-1,1,1,0,0,1,0,3,-6,-2,2,0,0,0,0,0,-1,-6,-2,0,9,-9,5,3
,0,0,-84,4,12,0,-81,-594,-3240],[3240,-594,81,84,4,-12,0,-9,-
9,-5,3,0,0,-6,2,2,0,0,0,0,-6,2,-1,0,0,-3,-3,0,1,1,0,1,-
1,-1,1,0,1,1,0,0,0,-3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1
,0,0,0,0,0,0,-1,0,0,0,0,0,3,3,0,-1,-1,0,-1,0,-1,1,1,-1
,0,0,-1,0,3,6,-2,-2,0,0,0,0,0,1,6,-2,0,9,9,5,-3,0,0,-84
,-4,12,0,-81,594,-3240],[3360,336,-6,16,-16,-16,60,6,-18,5,0,
0,6,-16,0,0,-6,2,0,2,0,0,-1,-2,0,-2,-2,2,1,-1,0,-4,0,0,0
,0,0,1,0,-2,-12,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0
,0,0,0,0,1,0,0,0,0,0,2,2,0,1,-1,-2,0,0,0,0,0,4,12,2,-1
,0,-2,16,0,0,-2,0,-2,6,2,1,0,0,0,18,-6,-5,0,-6,0,-16,16,
16,-60,6,-336,-3360],[3360,-336,-6,16,16,-16,60,-6,-18,5,0,0,-
6,16,0,0,-6,-2,0,2,0,0,1,-2,0,-2,2,2,-1,-1,0,4,0,0,0,0,
0,-1,0,2,-12,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,
0,0,0,1,0,0,0,0,0,-2,2,0,1,1,-2,0,0,0,0,0,-4,12,-2,1,0
,-2,-16,0,0,-2,0,2,6,2,-1,0,0,0,18,6,-5,0,6,0,-16,-16,16
,-60,6,336,-3360],[4096,512,64,0,0,0,64,-16,0,-4,0,0,8,0,0,
0,-8,0,1,0,0,0,-1,0,-4,0,0,0,2,0,0,0,1,1,0,-1,0,0,1,0,
-8,0,0,0,0,0,0,0,-4,-1,0,0,0,0,0,-1,-4,-1,0,0,0,0,0,0,-
1,0,0,0,0,0,0,0,0,0,2,0,0,-1,0,1,1,0,-8,0,0,1,0,0,0,0
,0,1,0,-8,0,-1,0,0,-4,0,-16,-4,0,8,0,0,0,0,64,64,512,4096
],[4096,-512,64,0,0,0,64,16,0,-4,0,0,-8,0,0,0,-8,0,1,0,0,
0,1,0,-4,0,0,0,-2,0,0,0,1,-1,0,-1,0,0,1,0,-8,0,0,0,0,0
,0,0,4,1,0,0,0,0,0,1,4,-1,0,0,0,0,0,0,-1,0,0,0,0,0,0,
0,0,0,-2,0,0,-1,0,-1,1,0,-8,0,0,1,0,0,0,0,0,1,0,-8,0,1
,0,0,-4,0,16,-4,0,-8,0,0,0,0,64,64,-512,4096],[4096,-512,64
,0,0,0,-64,16,0,-4,0,0,8,0,0,0,-8,0,1,0,0,0,-1,0,4,0,0
,0,-2,0,0,0,1,1,0,-1,0,0,-1,0,8,0,0,0,0,0,0,0,4,-1,0,0
,0,0,0,1,-4,-1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,2,0,0,
1,0,-1,-1,0,-8,0,0,1,0,0,0,0,0,-1,0,8,0,1,0,0,-4,0,-16,
4,0,-8,0,0,0,0,64,-64,512,-4096],[4096,512,64,0,0,0,-64,-16,
0,-4,0,0,-8,0,0,0,-8,0,1,0,0,0,1,0,4,0,0,0,2,0,0,0,1,-
1,0,-1,0,0,-1,0,8,0,0,0,0,0,0,0,-4,1,0,0,0,0,0,-1,4,-1
,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-2,0,0,1,0,1,-1,0,-8
,0,0,1,0,0,0,0,0,-1,0,8,0,-1,0,0,-4,0,16,4,0,8,0,0,0,0
,64,-64,-512,-4096],[4200,-420,-30,40,20,-8,-30,0,-6,0,-2,3,3
,-4,0,0,-3,-1,0,1,4,0,0,-2,0,-2,-4,1,0,0,1,2,0,0,0,0,0
,0,0,-1,15,-2,8,-24,4,4,0,-2,-3,0,-1,3,3,-1,0,0,-3,0,0,-
8,0,0,0,0,0,1,1,-4,-6,40,-4,-2,1,0,0,1,0,0,0,0,0,2,15,-
1,0,0,-2,-4,0,0,1,0,-1,-3,-2,0,4,0,0,-6,0,0,-2,3,3,40,20
,-8,-30,-30,-420,4200],[4200,420,-30,40,-20,-8,-30,0,-6,0,2,3
,-3,4,0,0,-3,1,0,1,-4,0,0,-2,0,-2,4,1,0,0,-1,-2,0,0,0,0
,0,0,0,1,15,-2,8,-24,-4,-4,0,-2,3,0,-1,3,3,1,0,0,3,0,0,
-8,0,0,0,0,0,-1,1,4,-6,40,4,-2,-1,0,0,1,0,0,0,0,0,-2,15
,1,0,0,-2,4,0,0,1,0,1,-3,-2,0,-4,0,0,-6,0,0,2,-3,3,40,-
20,-8,-30,-30,420,4200],[4200,210,-75,20,-20,4,-60,15,3,0,3,3
,-3,-26,2,2,6,-2,0,-2,-2,-2,0,2,0,-7,5,1,0,0,-1,1,0,0,1
,0,1,0,0,-1,-15,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,-5,7,1,0,0,-1,-1,0,-1,0,0,-1,15,1
,0,0,-1,26,-2,-2,2,0,2,-6,-2,0,2,2,0,-3,-15,0,-3,3,-3,-20
,20,-4,60,75,-210,-4200],[4200,-210,-75,20,20,4,-60,-15,3,0,-3
,3,3,26,2,-2,6,2,0,-2,2,-2,0,2,0,-7,-5,1,0,0,1,-1,0,0,1
,0,-1,0,0,1,-15,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,5,7,-1,0,0,-1,-1,0,1,0,0,1,15,-1,
0,0,-1,-26,-2,2,2,0,-2,-6,-2,0,-2,2,0,-3,15,0,3,-3,-3,-20
,-20,-4,60,75,210,-4200],[4200,0,-120,40,0,24,15,0,8,0,0,-4,
0,0,0,0,6,0,0,-2,0,0,0,3,0,-8,0,0,0,0,0,0,0,0,0,0,0,0
,0,0,-12,0,8,104,0,0,0,-1,0,0,2,2,2,0,0,0,0,0,0,-8,-1,0
,0,-4,0,0,2,0,-1,-40,0,-8,0,0,0,0,0,0,0,0,0,0,-12,0,0,0
,0,0,0,0,-2,0,0,6,3,0,0,0,0,8,0,0,0,0,-4,40,0,24,15,-
120,0,4200],[4480,0,-80,0,0,0,-44,0,16,0,0,4,0,0,0,0,4,0,-
2,0,0,0,0,0,-5,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,-20,0,0,-
128,0,0,0,4,0,0,0,-2,4,0,-1,0,0,0,0,0,0,0,0,0,0,0,-2,0
,4,64,0,0,0,0,0,0,0,1,0,0,0,0,-20,0,0,1,0,0,0,0,0,-2,0
,4,0,0,0,0,-5,16,0,0,0,0,4,0,0,0,-44,-80,0,4480],[4536,-
378,-81,60,20,12,0,9,9,-4,-3,0,0,-30,-2,-2,0,0,0,0,-6,2,1
,0,-6,3,3,0,2,0,0,-1,0,0,-1,0,-1,0,0,0,0,3,0,0,0,0,0,0
,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,-3,-3,0
,0,-2,0,1,0,1,0,0,1,0,0,0,0,-3,30,2,2,0,0,0,0,0,-1,6,-
2,6,-9,-9,4,3,0,0,-60,-20,-12,0,81,378,-4536],[4536,378,-81,
60,-20,12,0,-9,9,-4,3,0,0,30,-2,2,0,0,0,0,6,2,-1,0,-6,3,
-3,0,-2,0,0,1,0,0,-1,0,1,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0
,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,3,-3,0,0,2,0,1,
0,-1,0,0,-1,0,0,0,0,-3,-30,2,-2,0,0,0,0,0,1,-6,-2,6,-9,9
,4,-3,0,0,-60,20,-12,0,81,-378,-4536],[4536,0,0,-72,0,-24,81
,0,0,6,0,0,0,0,0,0,0,0,0,0,0,0,0,-3,1,0,0,0,0,-2,0,0,
0,0,0,1,0,0,0,0,0,0,24,-72,0,0,0,-3,0,0,0,0,0,0,1,0,0,
0,0,8,-1,0,0,-4,0,0,0,0,9,-24,0,0,0,-2,0,0,0,1,0,0,0,0
,0,0,0,0,0,0,0,0,0,0,0,0,-3,0,0,0,1,0,0,6,0,0,0,-72,0
,-24,81,0,0,4536],[5600,0,-10,-80,0,16,100,0,-6,0,0,-6,0,0,
8,0,-4,0,2,-4,0,0,0,2,0,-2,0,-2,0,0,0,0,0,0,-2,0,0,0,1
,0,-2,-2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
0,0,0,0,0,0,0,2,0,0,0,2,2,0,0,0,0,0,2,0,0,-1,2,0,-8,0
,4,-2,0,4,-2,0,0,0,0,6,0,0,0,0,6,80,0,-16,-100,10,0,-5600
],[5600,280,20,-80,0,-16,-20,-20,12,0,0,3,-1,8,0,0,2,2,-1,2
,-8,0,0,-2,0,4,4,-1,0,0,-1,0,0,0,0,0,0,0,1,1,-11,-4,0,0
,0,0,0,0,4,-1,0,0,0,0,0,1,-4,0,0,0,0,0,0,0,0,0,0,0,0,
0,-4,-4,1,0,0,1,0,0,0,0,0,0,11,-1,0,-1,4,-8,0,0,-2,1,-2
,-2,2,0,8,0,0,-12,20,0,0,1,-3,80,0,16,20,-20,-280,-5600],[
5600,-280,20,-80,0,-16,-20,20,12,0,0,3,1,-8,0,0,2,-2,-1,2,8
,0,0,-2,0,4,-4,-1,0,0,1,0,0,0,0,0,0,0,1,-1,-11,-4,0,0,0
,0,0,0,-4,1,0,0,0,0,0,-1,4,0,0,0,0,0,0,0,0,0,0,0,0,0,
4,-4,-1,0,0,1,0,0,0,0,0,0,11,1,0,-1,4,8,0,0,-2,1,2,-2,2
,0,-8,0,0,-12,-20,0,0,-1,-3,80,0,16,20,-20,280,-5600],[5670,
0,0,-90,0,6,-81,0,0,0,0,0,0,0,6,0,0,0,0,0,0,-2,0,3,5,0
,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,6,-90,0,0,-2,3,0,0,0
,0,0,0,1,0,0,0,-2,6,1,0,0,-2,0,0,0,0,-9,6,0,0,0,0,0,0
,0,-1,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,3,0,0,-2,5,0,0,
0,0,0,0,-90,0,6,-81,0,0,5670],[6075,-405,0,-45,15,-9,0,0,0,
0,0,0,0,27,3,-3,0,0,0,0,-9,1,0,0,0,0,0,0,0,0,0,0,-1,1,
0,0,0,0,0,0,0,0,-21,27,-1,-9,1,0,0,0,0,0,0,0,0,0,0,0,3
,-1,0,-1,1,3,0,0,0,3,0,-45,0,0,0,0,0,0,0,0,0,1,-1,0,0,
0,0,0,0,27,3,-3,0,0,0,0,0,0,-9,1,0,0,0,0,0,0,0,-45,15,-
9,0,0,-405,6075],[6075,405,0,-45,-15,-9,0,0,0,0,0,0,0,-27,3
,3,0,0,0,0,9,1,0,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,
0,0,-21,27,1,9,1,0,0,0,0,0,0,0,0,0,0,0,3,-1,0,1,-1,3,0
,0,0,-3,0,-45,0,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,-27,3
,3,0,0,0,0,0,0,9,1,0,0,0,0,0,0,0,-45,-15,-9,0,0,405,6075
],[7168,0,-128,0,0,0,-16,0,0,8,0,0,0,0,0,0,-8,0,-2,0,0,0
,0,0,2,0,0,0,0,0,0,0,0,0,0,1,0,0,-1,0,32,0,0,0,0,0,0,
0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,-2,0,0,0,0,0,0,0,0,
0,0,0,0,-1,0,0,0,0,-32,0,0,1,0,0,0,0,0,2,0,8,0,0,0,0,-
2,0,0,-8,0,0,0,0,0,0,16,128,0,-7168]]:
#
`coxeter/class_size/F4` := [1,12,12,36,16,32,32,96,96
,72,18,96,72,72,144,12,96,96,32,32,12,12,36,16,1]:
`coxeter/class_rep/F4` := [[],[1],[3],[2,3],[1,2,1,3,2,3,4,3],[1,
2],[3,4],[1,2,3],[2,3,4],[1,3],[2,3,2,3],[1,2,3,4],[1,2,3
,2,3],[2,3,2,3,4],[1,2,3,2,3,4],[1,2,1,3,2,1,3,4,3,2,3
,4],[1,2,4],[1,3,4],[1,2,3,2,3,4,3,2,3,4],[1,2,1,3,2,1,
3,2,3,4],[2,3,2,3,4,3,2,3,4],[1,2,1,3,2,1,3,2,3],[1,2,1
,3,2,1,3,2,3,4,3,2,3,4],[1,2,1,3,2,1,3,4,3,2,1,3,2,3,
4,3],[1,2,1,3,2,1,3,2,3,4,3,2,1,3,2,3,4,3,2,1,3,2,3,4]
]:
`coxeter/mytype/F4` := [1,-`2`^3,`2`^4*`-1`,-`4`^2*`-2`,-`-3`^4,-`3`^3,
`3`^4,-`6`*`2`*`-3`,-`6`*`-3`^2,-`2`^5*`-1`,`2`^4*`-1`^2,-`-6`^2,-`4`^2*
`-2`*`-1`,`4`^2*`2`*`-2`,-`-4`^3,`-2`^6,-`2`*`3`*`6`*`-1`,-`3`^2*`6`,-
`-3`^3*`-1`^3,`-3`^4,-`2`^3*`-1`^6,`2`^4*`-1`^3,-`4`^2*`-2`*`-1`^2,-`3`^4
,`-1`^12]:
`coxeter/bpermrep/F4` := {s1 = [[4,5],[6,7],[8,10]],s2 = [[3,4],[7,9]
,[10,12]],s3 = [[2,3],[4,6],[5,7],[9,11],[-12]],s4 = [[1,2],[6,8
],[7,10],[9,12],[-11]]},12:
`coxeter/irr_chars/F4` := [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
1,1,1,1,1,1,1,1],[1,-1,-1,1,1,1,1,-1,-1,1,1,1,-1,-1,1,1
,-1,-1,1,1,-1,-1,1,1,1],[1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,-
1,-1,1,-1,1,1,1,1,-1,-1,1,1],[1,-1,1,-1,1,1,1,1,-1,-1,1,
1,-1,1,-1,1,1,-1,1,1,-1,1,-1,1,1],[2,-2,0,0,-1,2,-1,0,1,
0,2,-1,-2,0,0,2,0,1,2,-1,-2,0,0,-1,2],[2,2,0,0,-1,2,-1,0
,-1,0,2,-1,2,0,0,2,0,-1,2,-1,2,0,0,-1,2],[2,0,-2,0,-1,-1
,2,1,0,0,2,-1,0,-2,0,2,1,0,-1,2,0,-2,0,-1,2],[2,0,2,0,-
1,-1,2,-1,0,0,2,-1,0,2,0,2,-1,0,-1,2,0,2,0,-1,2],[4,-2,-
2,2,2,1,1,-1,-1,0,0,0,0,0,0,0,1,1,-1,-1,2,2,-2,-2,-4],[4
,-2,2,-2,2,1,1,1,-1,0,0,0,0,0,0,0,-1,1,-1,-1,2,-2,2,-2,
-4],[4,2,-2,-2,2,1,1,-1,1,0,0,0,0,0,0,0,1,-1,-1,-1,-2,2,
2,-2,-4],[4,2,2,2,2,1,1,1,1,0,0,0,0,0,0,0,-1,-1,-1,-1,-2
,-2,-2,-2,-4],[4,0,0,0,1,-2,-2,0,0,0,4,1,0,0,0,4,0,0,-2
,-2,0,0,0,1,4],[6,0,0,-2,3,0,0,0,0,2,-2,-1,0,0,0,2,0,0
,0,0,0,0,-2,3,6],[6,0,0,2,3,0,0,0,0,-2,-2,-1,0,0,0,2,0
,0,0,0,0,0,2,3,6],[8,4,0,0,-2,2,-1,0,-1,0,0,0,0,0,0,0,
0,1,-2,1,-4,0,0,2,-8],[8,-4,0,0,-2,2,-1,0,1,0,0,0,0,0,0
,0,0,-1,-2,1,4,0,0,2,-8],[8,0,-4,0,-2,-1,2,1,0,0,0,0,0,
0,0,0,-1,0,1,-2,0,4,0,2,-8],[8,0,4,0,-2,-1,2,-1,0,0,0,0
,0,0,0,0,1,0,1,-2,0,-4,0,2,-8],[9,3,-3,-1,0,0,0,0,0,-1,
1,0,-1,1,1,-3,0,0,0,0,3,-3,-1,0,9],[9,3,3,1,0,0,0,0,0,1
,1,0,-1,-1,-1,-3,0,0,0,0,3,3,1,0,9],[9,-3,-3,1,0,0,0,0,
0,1,1,0,1,1,-1,-3,0,0,0,0,-3,-3,1,0,9],[9,-3,3,-1,0,0,0
,0,0,-1,1,0,1,-1,1,-3,0,0,0,0,-3,3,-1,0,9],[12,0,0,0,-3
,0,0,0,0,0,-4,1,0,0,0,4,0,0,0,0,0,0,0,-3,12],[16,0,0,0
,2,-2,-2,0,0,0,0,0,0,0,0,0,0,0,2,2,0,0,0,-2,-16]]:
#
`coxeter/class_size/H3` := [1,12,15,12,20,20,12,15,12,1]:
`coxeter/class_rep/H3` := [[],[1,2],[1],[1,2,3],[2,3],[1,2,1,2,3]
,[1,2,1,2],[1,3],[1,2,1,2,3,2,1,2,3],[1,2,1,2,1,3,2,1,2
,1,3,2,1,2,3]]:
`coxeter/mytype/H3` := [-16,-33,32,41,-52,56,-59,-64,95,128]:
H3 := [[1,1,1,1,1,1,1,1,1,1],[1,1,-1,-1,1,-1,1,1,-1,-1],[3
,-A1+1,-1,A1,0,0,A1,-1,-A1+1,3],[3,-A1+1,1,-A1,0,0,A1,-1,A1-1,-3],[
3,A1,-1,-A1+1,0,0,-A1+1,-1,A1,3],[3,A1,1,A1-1,0,0,-A1+1,-1,-A1,-3],
[4,-1,0,-1,1,1,-1,0,-1,4],[4,-1,0,1,1,-1,-1,0,1,-4],[5,0,
1,0,-1,-1,0,1,0,5],[5,0,-1,0,-1,1,0,1,0,-5]]:
`coxeter/irr_chars/H3`:=subs(A1=(sqrt(5)+1)/2,H3): H3:='H3':
#
`coxeter/class_size/H4` := [1,24,144,60,40,720,480,24,400,1200,288,
480,720,144,450,720,1200,1800,720,60,144,480,720,1200,400,288,24
,720,480,60,40,144,24,1]:
`coxeter/class_rep/H4` := [[],[1,2,1,2,1,3,2,1,2,3,4,3],[1,2],[1
],[1,2,1,2,1,3,2,1,2,1,3,2,1,4,3,2,1,2,3,4],[1,2,3],[1
,2,1,2,3,2,4,3],[1,2,1,2,1,3,2,1,2,1,3,4,3,2,1,2,1,3,
2,1,2,3,4,3],[2,3],[1,2,1,2,3],[1,2,1,2,1,3,2,1,2,1,3,2
,4,3],[1,2,3,4],[1,2,4],[1,2,1,2],[1,3],[1,2,1,2,3,4],[1
,2,1,2,3,2,1,2,3,4],[2,3,4],[1,2,1,2,1,3,2,1,2,1,3,4,3
,2,1,2,3,4],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,4,3,2,1,2,1,
3,2,1,4,3,2,1,2,3,4],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,3,4]
,[1,2,1,2,3,2,1,2,3,4,3,2,1,2,3,4],[1,2,1,2,3,2,1,2,3]
,[1,3,4],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,3,4]
,[1,2,1,2,1,3,2,1,2,1,3,2,1,4,3,2,1,2,1,3,2,1,2,3,4,3
],[1,2,1,2,1,3,2,1,2,1,3,2,4,3,2,1,2,1,3,2,1,2,3,4,3,
2,1,2,1,3,2,1,2,3,4,3],[1,2,1,2,4],[1,2,1,2,1,3,2,1,2,
3,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,3,4],[1,2,1,2,1,3,2
,1,2,1,3,2,1,2,3],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,4,3,2,
1,2,1,3,2,1,2,4,3,2,1,2,1,3,2,1,4,3,2,1,2,3,4],[1,2,1
,2,1,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,1,3,2,1,2,3,4,3,2
,1,2,1,3,2,1,2,3,4],[1,2,1,2,1,3,2,1,2,1,3,2,1,4,3,2,
1,2,1,3,2,1,2,3,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,1,3,2
,1,2,3,4,3],[1,2,1,2,1,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,
1,3,2,1,2,3,4,3,2,1,2,1,3,2,1,2,3,4,3,2,1,2,1,3,2,1,2
,3,4,3,2,1,2,1,3,2,1,2,3,4]]:
`coxeter/mytype/H4` := [32,53,65,-64,98,-82,96,133,104,-112,122,116
,-130,119,128,135,146,-160,175,200,165,190,-190,-208,224,242,281
,-238,290,-256,338,379,441,512]:
H4 := [[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
,1,1,1,1,1,1,1,1,1,1],[1,1,1,-1,1,-1,1,1,1,-1,1,1,-1,1
,1,1,1,-1,1,1,1,1,-1,-1,1,1,1,-1,1,-1,1,1,1,1],[4,-2*A1+2
,-A1+2,-2,2,A1-1,-A1+1,-2*A1,1,-1,1,-A1,A1,A1+1,0,0,0,0,0,0,-A1-1,A1
,-A1,1,-1,-1,2*A1,-A1+1,A1-1,2,-2,A1-2,2*A1-2,-4],[4,-2*A1+2,-A1+2,2,2
,-A1+1,-A1+1,-2*A1,1,1,1,-A1,-A1,A1+1,0,0,0,0,0,0,-A1-1,A1,A1,-1,-1
,-1,2*A1,A1-1,A1-1,-2,-2,A1-2,2*A1-2,-4],[4,2*A1,A1+1,2,2,A1,A1,2*A1-2,
1,1,1,A1-1,A1-1,-A1+2,0,0,0,0,0,0,A1-2,-A1+1,-A1+1,-1,-1,-1,-2*A1+2
,-A1,-A1,-2,-2,-A1-1,-2*A1,-4],[4,2*A1,A1+1,-2,2,-A1,A1,2*A1-2,1,-1,1
,A1-1,-A1+1,-A1+2,0,0,0,0,0,0,A1-2,-A1+1,A1-1,1,-1,-1,-2*A1+2,A1,-A1,
2,-2,-A1-1,-2*A1,-4],[6,-A1+4,-2*A1+2,0,3,0,-A1+1,A1+3,0,0,1,A1,0,2*
A1,-2,A1-1,-1,0,-A1,2,2*A1,A1,0,0,0,1,A1+3,0,-A1+1,0,3,-2*A1+2,-A1+4
,6],[6,A1+3,2*A1,0,3,0,A1,-A1+4,0,0,1,-A1+1,0,-2*A1+2,-2,-A1,-1,0,
A1-1,2,-2*A1+2,-A1+1,0,0,0,1,-A1+4,0,A1,0,3,2*A1,A1+3,6],[8,3,-2,0
,5,0,0,3,2,0,-2,0,0,-2,0,-1,1,0,-1,4,-2,0,0,0,2,-2,3,0
,0,0,5,-2,3,8],[8,2,-2,0,4,0,1,-2,2,0,-3,-1,0,-2,0,0,0,
0,0,0,2,1,0,0,-2,3,2,0,-1,0,-4,2,-2,-8],[9,-3*A1+3,-A1+2,3,0
,-A1+1,0,3*A1,0,0,-1,0,A1,A1+1,1,-A1,0,-1,A1-1,-3,A1+1,0,A1,0,0,-1
,3*A1,-A1+1,0,3,0,-A1+2,-3*A1+3,9],[9,-3*A1+3,-A1+2,-3,0,A1-1,0,3*A1,0
,0,-1,0,-A1,A1+1,1,-A1,0,1,A1-1,-3,A1+1,0,-A1,0,0,-1,3*A1,A1-1,0,-
3,0,-A1+2,-3*A1+3,9],[9,3*A1,A1+1,-3,0,-A1,0,-3*A1+3,0,0,-1,0,A1-1,-
A1+2,1,A1-1,0,1,-A1,-3,-A1+2,0,A1-1,0,0,-1,-3*A1+3,-A1,0,-3,0,A1+1,3
*A1,9],[9,3*A1,A1+1,3,0,A1,0,-3*A1+3,0,0,-1,0,-A1+1,-A1+2,1,A1-1,0,-
1,-A1,-3,-A1+2,0,-A1+1,0,0,-1,-3*A1+3,A1,0,3,0,A1+1,3*A1,9],[10,5,0
,0,4,0,-1,5,-2,0,0,-1,0,0,2,1,0,0,1,6,0,-1,0,0,-2,0,5,
0,-1,0,4,0,5,10],[16,-4,1,4,4,-1,-1,-4,1,1,1,-1,-1,1,0,0
,0,0,0,0,1,-1,-1,1,1,1,-4,-1,-1,4,4,1,-4,16],[16,-4,1,-4
,4,1,-1,-4,1,-1,1,-1,1,1,0,0,0,0,0,0,1,-1,1,-1,1,1,-4,1
,-1,-4,4,1,-4,16],[16,4,1,-4,-4,-1,-1,-4,1,1,-1,1,1,1,0,0
,0,0,0,0,-1,-1,-1,-1,-1,1,4,1,1,4,4,-1,-4,-16],[16,4,1,4
,-4,1,-1,-4,1,-1,-1,1,-1,1,0,0,0,0,0,0,-1,-1,1,1,-1,1,4
,-1,1,-4,4,-1,-4,-16],[16,2+4*A1,2*A1,0,2,0,-A1+1,4*A1-6,-2,0,-1
,-A1,0,-2*A1+2,0,0,0,0,0,0,2*A1-2,A1,0,0,2,1,-4*A1+6,0,A1-1,0,-2
,-2*A1,-2-4*A1,-16],[16,-4*A1+6,-2*A1+2,0,2,0,A1,-2-4*A1,-2,0,-1,A1-1,
0,2*A1,0,0,0,0,0,0,-2*A1,-A1+1,0,0,2,1,2+4*A1,0,-A1,0,-2,2*A1-2,4
*A1-6,-16],[18,3,-2,0,0,0,0,3,0,0,3,0,0,-2,2,-1,0,0,-1,-6,
-2,0,0,0,0,3,3,0,0,0,0,-2,3,18],[24,1-4*A1,2*A1-2,0,3,0,-A1+1
,-3+4*A1,0,0,-1,A1,0,-2*A1,0,1,-1,0,1,-4,-2*A1,A1,0,0,0,-1,-3+4*
A1,0,-A1+1,0,3,2*A1-2,1-4*A1,24],[24,-3+4*A1,-2*A1,0,3,0,A1,1-4*A1,0,0
,-1,-A1+1,0,2*A1-2,0,1,-1,0,1,-4,2*A1-2,-A1+1,0,0,0,-1,1-4*A1,0,A1
,0,3,-2*A1,-3+4*A1,24],[24,4-6*A1,2*A1-2,0,6,0,-1,2-6*A1,0,0,1,1,0
,-2*A1,0,0,0,0,0,0,2*A1,-1,0,0,0,-1,-2+6*A1,0,1,0,-6,-2*A1+2,-4
+6*A1,-24],[24,-2+6*A1,-2*A1,0,6,0,-1,-4+6*A1,0,0,1,1,0,2*A1-2,0,0
,0,0,0,0,-2*A1+2,-1,0,0,0,-1,4-6*A1,0,1,0,-6,2*A1,2-6*A1,-24],[
25,0,0,5,-5,0,0,0,1,-1,0,0,0,0,1,0,-1,1,0,5,0,0,0,-1,1
,0,0,0,0,5,-5,0,0,25],[25,0,0,-5,-5,0,0,0,1,1,0,0,0,0,1
,0,-1,-1,0,5,0,0,0,1,1,0,0,0,0,-5,-5,0,0,25],[30,-5*A1+5,0
,0,-3,0,A1-1,5*A1,0,0,0,-A1,0,0,-2,A1,1,0,-A1+1,-2,0,-A1,0,0,0
,0,5*A1,0,A1-1,0,-3,0,-5*A1+5,30],[30,5*A1,0,0,-3,0,-A1,-5*A1+5,0,
0,0,A1-1,0,0,-2,-A1+1,1,0,A1,-2,0,A1-1,0,0,0,0,-5*A1+5,0,-A1,0,-
3,0,5*A1,30],[36,-6,1,-6,0,1,0,6,0,0,-1,0,-1,1,0,0,0,0,0,
0,-1,0,1,0,0,1,-6,-1,0,6,0,-1,6,-36],[36,-6,1,6,0,-1,0,6
,0,0,-1,0,1,1,0,0,0,0,0,0,-1,0,-1,0,0,1,-6,1,0,-6,0,-1
,6,-36],[40,-5,0,0,1,0,1,-5,-2,0,0,1,0,0,0,-1,1,0,-1,4,0
,1,0,0,-2,0,-5,0,1,0,1,0,-5,40],[48,2,-2,0,-6,0,1,-2,0,0
,2,-1,0,-2,0,0,0,0,0,0,2,1,0,0,0,-2,2,0,-1,0,6,2,-2,-48
]]:
`coxeter/irr_chars/H4`:=subs(A1=(sqrt(5)+1)/2,H4): H4:='H4':
#
map(x->assign(coxeter[x],cat(`coxeter/`,x)),
 ['('base')', '('cartan_matrix')', '('char_poly')', '('class_rep')',
  '('class_size')', '('co_base')', '('cox_matrix')', '('cox_number')',
  '('cprod')', '('degrees')', '('diagram')', '('exponents')',
  '('highest_root')', '('index')', '('induce')', '('interior_pt')',
  '('iprod')', '('irr_chars')', '('length_gf')', '('longest_elt')',
  '('multperm')', '('name_of')','('names_of')', '('num_refl')', '('orbit')',
  '('orbit_size')', '('perm2word')', '('perm_char')', '('perm_rep')',
  '('pos_roots')', '('presentation')', '('rank')', '('reduce')',
  '('reflect')', '('restrict')', '('root_coords')', '('size')',
  '('stab_chain')', '('vec2fc')']):
map(x->assign(weyl[x],cat(`weyl/`,x)),
 ['('branch')', '('co_rho')', '('minuscule')', '('rho')', '('tensor')',
  '('toM')', '('toX')', '('weight_coords')', '('weight_mults')',
  '('weight_sys')', '('weights')', '('weyl_dim')']):
if [op(I)]=[1] then # we are in maple6 and lprint is broken
  printf(`%0.70s\n`, `coxeter and weyl 2.4v loaded.`);
  printf(`%0.70s\n`,
         `Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.`)
else
  lprint(`coxeter and weyl 2.4v loaded.`);
  lprint(`Run 'withcoxeter()' or 'withweyl()' to use abbreviated names.`)
fi;
