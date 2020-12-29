function [Ainvfun] = genAinvFunLocal3D(v, reach, spread, allvmode, boxsize, h, n, nlocal, sigma)

if allvmode==1
    numelallv=n/4;
elseif allvmode==2
    numelallv=n^2;
end

v = (v - sigma);

%t_stencil_setup=tic;
[allv,stencil]=allvstencil3D(v,h,nlocal,numelallv,reach,spread);
%t_stencil_setup=toc(t_stencil_setup)

%t_matrix_setup=tic;
[P,C]=buildPC3D(v,h,boxsize,allv,stencil,reach,spread);
%t_matrix_setup=toc(t_matrix_setup)

%tsetup=tic;
%invP=setupinvP3D(P,n,boxsize);
PMF=Multifrontal(P);
%tsetup=toc(tsetup)

funcA=@(u)applyA3D(u,v,h);
%funcM=@(f)solveinvP3D(invP,C*f);
funcM=@(f)PMF\(C*f);

restart=100;
tol=1e-14;
maxit=100;

Ainvfun = @(f)gmres_multicol(funcA,funcM,f,restart,tol,maxit);

end