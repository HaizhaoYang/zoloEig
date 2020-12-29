function [Ainvfun] = genAinvFunLocal2D(v, reach, spread, allvmode, boxsize, h, n, nlocal, sigma, ecut)

if nargin < 10
    ecut = Inf;
end

if allvmode==1
    numelallv=n/16;
elseif allvmode==2
    numelallv=n^2;
end

v = v - sigma;

%t_stencil_setup=tic;
[allv,stencil]=allvstencil2D(v,h,nlocal,numelallv,reach,spread,ecut);
%t_stencil_setup=toc(t_stencil_setup)

%t_matrix_setup=tic;
[P,C]=buildPC2D(v,h,boxsize,allv,stencil,reach,spread);
%t_matrix_setup=toc(t_matrix_setup)

%tsetup=tic;
PMF=Multifrontal(P);
%invP=setupinvP2D(P,n,boxsize);
%tsetup=toc(tsetup)

funcA=@(u)applyA2D(u,v,h,ecut);
funcM=@(f)PMF\(C*f);

restart = 100;
tol = 1e-14;
maxit = 100;

Ainvfun = @(f)gmres_multicol(funcA,funcM,f,restart,tol,maxit);

end