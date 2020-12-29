close all;
clear all;

warning('off','all');

load test.mat;

opt = [];
opt.verbose = true;
opt.nc = 20;
opt.r = [8,8];
optEqnSol.restart = 2; optEqnSol.maxit = 50; optEqnSol.tol = 1e-15;
GenAinvfun = @(sigma) AinvfuncRPA(Afun,sigma,diagA,n1,n2,optEqnSol);
% GenAinvfun(shift) returns a function handle equivalent to
% (Afun-shift*Bfun)^{-1}

%a = [a1+(b1-a1)/4,b1-(b1-a1)/4]; b = [a2+(b2-a2)/4,b2-(b2-a2)/4];
[eVecOut,eValOut,relerrs] = zologeigsdense(Afun,GenAinvfun,Bfun,n1+n2,a,b,opt);
%[eVecOut,eValOut,relerrs] = zologeigsdense(Afun,GenAinvfun,Bfun,n1+n2,SIGMA1,SIGMA2,opt);
AV = zeros(size(eVecOut));
for cntt = 1:size(eVecOut,2)
    AV(:,cntt) = Afun(eVecOut(:,cntt));
end
fprintf('error of zolo is %1.3e\n',norm(AV-Bfun(eVecOut)*eValOut));










