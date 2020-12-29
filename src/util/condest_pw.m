function [condNum,norm_1_H,norm_1_H_inv] = condest_pw(H,That,V,isSP)
% estimate the condition number of the function handle by planewave
% discretization

opt.restart = 40;
opt.maxit = 50;
opt.tol = 1e-13;
opt.isSP = isSP;


%% estimate the l1-norm of H
% n is the number of grid points per dimension
n = sqrt(numel(V));
x = randn(n^2,1);
x = x/sum(abs(x));

maxIter = 5;
numIter = 0;
while (1)
    numIter = numIter + 1;
    y = H(x);
    xi = ones(size(y));
    pos = find(y<0);
    xi(pos) = -1;
    z = H(xi);
    [val,idx] = max(abs(z));
    if val<=transpose(z)*x
        isNewy = 0;
        break;
    else
        x = zeros(n^2,1);
        x(idx) = 1;
    end
    if numIter>maxIter
        isNewy = 1;
        break;
    end
end
if isNewy
    y = H(x);
end
norm_1_H = sum(abs(y));

%% estimate the l1-norm of H^{-1}
% n is the number of grid points per dimension
n = sqrt(numel(V));
x = randn(n^2,1);
x = x/sum(abs(x));

maxIter = 5;
numIter = 0;
while (1)
    numIter = numIter + 1;
    [y,opt] = solveEqn(That,V,n,0,x,opt);
    xi = ones(size(y));
    pos = find(y<0);
    xi(pos) = -1;
    [z,opt] = solveEqn(That,V,n,0,xi,opt);
    [val,idx] = max(abs(z));
    if val<=transpose(z)*x
        isNewy = 0;
        break;
    else
        x = zeros(n^2,1);
        x(idx) = 1;
    end
    if numIter>maxIter
        isNewy = 1;
        break;
    end
end
if isNewy
    [y,opt] = solveEqn(That,V,n,0,x,opt);
end
norm_1_H_inv = sum(abs(y));

condNum = norm_1_H*norm_1_H_inv;

