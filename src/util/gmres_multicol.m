function [x,flag,relres,j] = gmres_multicol ...
    (Afun0,preAfun0,b,restart,tol,maxiter,x0)

% solves the systems (A-omega I)x = rhs
% ||(A-omega I)x - rhs||\(||rhs||) < e for all omega

% INPUTS
% A:     a function handle to apply a NxN matrix
% rhs:   a Nx1 vector
% omega: a vector containing the values of omega
% x0:    initial guess of the solution

%

% OUTPUTS
% x:       a NxNomega vector representing the solution, col(k) is for omega(k)
% iter:    The number of iterations needed, will be less than maxiter
% err:     ||(A-omega I)x - rhs||\(||n||), for each shift. err(k) is for omega(k)


if nargin < 7 || isempty(x0)
    x0 = zeros(size(b));
end

if nargin < 6 || isempty(maxiter)
    maxiter = 100;
end

if nargin < 5 || isempty(tol)
    tol = 1e-6;
end

if nargin < 4 || isempty(restart)
    restart = 30;
end

norms = @(X)sqrt(sum(abs(X).^2,1));

[N,Ncol] = size(b);
Nrst = restart;

if ~isa(Afun0,'function_handle')
    A = Afun0;
    Afun = @(x)A*reshape(x,N,[]);
else
    Afun = @(x)Afun0(reshape(x,N,[]));
end

if isempty(preAfun0)
    preAfun0 = @(x)x;
end

normb = norms(b);

converged = zeros(Ncol,1);
x = zeros(N,Ncol);
relres = zeros(Ncol,1);

Q = zeros(N,Ncol,Nrst);
H = zeros(Nrst+1,Nrst,Ncol);
vh10 = zeros(Nrst+1,Ncol);
z = zeros(Nrst,Ncol);

if isa(preAfun0,'function_handle')
    r0 = preAfun0(b);
else
    r0 = preAfun0*b;
end
h10 = norms(r0);
vh10(1,:) = h10;

Q(:,:,1) = bsxfun(@rdivide,r0,h10);

for j = 1:maxiter
    
    k = mod(j-1,Nrst)+1;
    
    if isa(preAfun0,'function_handle')
        r = preAfun0(Afun(Q(:,:,k)));
    else
        r = preAfun0*Afun(Q(:,:,k));
    end
    for i = 1:k
        H(i,k,:) = sum(conj(Q(:,:,i)).*r);
        r = r - bsxfun(@mtimes,reshape(H(i,k,:),1,[]),Q(:,:,i));
    end
    H(k+1,k,:) = norms(r);
    
    for ic = 1:Ncol
        if ~ converged(ic)
            Homega = H(1:k+1,1:k,ic);
            if cond(Homega) > 1/tol
                converged(ic) = 1;
                continue;
            end
            y = Homega\vh10(1:k+1,ic);
            z(1:k+1,ic) = vh10(1:k+1,ic)-Homega*y;
            relres(ic) = norms(z(:,ic))./normb(ic);
            x(:,ic) = x0(:,ic) + reshape(Q(:,ic,1:k),N,[])*y;
            if relres(ic) < tol
                converged(ic) = 1;
            end
        end
    end
    
    if maxel(relres) < tol
        flag = 0;
        return;
    end
    
    Q(:,:,k+1) = bsxfun(@rdivide,r,reshape(H(k+1,k,:),1,[]));
    
    % Restart every Nrst steps
    if k == Nrst
        x0 = x;
        r0 = b - Afun(x0(:,:));
        h10 = norms(r0);
        vh10(1,:) = h10;
        Q(:,:,1) = bsxfun(@rdivide,r0,h10);
    end
    
end

flag = 1;

end