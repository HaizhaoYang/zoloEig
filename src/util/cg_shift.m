function [x,flag,relres,it] = cg_shift(Afun,rhs,omega,tol,maxiter)

% solves the systems (A-omega I)x = rhs
% ||(A-omega I)x - rhs||\(||b||) < e for all omega

% Based on code given in Algorithm 9.3.1 in Golub and Van Loan

% INPUTS
% A:     a function handle to apply a NxN matrix, must be positive definite
% rhs:   a Nx1 vector
% omega: a vector containing the values of omega
% e:     a scaler

% OUTPUTS
% x:       a NxNomega vector representing the solution, col(k) is for omega(k)
% iter:    The number of iterations needed, will be less than maxiter
% err:     ||(A-omega I)x - rhs||\(||n||), for each shift. err(k) is for omega(k)

if nargin < 5 || isempty(maxiter)
    maxiter = 100;
end

if nargin < 4 || isempty(tol)
    tol = 1e-6;
end

if nargin < 3 || isempty(omega)
    omega = 0;
end

if ~isa(Afun,'function_handle')
    A = Afun;
    Afun = @(x)A*x;
end

%norms = @(X)sqrt(sum(abs(X).^2,1));
norms = @(X)sqrt(sum(abs(X).^2,1));

[N,Ncol] = size(rhs);
Nomega = length(omega);

converged = zeros(1,Ncol,Nomega);
u = zeros(1,Ncol,Nomega);
d = zeros(1,Ncol,Nomega);
c = zeros(N,Ncol,Nomega);
x = zeros(N,Ncol,Nomega);
rho = zeros(1,Ncol,Nomega);
err = zeros(1,Ncol,Nomega);

Nb = norms(rhs);
r0 = rhs;
beta0 = norms(r0);
q0 = zeros(size(rhs));

for it = 1:maxiter
    
    % q = r0/beta0;
    q = bsxfun(@rdivide,r0,beta0);
    z = Afun(q);
    
    % alpha = q'*z;
    alpha = dot(q,z);
    
    % r = z - alpha*q-beta0*q0;
    r = z - bsxfun(@mtimes,alpha,q)-bsxfun(@mtimes,beta0,q0);
    beta = norms(r);
    
    for m = 1:Nomega
        if ~ converged(m)
            if it == 1
                d(1,:,m) = alpha-omega(m);
                c(:,:,m) = q;
                rho(1,:,m) = beta0./d(1,:,m);
                % x(:,m) = rho(m)*q;
                x(:,:,m) = bsxfun(@mtimes,rho(1,:,m),q);
            else
                u(1,:,m) = beta0./d0(1,:,m);
                d(1,:,m) = alpha-omega(m) - beta0.*u(1,:,m);
                c(:,:,m) = q-bsxfun(@mtimes,u(1,:,m),c(:,:,m));
                rho(1,:,m) = -u(1,:,m).*d0(1,:,m).*rho(1,:,m)./d(1,:,m);
                x(:,:,m) = x(:,:,m) + bsxfun(@mtimes,rho(1,:,m),c(:,:,m));
            end
            
            err(1,:,m) = abs(beta.*rho(1,:,m))./Nb;
            
            if max(err(1,:,m)) < tol
                converged(m) = 1;
            end
        end
    end
    
    if maxel(err) < tol
        flag = 0;
        relres = err;
        return;
    end
    
    r0 = r;
    q0 = q;
    beta0 = beta;
    d0 = d;
end

flag = 1;
relres = err;

end