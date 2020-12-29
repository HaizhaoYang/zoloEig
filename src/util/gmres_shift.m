function [x,flag,relres,j] = gmres_shift ...
    (Afun0,b,omega,restart,tol,maxiter)

% solves the systems (A-omega I)x = rhs
% ||(A-omega I)x - rhs||\(||rhs||) < e for all omega

% INPUTS
% A:     a function handle to apply a NxN matrix
% rhs:   a Nx1 vector
% omega: a vector containing the values of omega

% OUTPUTS
% x:       a NxNomega vector representing the solution, col(k) is for omega(k)
% iter:    The number of iterations needed, will be less than maxiter
% err:     ||(A-omega I)x - rhs||\(||n||), for each shift. err(k) is for omega(k)



if nargin < 6 || isempty(maxiter)
    maxiter = 100;
end

if nargin < 5 || isempty(tol)
    tol = 1e-6;
end

if nargin < 4 || isempty(restart)
    restart = 30;
end

if nargin < 3 || isempty(omega)
    omega = 0;
end

norms = @(X)sqrt(sum(abs(X).^2,1));

[N,Ncol] = size(b);
Nomega = length(omega);
Nrst = restart;
normb = norms(b);

if ~isa(Afun0,'function_handle')
    A = Afun0;
    Afun = @(x)A*reshape(x,N,[]);
else
    Afun = @(x)Afun0(reshape(x,N,[]));
end

% index for reference of omega
[~,soidx] = min(abs(omega));

converged = zeros(Ncol,Nomega);
x = zeros(N,Ncol,Nomega);
relres = zeros(Ncol,Nomega);

Q = zeros(N,Ncol,Nrst);
H = zeros(Nrst+1,Nrst,Ncol);
vh10 = zeros(Nrst+1,Ncol);
beta = ones(Ncol,Nomega);
z = zeros(Nrst,Ncol);

x0 = zeros([size(b),Nomega]);
r0 = b;
h10 = norms(r0);
vh10(1,:) = h10;
beta0 = beta;


Q(:,:,1) = bsxfun(@rdivide,r0,h10);

for j = 1:maxiter
    
    k = mod(j-1,Nrst)+1;
    
    r = Afun(Q(:,:,k)) - omega(soidx)*Q(:,:,k);
    for i = 1:k
        H(i,k,:) = sum(conj(Q(:,:,i)).*r);
        r = r - bsxfun(@mtimes,reshape(H(i,k,:),1,[]),Q(:,:,i));
    end
    H(k+1,k,:) = norms(r);
    
    for ic = 1:Ncol
        if ~ converged(ic,soidx)
            Homega = H(1:k+1,1:k,ic);
            y = Homega\vh10(1:k+1,ic);
            z(1:k+1,ic) = vh10(1:k+1,ic)-Homega*y;
            relres(ic,soidx) = norms(z(:,ic))./normb(ic);
            x(:,ic,soidx) = x0(:,ic,soidx) + reshape(Q(:,ic,1:k),N,[])*y;
            if relres(ic,soidx) < tol
                converged(ic,soidx) = 1;
            end
        end
    end
    
    for io = 1:Nomega
        if ~ all(converged(:,io)) && io ~= soidx
            Hshift = zeros(k+1,k);
            Hshift((0:k-1)*(k+1)+(1:k)) = omega(io)-omega(soidx);
            for ic = 1:Ncol
                if ~ converged(ic,io)
                    Homega = [H(1:k+1,1:k,ic)-Hshift,z(1:k+1,ic)];
                    if cond(Homega) > 1/tol
                        converged(ic,io) = 1;
                        continue;
                    end
                    y = Homega\(vh10(1:k+1,ic)*beta0(ic,io));
                    beta(ic,io) = y(end);
                    relres(ic,io) = ...
                        norms(vh10(1:k+1,ic) - Homega(:,1:k)*y(1:k)) ...
                        ./normb(ic);
                    x(:,ic,io) = x0(:,ic,io) ...
                        + reshape(Q(:,ic,1:k),N,[])*y(1:k);
                    
                    if relres(ic,io) < tol
                        converged(ic,io) = 1;
                    end
                end
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
        r0 = b - Afun(x0(:,:,soidx))...
            + omega(soidx)*reshape(x0(:,:,soidx),N,[]);
        h10 = norms(r0);
        vh10(1,:) = h10;
        beta0 = beta;
        Q(:,:,1) = bsxfun(@rdivide,r0,h10);
    end
    
end

flag = 1;

end