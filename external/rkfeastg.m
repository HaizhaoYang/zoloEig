function [V,D,relerrs] = rkfeastg(A,B,lmin,lmax,m,r_feast,opt)

if nargin < 6
    opt.maxiter = 100;
    opt.reltol = 1e-10;
    opt.verbose = 0;
end

if ~isfield(opt,'maxiter')
    maxiter = 100;
else
    maxiter = opt.maxiter;
end

if ~isfield(opt,'reltol')
    reltol = 1e-10;
else
    reltol = opt.reltol;
end

if ~isfield(opt,'verbose')
    verbose = false;
else
    verbose = opt.verbose;
end

if verbose
    fprintf('-----------------------------------------\n');
    tottime = cputime;
    partime = cputime;
    fprintf('RK Construction ...');
end

N = size(A,1);

x = rkfun();                                   % independent variable x
t = 2/(lmax-lmin)*x - (lmin+lmax)/(lmax-lmin); % map [lmin,lmax] -> [-1,1]
s = rkfun('step', r_feast);
r = s(t);
r = factorization(r,A,B);

if verbose
    partime = cputime-partime;
    fprintf('\b\b\b        %12.2f secs\n',partime);
end

V = randn(N, m+5);
relerrs = zeros(1,maxiter);
for iter = 1:maxiter
    if verbose
        fprintf('Iteration %4d ...',iter);
        partime = cputime;
    end
    % Apply rational filter to V.
    V = r(A,V);
    [V,~] = qr(V,0);
    % Compute and sort Ritz pairs.
    Am = V'*A*V; Bm = V'*B*V;
    [W, D]   = eig(Am, Bm);
    [D, ind] = sort(diag(D)); W = W(:, ind);
    % B-normalize W.
    nrm = sqrt(diag(W'*Bm*W)); W = W/diag(nrm);
    % Form approximate eigenvectors.
    V = V*W;
    if verbose
        partime = cputime-partime;
        fprintf('\b\b\b         %12.2f secs\n',partime);
    end
    % Check residuals and number of eigenpairs inside
    % search iterval.
    Di = diag(D(lmin < D & D < lmax));
    Vi = V(:, lmin < D & D < lmax);
    relerrs(iter) = norm(A*Vi - B*Vi*Di, 'fro')/norm(A*Vi);
    if relerrs(iter) < reltol ...
            || ( iter > 1 && relerrs(iter) > relerrs(iter-1) )
        break;
    end
end

if verbose
    tottime = cputime-tottime;
    fprintf('Total time              %12.2f secs\n',tottime);
    fprintf('-----------------------------------------\n');
end

D = diag(Di);

if nargout == 1
    V = D(1:m);
    return;
end

V = V(:,1:m);

D = D(1:m);

relerrs = relerrs(1:iter);

end