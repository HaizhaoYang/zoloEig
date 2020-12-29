function [V,D,relerrs] = rkfeast(A,lmin,lmax,m,r_feast,opt)

if nargin < 6
    opt.itmax = 100;
    opt.reltol = 1e-10;
    opt.verbose = 0;
end

if ~isfield(opt,'itmax')
    itmax = 100;
else
    itmax = opt.itmax;
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
r = factorization(r,A);

if verbose
    partime = cputime-partime;
    fprintf('\b\b\b        %12.2f secs\n',partime);
end

V = randn(N, m);
relerrs = zeros(1,itmax);
for iter = 1:itmax
    if verbose
        fprintf('Iteration %4d ...',iter);
        partime = cputime;
    end
    % Apply rational filter to V.
    V = r(A,V);
    [V,~] = qr(V,0);
    % Compute and sort Ritz pairs.
    Am = V'*A*V;
    [W, D]   = eig(Am);
    [D, ind] = sort(diag(D));
    W = W(:, ind);
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
    relerrs(iter) = norm(A*Vi - Vi*Di, 'fro')/norm(A*Vi);
    if relerrs(iter) < reltol
        break;
    end
end

if verbose
    tottime = cputime-tottime;
    fprintf('Total time              %12.2f secs\n',tottime);
    fprintf('-----------------------------------------\n');
end

if length(diag(Di)) == m
    D = diag(Di);
end

if nargout == 1
    V = D(1:m);
    return;
end

V = V(:,1:m);

D = D(1:m);

relerrs = relerrs(1:iter);

end