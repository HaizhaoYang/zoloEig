function [Uout,eigvals,relerrs] = zoloeigs(A,a,b,varargin)

[nc,r,reltol,itmax,verbose,itsolacc] = zolooptions(varargin{:});

n = size(A,1);

%========================================================================
% Zolo Eigs
if verbose
    fprintf('-----------------------------------------\n');
    tottime = cputime;
    partime = cputime;
    fprintf('Setup Calculation ...');
end

if length(a) == 2 && length(b) == 2
    a1 = a(1);
    a2 = a(2);
    b1 = b(1);
    b2 = b(2);
else
    [a1,a2,b1,b2] = eiggaps(A,a,b);
end

[a,b,gammaval,l1] = MobiusFactor(a1,a2,b1,b2);

if isempty(r)
    r = zolopara(l1,reltol);
end

[ccoef1,acoef1] = zolocoef(l1,r(1),1);

if verbose
    partime = cputime-partime;
    fprintf('\b\b\b      %12.2f secs\n',partime);
    fprintf('    Relative Gap size            %.2e\n', ...
        min(a2-a1,b2-b1)/(b2-a1));
    if length(r) == 1
        fprintf('    Zolo parameters           %.2e,%2d\n',l1,r);
    else
        l2 = zoloeval(l1,l1,r(1),1);
        fprintf('    Zolo1 parameters          %.2e,%2d\n',l1,r(1));
        fprintf('    Zolo2 parameters          %.2e,%2d\n',l2,r(2));
    end
    fprintf('PreFactorization ...');
    partime = cputime;
end

prefact = cell(r(1),1);
SMF = SymbolMF(A+1i*speye(n));
for it = 1:r(1)
    ai = acoef1(it);
    ci = ccoef1(2*it-1);
    sigmai = (gammaval*b+1i*sqrt(ci)*a)/(gammaval+1i*sqrt(ci));
    wi = ai*(sigmai-a)/(2*(gammaval+1i*sqrt(ci)));
    prefact{it} = Multifrontal((A-sigmai*speye(n))/wi,SMF);
end
cc = sum(acoef1.*gammaval./(gammaval^2+ccoef1(1:2:end)));

if verbose
    partime = cputime-partime;
    fprintf('\b\b\b       %12.2f secs\n',partime);
end

R = @(X)applyZoloMat(A,speye(n),cc,prefact,X);

if length(r) > 1
    l2 = zoloeval(l1,l1,r(1),1);
    [ccoef2,acoef2] = zolocoef(l2,r(2),0);
    RR = @(X)applyZoloFunc(R,ccoef2,acoef2,X,itsolacc);
else
    RR = R;
end

if isreal(A)
    X = randn(n,nc);
else
    X = randn(n,nc)+1i*randn(n,nc);
end
relerrs = zeros(1,itmax);
trunctol = max(1e-8,abs(zoloeval(l1,l1,r)-1));
for it = 1:itmax
    
    if verbose
        fprintf('Iter %2d ...', it);
    end
    
    % Apply rational filter to V.
    [Y,niter] = RR(X);
    X = (real(Y)+X)/2;
    
    [X,Rfac,~] = svd(X,'econ');
    addOne = 1;
    numAdded = 0;
    Xbk = X;
    while Rfac(end,end)/Rfac(1,1) > trunctol && size(X,2) < size(X,1)
        if addOne == 1
            Xtmp = randn(size(X,1),1);
            numAdded = numAdded+1;
        else
            Xtmp = randn(size(X));
            numAdded = numAdded+size(X,2);
        end
        addOne = -1*addOne;
        [Y,niter] = RR(Xtmp);
        Xtmp = (real(Y)+Xtmp)/2;
        X = [X Xtmp];
        [X,Rfac,~] = svd(X,'econ');
    end
    if numAdded == 1
        X = Xbk;
        [X,Rfac,~] = svd(X,'econ');
    end
    idx = abs(diag(Rfac)) > trunctol*Rfac(1,1);
    X = X(:,idx);
    
    
    % Compute and sort Ritz pairs.
    Am = X'*A*X;
    Am = (Am+Am')/2;
    
    [W,D] = eig(Am);
    X = X*W;
    
    idx = diag(D) > a & diag(D) < b;
    Xi = X(:,idx);
    D = D(idx,idx);
    relerr = norm(A*Xi - Xi*D, 'fro')/norm(A*Xi, 'fro');
    relerrs(it) = relerr;
    
    if verbose
        fprintf('\b\b\b             %3d/%3d cols\n', size(D,1), nc);
        fprintf('    Relative error               %.2e\n', relerr);
        fprintf('    GMRES Iteration number    %11d\n',niter);
    end
    
    if relerr < reltol
        break;
    end
    
end

if verbose
    tottime = cputime-tottime;
    fprintf('Total time              %12.2f secs\n',tottime);
    fprintf('-----------------------------------------\n');
end

if nargout == 1
    Uout = diag(D);
    return;
end

Uout = Xi;
eigvals = D;
relerrs = relerrs(1:it);

end