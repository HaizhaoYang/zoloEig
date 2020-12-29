function [Uout,eigvals,relerrs] = zoloeigsdense(Afun,GenAinvfun,n,a,b,varargin)

[nc,r,reltol,itmax,verbose,itsolacc] = zolooptions(varargin{:});

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
    [a1,a2,b1,b2] = eiggaps(Afun,n,a,b);
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

Ainvzolo = cell(r(1),1);
for it = 1:r(1)
    ai = acoef1(it);
    ci = ccoef1(2*it-1);
    sigmai = (gammaval*b+1i*sqrt(ci)*a)/(gammaval+1i*sqrt(ci));
    wi = ai*(sigmai-a)/(2*(gammaval+1i*sqrt(ci)));
    Ainvfun = GenAinvfun(sigmai);
    Ainvzolo{it} = @(x) wi*Ainvfun(x);
end
cc = sum(acoef1.*gammaval./(gammaval^2+ccoef1(1:2:end)));

if verbose
    partime = cputime-partime;
    fprintf('\b\b\b       %12.2f secs\n',partime);
end

R = @(X)applyZoloDenseMat(cc,Ainvzolo,X);

if length(r) > 1
    l2 = zoloeval(l1,l1,r(1),1);
    [ccoef2,acoef2] = zolocoef(l2,r(2),0);
    RR = @(X)applyZoloFunc(R,ccoef2,acoef2,X,itsolacc);
else
    RR = R;
end

X = randn(n,nc);
relerrs = zeros(1,itmax);
for it = 1:itmax
    
    if verbose
        fprintf('Iter %2d ...', it);
    end
    
    % Apply rational filter to V.
    [Y,niter] = RR(X);
    X = (Y+X)/2;
    [X,~] = qr(X,0);
    % Compute and sort Ritz pairs.
    Am = X'*Afun(X);
    Am = (Am+Am')/2;
    
    [W,D] = eig(Am);
    X = X*W;
    
    idx = diag(D) > a & diag(D) < b;
    Xi = X(:,idx);
    D = D(idx,idx);
    relerr = norm(Afun(Xi) - Xi*D, 'fro')/norm(Afun(Xi), 'fro');
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