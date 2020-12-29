function [ Y, Ys ] = chebJacSlice( H, X, m, a, b, a0, b0, isJac, Ys )
% 1) Interplate a step function defined on [a0,b0] with a support on [a,b]
% in [a0,b0] using Chebyshev series with highest order m. 2) Filter X with
% the filter constructed above to get Y. 3) Y is the approximate eigenspace
% corresponding to those eigenvalues of H located in [a,b].

if isempty(Ys)
    isNew = 1;
else
    isNew = 0;
end

alpha = pi/(m+2);

% scale
a = (a-(b0+a0)/2)/((b0-a0)/2);
b = (b-(b0+a0)/2)/((b0-a0)/2);
e = (b0-a0)/2;
c = (b0+a0)/2;

[sz1,sz2] = size(X);
Y = sparse(sz1,sz2);

if isa(H,'function_handle')
    Hfunc = H;
else
    Hfunc = @(x) H*x;
end

if isNew
    Ys = cell(m+1);
    Ys{1} = X;
    Ys{2} = ( Hfunc(X) - c*X ) / (c-a0);
    sigma = e / (c-a0);
end

for cnt = 0 : m
    
    if cnt == 0
        gamma = (acos(a)-acos(b))/pi;
    else
        gamma = 2*(( sin(cnt*acos(a)) - sin(cnt*acos(b)) )/cnt)/pi;
    end
    
    if isJac
        % use Chebyshev-Jackson approximation to smooth out the
        % approximation by Chebshev series
        g = ( (1 - cnt / (m + 2)) * sin(alpha) * cos(cnt*alpha) ...
            + cos(alpha) * sin(alpha) / (m+2) ) / sin(alpha);
    else
        % use original Chebyshev series
        g = 1;
    end
    
    if isNew && cnt>1
        sigma1 = 1/(2/e*(c-a0)-sigma);
        Ys{cnt+1} = (Hfunc(Ys{cnt})-c*Ys{cnt})*(2*sigma1/e) ...
            - (sigma*sigma1)*Ys{cnt-1};
        sigma = sigma1;
    end
    
    Y = Y + (gamma*g) * Ys{cnt+1};
    
end

end