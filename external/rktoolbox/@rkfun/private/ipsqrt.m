function ratfun = ipsqrt(a1,b2,m,h)
%IPSQRT   Balanced Remez approximant for sqrt(x+(h*x/2)^2) 
% on [a1,b2] of degree (m,m-1), where a1 < 0 < b2.
% The approximant is generated from pre-tabulated interpolation 
% nodes for 1-interval uniform approximants of sqrt(z).

load ipsqrt

if a1 == 0,
    k1 = 0;
end
if b2 == 0,
    k1 = m;
end
if a1 < 0 && b2 > 0,
    c = log(sqrt(abs(b2/a1)))/(sqrt(2)*pi);
    k1 = round(.5*(m - sqrt(2*c^2*m-c^4)));
end
k2 = m - k1;

if k1 == 0,
    ip1 = [];
else
    ip1 = ip{k1};
end
if k2 == 0,
    ip2 = [];
else
    ip2 = ip{k2};
end

x = [ a1*ip1 ; b2*ip{k2} ];
A = diag(x);
F = diag(sqrt(x+(h*x/2).^2));
b = ones(length(x),1);
xi = inf(1, m-1);
[xi, ratfun, misfit, out] = rkfit(F, A, b, xi, ...
                                    struct('tol', 1e-15, ...
                                           'maxit', 20, ...
                                           'k', +1, ...
                                           'real', 0,...
                                           'reduction', 0));
                                       
return
