function [k,e] = mellipke(alpha,tol)
%ELLIPKE Complete elliptic integral. Modified from Matlab's built-in code
%for improved accuracy.
if nargin<1
    error(message('MATLAB:ellipke:NotEnoughInputs'));
end

a0 = 1;
i1 = 0;
mm = 1;

m = sin(alpha) * sin(alpha);

if nargin<2
    tol = eps;
end
if ~isreal(m),
    error(message('MATLAB:ellipke:ComplexInputs'))
end
if isempty(m)
    k = zeros(size(m));
    e = k;
    return;
end
if any(m(:) < 0) || any(m(:) > 1),
    error(message('MATLAB:ellipke:MOutOfRange'));
end

b0 = cos(alpha);
s0 = m;
while mm > tol
    a1 = (a0+b0)/2;
    b1 = sqrt(a0.*b0);
    c1 = (a0-b0)/2;
    i1 = i1 + 1;
    w1 = 2^i1*c1.^2;
    mm = max(w1(:));
    s0 = s0 + w1;
    a0 = a1;
    b0 = b1;
end
k = pi./(2*a1);
e = k.*(1-s0/2);



end