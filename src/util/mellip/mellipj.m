function [sn,cn,dn] = mellipj(u,alpha,tol)
%mELLIPJ Jacobi elliptic functions. MATLAB's built-in code, modified for improved accuracy.
%   [SN,CN,DN] = ELLIPJ(U,M) returns the values of the Jacobi elliptic
%   functions Sn, Cn and Dn, evaluated for corresponding elements of
%   argument U and parameter M.  U and M must be arrays of the same
%   size or either can be scalar.  As currently implemented, M is
%   limited to 0 <= M <= 1.
%
%   [SN,CN,DN] = ELLIPJ(U,M,TOL) computes the elliptic functions to
%   the accuracy TOL instead of the default TOL = EPS.
%

if nargin<2
    error(message('MATLAB:ellipj:NotEnoughInputs'));
end

u=real(u);
alpha=real(alpha);

m = sin(alpha) * sin(alpha);

%classin = superiorfloat(u,m);
classin='double';
if nargin<3
    tol = eps(classin);
end


if ~isreal(u) || ~isreal(m)
    error(message('MATLAB:ellipj:ComplexInputs'))
end

if isscalar(m)
    m = m(ones(size(u)));
end
if isscalar(u)
    u = u(ones(size(m)));
end
if ~isequal(size(m),size(u))
    error(message('MATLAB:ellipj:InputSizeMismatch'));
end

mmax = numel(u);

cn=u;
sn = cn;
dn = sn;

m = m(:).';
u = u(:).';

if any(m < 0) || any(m > 1)
    error(message('MATLAB:ellipj:MOutOfRange'));
end

% pre-allocate space and augment if needed
chunk = 10;
c(1) = real(sin(alpha));
b(1) = real(cos(alpha));
a(1) = c(1);
a = [a; zeros(chunk,mmax)];
b = [b; zeros(chunk,mmax)];
c = [c; zeros(chunk,mmax)];
a(1)=1;

n = zeros(1,mmax);
i = 1;
while any(abs(c(i,:)) > tol) && i<1000
    i = i + 1;
    if i > size(a,1)
        a = [a; zeros(chunk,mmax)];
        b = [b; zeros(chunk,mmax)];
        c = [c; zeros(chunk,mmax)];
    end
    a(i,:) = 0.5 * (a(i-1,:) + b(i-1,:));
    b(i,:) = sqrt(a(i-1,:) .* b(i-1,:));
    c(i,:) = 0.5 * (a(i-1,:) - b(i-1,:));
    in = find((abs(c(i,:)) <= tol) & (abs(c(i-1,:)) > tol));
    if ~isempty(in)
        [mi,ni] = size(in);
        n(in) = repmat((i-1), mi, ni);
    end
end
%phin = zeros(i,mmax,classin);
phin(1)=u;
phin = [phin;zeros(i-1,mmax)];


phin(i,:) = (2 .^ n).*a(i,:).*u;
while i > 1
    i = i - 1;
    in = find(n >= i);
    phin(i,:) = phin(i+1,:);
    if ~isempty(in)
        % phin(i,in) = 0.5 * ...
        %     (asin(c(i+1,in).*sin(rem(phin(i+1,in),2*pi))./a(i+1,in)) ...
        %     + phin(i+1,in));
        % phin(i,in) = 0.5 * ...
        %     (asin(c(i+1,in).*sin(rem(real(phin(i+1,in)),2*pi)) ...
        %     ./a(i+1,in)) + phin(i+1,in));
        phin(i,in) = 0.5 * ...
            (asin(c(i+1,in).*sin(real(phin(i+1,in)))./a(i+1,in)) ...
            + phin(i+1,in));
    end
end
%sn(:) = sin(rem(real(phin(1,:)),2*pi));
%cn(:) = cos(rem(real(phin(1,:)),2*pi));
sn(:) = sin(real(phin(1,:)));
cn(:) = cos(real(phin(1,:)));
dn(:) = sqrt(1 - m .* (sn(:).').^2);

% special case m = 1
%m1 = find(m==1);
%sn(m1) = tanh(u(m1));
%cn(m1) = sech(u(m1));
%dn(m1) = sech(u(m1));
% special case m = 0
%dn(m==0) = 1;
end