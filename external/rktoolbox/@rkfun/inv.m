function obj = inv(a)
%INV    Inverse RKFUN corresponding to a Moebius transform.

if length(a.coeffs) > 2
    error(['INV: Function inversion is currently only possible for' ...
        ' rkfuns of type at most [1,1].']);
end

% Get Moebius transform parameters (rho*z - eta)/(nu*z - mu).
if length(a.coeffs) == 1
    rho = 0; eta = a.coeffs(1); mu = 1; nu = 0;
elseif length(a.coeffs) == 2
    rho = a.coeffs(1)*a.K(2, 1) - a.coeffs(2)*a.K(1, 1);
    eta = a.coeffs(1)*a.H(2, 1) - a.coeffs(2)*a.H(1, 1);
    mu  = a.H(2, 1);
    nu  = a.K(2, 1);
end

obj = gallery('moebius', [-mu, eta, -nu, rho]);
end