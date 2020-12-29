function obj = factorization( obj, A, B )

if nargin < 3
    B = speye(size(A));
end

xi = util_pencil_poles(obj.K, obj.H);
[nu,mu] = poles_to_moebius(xi);

obj.prefact.nu = nu;
obj.prefact.mu = mu;

SMF = SymbolMF(A+1i*B);
obj.prefact.MF = cell(1,length(nu));

obj.prefact.A = A;
obj.prefact.B = B;

obj.prefact.MFB = Multifrontal(B,SMF);

for it = 1:length(nu)
    obj.prefact.MF{it} = Multifrontal((nu(it)*A-mu(it)*B),SMF);
end

    function [nu,mu] = poles_to_moebius(xi)
        % POLES_TO_MOEBIUS Moebius transformation with poles xi.
        %
        % Finite xi is replaced with (nu, mu) := (1, xi) and
        % (rho, eta) := (0, 1),  and xi = inf is replaced by
        % (nu, mu) := (0, 1) and (rho, eta) := (1, 0).
        
        nu  = ones(1, length(xi));
        mu  = xi;
        
        %nu(isinf(xi))  = 0; % does not work with vpa
        nu(abs(xi)==inf) = 0;
        %mu(isinf(xi))  = 1; % does not work with vpa
        mu(abs(xi)==inf) = 1;
        
    end

end

