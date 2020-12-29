function x = util_continuation_fom(AB, nu, mu, x, param)  
% Approximates (nu*A-mu*B)\x using param.continuation_m iterations of FOM.
% AB is a pencil structure with fields 'isreal', 'multiply', and 'solve'. See
% rat_krylov for more details.
  
  Pencil.isreal   = AB.isreal && isreal(nu) && isreal(mu);
  Pencil.multiply = @(rho, eta, x) rho*AB.multiply(nu, mu, x) - eta*x;  
  Pencil.solve    = @(nu,  mu,  x) -mu*x;
  
  if ~Pencil.isreal, param.real = 0; end
  
  param.waitbar = false;
  param.continuation = 'last';
  param.continuation_root = NaN;
  param.p = 1;
  
  [V, K, H] = rat_krylov(Pencil, x, inf(1, param.continuation_m), param);  
  
  H = H/K(1:end-1, :);
      
  x = V(:, 1:end-1)*(H(1:end-1, :)\(param.inner_product(x, V(:, 1:end-1))));  
end
  
