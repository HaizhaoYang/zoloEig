function check = test_siso_red()
  n = 30;
  N = n+1;
  z = sin(pi*(n:-2:-n)/(2*n)).';
  A = spdiags(z, 0, N, N) ;
  v = sqrt(1:N).';
  tol = 1e-13;

  AB.isreal   = 1;
  AB.multiply = @(rho, eta, x) rho*(z.*x)-eta*x;     
  AB.solve    = @(nu,  mu,  x) x./(nu*z-mu);
  
  F = eye(N)/(A+1.1*speye(N));
  F = F*A;
  F = F/(A-1.1*speye(N));

  xi = [-sin(pi/N), sin(pi/N)];
  [xi, ratfun, misfit, out] = rkfit(F, AB, v, xi, 5, tol, 'real');
  
  n1 = norm(sort(xi)-sort([-1.1 1.1]), inf);
  f1 = out.k_reduced == -1;
  
  V = rat_krylov(AB, v, xi);  
  Fv = F*v;
  n2 = norm(ratfun(AB, v)-Fv)/norm(Fv);
  
  check = n1 < sqrt(tol) && n2 < tol && f1;
end
