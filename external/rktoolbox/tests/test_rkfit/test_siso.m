function check = test_siso()
  n = 30;
  N = n+1;
  A = spdiags(linspace(-1, 1, N).', 0, N, N) ;
  v = sqrt(1:N).';
  v = v/norm(v);
  tol = 1e-13;

  F = eye(N)/(A+1.1*speye(N));
  F = F*A;
  F = F/(A-1.1*speye(N));

  xi = [1+2i, 1-2i];
  
  param.real      = 1;
  param.reduction = 0;
  param.maxit     = 4;
  param.tol       = tol;
  
  [xi, ratfun, misfit] = rkfit(F, A, v, xi, param);
  
  n1 = norm(sort(xi)-sort([-1.1 1.1]), inf);

  V = rat_krylov(A, v, xi);  
  Fv = F*v;
  sc = norm(Fv);
  n2 = norm(ratfun(A, v)-Fv)/sc;
  n3 = norm(ratfun(A, v)-V*(V'*Fv))/sc;
    
  check = n1 < sqrt(tol) && n2 < tol && n3 < tol;
end
