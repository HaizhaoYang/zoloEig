function check = test_siso_weight()
  N = 30;
  A = spdiags(linspace(-2, 2, N).', 0, N, N);
  v = sqrt(1:N).';
  v = v/norm(v);
  tol = 1e-13;

  F = eye(N)/(A+2.1*speye(N));
  F = F*A;
  F = F/(A-2.05*speye(N));

  D = 1./(abs(diag(F)).^.5);
  D = spdiags(D, 0, N, N);
  
  xi = [1+2i, 1-2i];
  
  param.real      = 1;
  param.reduction = 0;
  param.maxit     = 4;
  param.tol       = tol;
  param.D         = D;
  
  [xi, ratfun, misfit] = rkfit(@(X) F*X, A, v, xi, param);

  n1 = norm(sort(xi)-sort([2.05 -2.1]), inf);
  
  [V, K, H] = rat_krylov(A, v, xi);    
  Fv = F*v;

  n2 = norm(D*(ratfun(A, v)-Fv))/norm(D*Fv);
  n3 = norm(ratfun(A, v)-V*(V'*Fv))/norm(Fv);

  %[n1 n2 n3]  
  check = n1 < sqrt(tol) && n2 < tol && n3 < tol;
end
