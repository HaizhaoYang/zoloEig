function check = test_subdiag()
  N = 25;
  A = gallery('tridiag', N);
  I = eye(N);
  v = eye(N, 1);

  tol = 1e-14;

  F = (A+2*I)\(A/(A+3*I));
  F = F*((A+13*I)\(A+4*I));
  F = F*((A- 8*I)\(A+6*I));

  param.k = -1;
  param.maxiter = 3;

  xi = inf(4, 1);
  [xi, ratfun, misfit] = rkfit(F, A, v, xi, param);

  m1 = misfit(end);
  m2 = norm(sort(xi)-sort([-2 -3 -13 8]))*1e-5;
  r1 = ratfun(1e100);
  r2 = ratfun(-1e150);

  check = abs([m1 m2 r1 r2]) < tol;
end
