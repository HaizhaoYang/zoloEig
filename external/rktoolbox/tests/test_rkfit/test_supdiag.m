function check = test_supdiag()
  N = 22;
  A = gallery('tridiag', N);
  I = eye(N);
  v = eye(N, 1);

  tol = 1e-13;

  F = ((A-12*I)/(A+2*I))*(A-5*I);
  F = F*(A/(A+3*I));
  F = F*((A+13*I)/(A+4*I));
  F = F*((A- 8*I)/(A+6*I));

  param.k = 1;
  param.maxiter = 10;
  param.real = 1;
  xi = inf(4, 1);
  [xi, ratfun, misfit] = rkfit(F, A, v, xi, param);

  m1 = misfit(end);
  m2 = norm(sort(xi)-sort([-2 -3 -4 -6]))*1e-6;
  r1 = 1./ratfun(1e150);
  r2 = 1./ratfun(-1e100);

  check = abs([m1 m2 r1 r2]) < tol;
end
