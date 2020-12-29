function check = test_pole()
  A = [3 1 0; 1 2-3i 0; 0 0 1];
  F = expm(A);
  v = [-1i; 1; 1];
  v = v/norm(F*v);

  tol = 1e-14;

  xi = inf;
  [xi, ratfun, misfit] = rkfit(F, A, v, xi);

  check(1) = (misfit/norm(F*v)) < tol;
  check(2) = abs(ratfun(xi+10*eps)) > 1/tol;
end
