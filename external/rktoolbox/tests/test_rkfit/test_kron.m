function check = test_kron()
  A = gallery('tridiag', 20);
  v = eye(20, 10);
  F = expm(A);

  tol = 1e-13;

  xi = inf(6, 1);
  [~, ratfun1, misfit1] = rkfit(F, A, v, xi, 'real');

  AA = kron(eye(10), A);
  FF = kron(eye(10), F);
  vv = v(:);
  [~, ratfun2, misfit2] = rkfit(FF, AA, vv, xi, 'real');

  m1 = misfit1(end);
  m2 = misfit2(end);

  r1 = norm(ratfun1(A, v) - ratfun2(A, v))/norm(ratfun1(A, v));
  r2 = norm(ratfun1(A, v)-F*v)/norm(ratfun1(A, v));

  check = [m1 m2 r1 r2] < tol;
end
