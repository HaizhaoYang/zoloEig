function check = test_residue()
  N = 65;
  
  tol = 1e-12;
  
  A = .01*N^2*gallery('tridiag',N);
  v = eye(N,3);
  fm = @(X) sqrtm(full(A));
  F = fm(A);
  exact = F*v;
  xi = inf(1, 8);
  
  param.k = -1;
  param.maxit = 5;
  
  [~, ratfun, misfit] = rkfit(F, A, v, xi, param);
  
  m1 = misfit(end);
  
  [res,xi,absterm] = residue(ratfun);
  teval = absterm*v;
  for j = 1:length(xi)
    teval = teval + res(j)*((A - xi(j)*speye(size(A)))\v);
  end
  m2 = norm(teval - ratfun(A, v))/norm(exact);

  check = m1*1e-5 < tol &&  m2 < tol;
end