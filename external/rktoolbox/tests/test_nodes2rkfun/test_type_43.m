function check = test_type_43()
  tol = 1e-12;

  r1 = rkfun.nodes2rkfun([5+6i, 5-6i], 7);
  r2 = rkfun.nodes2rkfun([1, -inf, 2], [3i, -3i]);
  
  r  = rkfun.nodes2rkfun([5+6i, 1, 5-6i, inf, 2], [NaN, 3i, 7, -3i]);
  
  ev  = [1.7, pi, 1+3i, 1-3i, 2.8, -2, 3, 5+15i, 0];
  rev = r(ev);
  
  n1 = 1./feval(r, 1e20);
  n2 = norm((r1(ev).*r2(ev)-rev)./(rev), inf);
  n3 = ~(isreal(r) && isreal(r1) && isreal(r2));
  
  check = [n1 n2 n3] < tol;
end
