function check = test_type_33()
  tol = 1e-12;

  r1 = rkfun.nodes2rkfun(1,  7);
  r2 = rkfun.nodes2rkfun(2,  8);
  r3 = rkfun.nodes2rkfun(-1, 0);
  
  r  = rkfun.nodes2rkfun([2, 1, -1], [0, 8, 7]);

  ev  = [1.7, pi, 1+3i, 1-3i, 2.8, -2, 3, 5+15i, 6];
  rev = r(ev);
  
  n1 = norm((r1(ev).*r2(ev).*r3(ev)-rev)./(rev), inf);
  n2 = ~(isreal(r) && isreal(r1) && isreal(r2));

  check = [n1 n2] < tol;
end
