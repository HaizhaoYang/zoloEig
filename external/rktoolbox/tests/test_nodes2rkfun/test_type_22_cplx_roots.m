function check = test_type_22_cplx_roots()
  tol = 5e-13;
  
  % Case abs(a2) > abs(a3).  
  r  = rkfun.nodes2rkfun([2-7i, 2+7i], [3, 3]);
  rts = roots(r);
  pls = poles(r);  
  n1 = abs(sum(rts)-4) + abs(prod(rts)-53);
  n2 = abs(sum(pls)-6) + abs(prod(pls)-9);
  n3 = abs(feval(r, 2) - 49)/50;
  
  % Case abs(a2) < abs(a3).  
  r  = rkfun.nodes2rkfun([1-4i, 1+4i], [160, -.1]);
  rts = roots(r);
  pls = poles(r);  
  n4 = abs(sum(rts)-2)/20 + abs(prod(rts)-17)/20;
  n5 = abs(sum(pls)-159.9)/200 + abs(prod(pls)+16)/20;
  n6 = abs(feval(r, -2) - 8.122157244964262e-02);

  check = [n1 n2 n3 n4 n5 n6] < tol;
end
