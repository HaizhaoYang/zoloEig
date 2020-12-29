function check = test_type_22_cplx_poles()
  tol = 5e-13;
  
  % Case abs(a2) < abs(a3). 
  r  = rkfun.nodes2rkfun([2-7i, 2+7i], [3+5i, 3-5i]);
  rts = roots(r);
  pls = poles(r);  
  n1 = abs(sum(rts)-4) + abs(prod(rts)-53);
  n2 = abs(sum(pls)-6) + abs(prod(pls)-34);
  n3 = abs(feval(r, 1) - 50/29);  
  
  % Case abs(a2) > abs(a3).  
  r  = rkfun.nodes2rkfun([2-4i, 2+4i], [3+5i, 3-5i]);
  rts = roots(r);
  pls = poles(r);  
  n4 = abs(sum(rts)-4) + abs(prod(rts)-20);
  n5 = abs(sum(pls)-6) + abs(prod(pls)-34);
  n6 = abs(feval(r, pi) - 6.915747505794193e-01);
  
  % Real roots.  
  r  = rkfun.nodes2rkfun([3, -1.5], [-1i, 1i]);
  rts = roots(r);
  pls = poles(r); 
  n7 = abs(sum(rts)-1.5) + abs(prod(rts)+4.5);
  n8 = abs(sum(pls)) + abs(prod(pls)-1);
  n9 = abs(feval(r, 1) + 2.5);
  
  check = [n1 n2 n3 n4 n5 n6 n7 n8 n9] < tol;
end
