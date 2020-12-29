function check = test_type_21()
  tol = 1e-14;
  
  % Case abs(a2) > abs(a3).  
  r   = rkfun.nodes2rkfun([13+7i, 13-7i], -5);
  rts = roots(r);
  n1  = 1./feval(r, 1e30);  
  n2  = abs(sum(rts)-26)/50 + abs(prod(rts)-218)/400;    
  n3  = abs(poles(r) + 5);
  n4  = abs(feval(r, 1) - 193/6)/20;

  % Case abs(a2) < abs(a3).
  r   = rkfun.nodes2rkfun([.1+.4i, .1-.4i], [inf, .2]);
  rts = roots(r);
  n5  = 1./feval(r, 1e30);
  n6  = abs(sum(rts)-.2) + abs(prod(rts)-.17);
  n7  = abs(poles(r) - .2);
  n8  = abs(feval(r, -1) + 1.141666666666667);

  check = [n1 n2 n3 n4 n5 n6 n7 n8] < tol;
end
