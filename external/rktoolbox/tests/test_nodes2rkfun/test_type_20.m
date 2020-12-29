function check = test_type_20()
  tol = 1e-14;
  
  % Case abs(a2) > abs(a3).  
  r   = rkfun.nodes2rkfun([13+7i, 13-7i], [inf]);
  rts = roots(r);
  n1  = 1./feval(r, 1e30);  
  n2  = abs(sum(rts)-26)/26 + abs(prod(rts)-218)/218;    
  n3  = abs(feval(r, 1) - 193)/193;

  % Case abs(a2) < abs(a3).
  r   = rkfun.nodes2rkfun([.1+4i, .1-4i], [inf inf]);
  rts = roots(r);
  n4  = 1./feval(r, 1e30);
  n5  = abs(sum(rts)-.2) + abs(prod(rts)-16.01)/16;
  n6  = abs(feval(r, 1) - 16.81)/16;
  
  check = [n1 n2 n3 n4 n5 n6] < tol;
end
