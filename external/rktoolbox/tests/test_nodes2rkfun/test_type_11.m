function check = test_type_11()
  tol = 1e-14;
  
  r  = rkfun.nodes2rkfun(2, 3);  
  n1 = abs(roots(r)-2);
  n2 = abs(poles(r)-3);  
  n3 = abs(feval(r, 5) - 1.5);
  
  r  = rkfun.nodes2rkfun([], 3);  
  n4 = feval(r, 1e100);
  n5 = abs(poles(r)-3);  
  n6 = abs(feval(r, 5) - 0.5);
     
  r  = rkfun.nodes2rkfun(2, inf);
  n7 = abs(roots(r)-2);
  n8 = 1/feval(r, 1/tol^2);  
  n9 = abs(feval(r, 5) - 3);
  
  check = [n1 n2 n3 n4 n5 n6 n7 n8 n9] < tol;
end
