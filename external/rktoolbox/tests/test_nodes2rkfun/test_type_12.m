function check = test_type_12()
  tol = 1e-14;
  
  % Case abs(a2) > abs(a3).       
  r  = rkfun.nodes2rkfun(2, [7+9i, 7-9i]);
  n1 = abs(roots(r)-2) + abs(feval(r, 1e50));
  n2 = norm(real(poles(r)) - [7; 7], inf);
  n3 = norm(sort(imag(poles(r))) - sort([9; -9]), inf);
  
  % Case zeta numerically equals real(xi).
  r  = rkfun.nodes2rkfun(7+2e-16, [7+9i, 7-9i]);
  n4 = abs(roots(r)-7-2e-16) + abs(feval(r, 1e50));
  n5 = norm(real(poles(r)) - [7; 7], inf);
  n6 = norm(sort(imag(poles(r))) - sort([9; -9]), inf);

  % Case abs(a3) > abs(a2).      
  r  = rkfun.nodes2rkfun(6, [2+.1i, 2-.1i]);   
  n7 = abs(roots(r)-6) + abs(feval(r, 1e50));
  n8 = norm(real(poles(r)) - [2; 2], inf);
  n9 = norm(sort(imag(poles(r))) - sort([.1; -.1]), inf);

  check = [n1 n2 n3 n4 n5 n6 n7 n8 n9] < tol;
end
