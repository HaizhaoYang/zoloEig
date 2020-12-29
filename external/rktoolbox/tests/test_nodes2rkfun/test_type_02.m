function check = test_type_02()
  tol = 1e-14;
  
  r  = rkfun.nodes2rkfun([], [5+2i, 5-2i]);
  n1 = feval(r, 1e100);
  n2 = norm(real(poles(r)) - [5; 5], inf);
  n3 = norm(sort(imag(poles(r))) - sort([2; -2]), inf);
  n4 = abs(feval(r, 1) - .05);

  check = [n1 n2 n3 n4] < tol;
end
