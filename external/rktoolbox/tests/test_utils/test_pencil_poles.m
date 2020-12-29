function check = test_pencil_poles()
  tol = 1e-13;
  
  A = diag(sin(pi*(10:-2:-10)/20));
  
  [~, K, H] = rat_krylov(A, ones(11, 1), [5, 6]);  
  pls = util_pencil_poles(K, H);
  n1 = norm(sort(pls)-[5 6], inf);
  
  [~, K, H] = rat_krylov(A, ones(11, 1), [inf, inf]);  
  flag = isinf(util_pencil_poles(K, H));

  [~, K, H] = rat_krylov(A, ones(11, 1), [1+2i, 1-2i, -10], 'real');  
  pls = util_pencil_poles(K, H);
  n2 = norm(sort(real(pls))-sort(real([1-2i, 1+2i, -10])), inf);
  n3 = norm(sort(imag(pls))-sort(imag([1-2i, 1+2i, -10])), inf);

  xi = [4i, -4i, 2-2i, 2+2i, -10];
  [~, K, H] = rat_krylov(A, ones(11, 1), xi, 'real');  
  pls = util_pencil_poles(K, H);
  n4 = norm(sort(real(pls))-sort(real(xi)), inf);
  n5 = norm(sort(imag(pls))-sort(imag(xi)), inf);

  check = [n1 n2 n3 n4 n5] < tol;
  check = [check flag];
end