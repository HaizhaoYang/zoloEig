function check = test_expl()
  N   = 60;
  tol = 1e-12;

  A = eye(N) + 1i*gallery('tridiag', N);
  b = sum(eye(N, 5), 2);

  xi = [-1, -2, -3, -3];
  [V, K, H] = rat_krylov(A, b, xi);
  
  [KT, HT, QT, ZT] = move_poles_expl(K, H, xi*1i);
    
  pls = util_pencil_poles(KT, HT);
  
  n1 = norm(sort(imag(pls(:)))-sort(xi(:)), inf) + norm(real(pls));
  n2 = norm(A*(V*QT')*KT-(V*QT')*HT)/(4*norm(KT)+norm(HT));
  
  [K, H] = move_poles_expl(KT, HT, 2+xi*1i);
  pls = util_pencil_poles(K, H)-2;
  n3 = norm(sort(imag(pls(:)))-sort(xi(:)), inf) + norm(real(pls));

  check = [n1 n2 n3] < tol;
end
