function check = test_impl()
  N   = 60;
  tol = 1e-14;

  A = gallery('tridiag', N);
  b = eye(N, 1);

  [V, K, H] = rat_krylov(A, b, [-1, -3, -6]);
  
  c = 0*H(:, 1); c(end) = 2;
  [KT, HT, QT, ZT] = move_poles_impl(K, H, c, 'real');
    
  pls = util_pencil_poles(KT, HT);
  rts = eig(H(1:end-1, :), K(1:end-1, :));
  
  n1 = norm(sort(pls(:))-sort(rts(:)), inf);
  n2 = norm(A*(V*QT')*KT-(V*QT')*HT)/(4*norm(KT)+norm(HT));
  
  [K, H] = move_poles_impl(KT, HT, QT(:, 1), 'real');
  pls = util_pencil_poles(K, H);
  ori = [-1, -3, -6];
  
  n3 = norm(sort(real(pls(:)))-sort(real(ori(:))), inf);

  check = [n1 n2 n3/200] < tol;
end
