function check = test_hh2th()
  N   = 50;
  tol = 1e-12;

  A = gallery('moler', N);
  b = eye(N, 1);
  nrm = 1e3;
  
  xi = [-1, -2, -3, -3i, 3i, 2+2i, 2-2i];  
  [V, K, H] = rat_krylov(A, b, xi, 'real');
    
  [Ke, He, Qe, Ze] = move_poles_expl(K, H, inf*xi);
  [KT, HT, QT, ZT] = util_hh2th(K, H);

  n1 = norm(util_pencil_poles(He, Ke));
  n2 = norm(util_pencil_poles(HT, KT));

  VT = V*QT'; Ve = V*Qe';
  n3 = norm(A*VT*KT-VT*HT)/(nrm*norm(KT)+norm(HT));  
  n3 = n3 + norm(A*Ve*Ke-Ve*He)/(nrm*norm(Ke)+norm(He));

  
  xi = [-1, -2, -3, inf];  
  [V, K, H] = rat_krylov(A, b, xi);
    
  [Ke, He, Qe, Ze] = move_poles_expl(K, H, inf*xi);
  [KT, HT, QT, ZT] = util_hh2th(K, H);

  n4 = norm(util_pencil_poles(He, Ke));
  n5 = norm(util_pencil_poles(HT, KT));

  VT = V*QT'; Ve = V*Qe';
  n6 = norm(A*VT*KT-VT*HT)/(nrm*norm(KT)+norm(HT));  
  n6 = n6 + norm(A*Ve*Ke-Ve*He)/(nrm*norm(Ke)+norm(He));

  check = [n1 n2 n3 n4 n5 n6] < tol;
end
