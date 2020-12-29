function check = test_recompute_poly()
  N   = 81;
  tol = 1e-13;

  A = gallery('tridiag', N, 1, 100, 0.45);
  B = speye(N);
  b = ones(N, 1); 

  b = b/norm(b); 
  
  xi = [5.4472e-01+1.9874e+00i, 6.4731e-01+4.3735e-01i, 5.4389e-01+ ...
        2.1160e-01i, 7.2105e-01+2.1939e-01i, 5.2250e-01+1.2718e-01i];  
  xi = xi - 110;
  xi = [xi conj(xi)];
  xi = util_cplxpair(xi);
  xi = [xi -50 inf];
  param.real = 1;
  
  [V, K, H] = rat_krylov(A, B, b, xi, param);
  scale = 100*norm(K)+norm(H);
  n1 = norm(A*V*K-V*H)/scale;
  n2 = norm(V'*V-eye(size(V, 2)))/length(xi);
  
  [KT, HT, QT, ZT] = util_hh2th(K, H);    
  [Vp, Kp, Hp] = rat_krylov(A, V*QT(1, :)', inf*xi);

  n3 = norm(V*(V'*Vp)-Vp)/(scale*length(xi)*10);

  check = [n1 n2 n3] < tol;
end
