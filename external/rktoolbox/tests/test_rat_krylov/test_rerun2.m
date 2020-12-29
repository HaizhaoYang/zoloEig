function check = test_rerun2()
  N   = 100;
  tol = 1e-14;

  A = gallery('tridiag', N);
  b = ones(N, 1); 
  b = b/norm(b);
  
  xi = repmat([1.3+2i, 1.3-2i, 7i, -7i, 0.4, 0.8], 1, 2);
  param.real = 1;
  [V, K, H] = rat_krylov(A, b, xi, param);
  scale = norm(K)+norm(H); scale = 5*scale;
  n1 = norm(A*V*K-V*H)/scale;
  n2 = norm(V'*V-eye(size(V, 2)))/length(xi);
  
  [Vr, Kr, Hr] = rat_krylov(A, b, K, H, param);
  n3 = norm(A*Vr*Kr-Vr*Hr)/scale;
  n4 = norm(V*(V'*Vr)-Vr)/scale;
  
  B = gallery('toeppen', N);
  scale = scale*21;
  [Vr, Kr, Hr] = rat_krylov(B*A, B, b, K, H, param);
  n5 = norm(B*A*Vr*K-B*Vr*H)/scale;
  n6 = norm(V*(V'*Vr)-Vr)/scale;

  check = [n1 n2 n3 n4 n5 n6] < tol;
end
