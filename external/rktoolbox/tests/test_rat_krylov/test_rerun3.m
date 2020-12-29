function check = test_rerun3()
  N   = 64;
  tol = 1e-14;

  A = gallery('poisson', sqrt(N));
  b = ones(N, 1); 
  b = b/norm(b);
  
  xi = [inf repmat([1 3.5 4.15 5], 1, 2)];

  [V, K, H] = rat_krylov(A, b, xi);

  scale = 8*norm(K)+norm(H);
  n1 = norm(A*V*K-V*H)/scale;
  n2 = norm(V'*V-eye(size(V, 2)))/length(xi);
  
  [Vr, Kr, Hr] = rat_krylov(A, b, K, H);
  
  n3 = norm(A*Vr*Kr-Vr*Hr)/scale;
  n4 = norm(V*(V'*Vr)-Vr)/scale;

  B = gallery('toeppen', N);
  scale = 8*norm(K)+21*norm(H);
  [Vr, Kr, Hr] = rat_krylov(B*A, B, b, K, H);
  scale = scale*norm(Vr);
  n5 = norm(B*A*Vr*K-B*Vr*H)/scale;
  
  check = [n1 n2 n3 n4 n5] < tol;
end
