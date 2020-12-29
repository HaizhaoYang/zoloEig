function check = test_unsym_III
 
  N = 50;
  A = [];
  x = linspace(1, 5, N/2);
  for j = 1:N/2
    A = blkdiag(A, [0 x(j); -x(j) 0]);
  end

  B = gallery('kms', N);  
  
  nrmA = 5;
  nrmB = 3;
  tol = 1e-14;
  
  b = B\ones(N, 1);
  
  xi = [-20 -10 10 20 1+2i 1-2i];
  m = length(xi);
  
  param.real = 1;
  param.inner_product = @(x, y) y'*(B*x);
  
  [Vo, Ko, Ho] = rat_krylov(A, B, b, xi, param);
  n1 = norm(A*Vo*Ko - B*Vo*Ho);
  n1 = n1/(norm(Vo)*(nrmA*norm(Ko) + nrmB*norm(Ho)));
  n2 = norm(param.inner_product(Vo, Vo) - eye(m+1));

  [V,  K,  H]  = rat_krylov(A, B, b, xi(1));
  for j = 2:m
    [V, K, H] = rat_krylov(A, B, V, K, H, xi(j));
  end
  n3 = norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n4 = norm(V'*V - eye(m+1)) + norm(V*(V'*Vo) - Vo);
  
  [V, K, H] = rat_krylov(A, B, b, xi(1:end-2));
  [V, K, H] = rat_krylov(A, B, V, K, H, xi(end-1:end));
  n3 = n3 + norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n4 = n4 + norm(V'*V - eye(m+1)) + norm(V*(V'*Vo) - Vo);
  
  n3 = n3/2; n4 = n4/4;

  param.real = 1;
  param.reorth = 1;
  param.inner_product = @(x, y) y'*x;
  xi = [xi 1+2i 1-2i];
  [Vo, Ko, Ho] = rat_krylov(struct('A', A, 'B', B), b, xi);
  [V,  K,  H]  = rat_krylov(A, B, b, xi(1:end-2), param);
  [V,  K,  H]  = rat_krylov(B\A,  V, K, H, xi(end-1:end));  
  n5 = norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n6 = norm(V'*V - eye(m+3)) + norm(V*(V'*Vo) - Vo);  
  n6 = n6/(2*m);
  
  check = [n1 n2 n3 n4 n5 n6] < tol;
end
