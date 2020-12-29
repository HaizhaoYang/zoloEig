function check = test_sym_III
  
  N = 51;
  A = sparse(diag(linspace(10, 15, N)));
  B = 2*speye(N);  
  B([1 26 51], [1 26 51]) = reshape([7 1:8], 3, 3);

  nrmA = 46;
  nrmB = 15;
  tol = 1e-14;
  
  b = ones(N, 1);
  b(1:5:N)  = 0.5;
  b(1:10:N) = 2;
  
  xi = [inf -20 -10 10 20 0];
  m = length(xi);
    
  [V, K, H] = rat_krylov(A, B, b, xi);
  n1 = norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n2 = norm(V'*V - eye(m+1));

  xi = [xi 1+2i 1-2i];
  [V, K, H] = rat_krylov(A, B, b, xi, 'real');   
  n1 = n1 + norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n2 = n2 + norm(V'*V - eye(m+3));

  param.reorth = 0;
  param.real = 0;
  param.refinement = 1;
  xi = xi(1:end-2);
  [V, K, H] = rat_krylov(A, B, b, xi, param);
  n1 = n1 + norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n2 = n2 + norm(V'*V - eye(m+1))*1e-1;
  
  n1 = n1/3; n2 = n2/3;
  
  [Vo, Ko, Ho] = rat_krylov(A, B, b, xi);
  [V,  K,  H]  = rat_krylov(A, B, b, xi(1));
  for j = 2:m
    [V, K, H] = rat_krylov(A, B, V, K, H, xi(j));
  end
  n3 = norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n4 = norm(V'*V - eye(m+1)) + norm(V*(V'*Vo) - Vo);
  
  [V, K, H] = rat_krylov(A, B, b, xi(1:end-1));
  [V, K, H] = rat_krylov(A, B, V, K, H, xi(end));
  n3 = n3 + norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n4 = n4 + norm(V'*V - eye(m+1)) + norm(V*(V'*Vo) - Vo);
  
  n3 = n3/2; n4 = n4/4;

  param.real = 1;
  param.reorth = 1;
  xi = [xi 1+2i 1-2i];
  [Vo, Ko, Ho] = rat_krylov(struct('A', A, 'B', B), b, xi, 'real');
  [V,  K,  H]  = rat_krylov(A, B, b, xi(1:end-2), param);
  [V,  K,  H]  = rat_krylov(B\A,  V, K, H, xi(end-1:end), 'real');  
  n5 = norm(A*V*K - B*V*H)/(nrmA*norm(K) + nrmB*norm(H));
  n6 = norm(V'*V - eye(m+3)) + norm(V*(V'*Vo) - Vo);
  
  n6 = n6/(2*m);
  
  check = [n1 n2 n3 n4 n5 n6] < tol;
end
