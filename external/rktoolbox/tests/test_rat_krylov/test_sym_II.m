function check = test_sym_II
  
  N = 50;
  A = gallery('lehmer', N);
  nrm = 28;
  tol = 1e-13;
  
  b = ones(N, 1);
  b(1:2:N) = 0;
  b = b/norm(b);
  
  xi = [inf, 2+1i, 2-1i, 5-1i, 5+1i, -1, -1];
  m = length(xi);
  
  % In the following we check the backward error of the produced
  % RAD and the orthonormality of the constructed basis. The first
  % four calls are all mathematically equivalent.  
  
  [V, K, H] = rat_krylov(A, b, xi);
  n1 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n2 = norm(V'*V - eye(m+1))/m;
      
  V = rat_krylov(A, b, K, H);
  n3 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n4 = norm(V'*V - eye(m+1))/m;

  [V, K, H] = rat_krylov(A, eye(N), b, xi, 'real');
  n1 = n1 + norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n2 = n2 + norm(V'*V - eye(m+1))/m;
  
  V = rat_krylov(A, b, K, H, 'real');
  n3 = n3 + norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n4 = n4 + norm(V'*V - eye(m+1))/(2*m);
  
  
  [V, K, H] = rat_krylov(A, eye(N), b, xi);
  n5 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n6 = norm(V'*V - eye(m+1))/(m);
  
  V = rat_krylov(struct('A', 2*A, 'B', 2*eye(N)), b/2, K, H);
  n7 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));  
  n8 = norm(4*V'*V - eye(m+1))/(2*m);

  AA = gallery('tridiag', N);
  [V, K, H] = rat_krylov(AA, eye(N), b, K, H);
  n7 = n7 + norm(AA*V*K - V*H)/(4*norm(K) + norm(H));
    
  % A slightly modified example.  
  E = AA;
  param.inner_product = @(x, y) y'*((10*eye(N) + E)*x);
  AB.multiply = @(rho, eta, x)  rho*A*x-eta*x;
  AB.solve    = @(nu,  mu,  x)  (nu*A-mu*speye(N))\x;
  AB.isreal   = 1;  
  
  [V, K, H] = rat_krylov(AB, b, xi, param);
  n9 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n0 = norm(param.inner_product(V, V) - eye(m+1))/(2*m);

  b = b/sqrt(param.inner_product(b, b));
  V = rat_krylov(AB, b, K, H, param);
  n9 = n9 + norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n0 = n0 + norm(param.inner_product(V, V) - eye(m+1))/(2*m);

  check = [n1 n2 n3 n4 n5 n6 n7 n8 n9 n0] < tol;
end