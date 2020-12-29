function check = test_unsym_II
  
    
  N = 60;
  A = gallery('dorr', N);
  nrm = 200;
  tol = 2e-13;
  
  b = ones(N, 1);
  b(2:3:N) = 0;
  b = b/norm(b);
  
  xi = 5:5:25; xi = -xi;
  m = length(xi);
  
  % In the following we check the backward error of the produced
  % RAD and the orthonormality of the constructed basis. The first
  % four calls are all mathematically equivalent.  
  
  [V, K, H] = rat_krylov(A, b, xi);
  n1 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n2 = norm(V'*V - eye(m+1))/m;
      
  V = rat_krylov(A, b, K, H);
  n3 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n4 = norm(V'*V - eye(m+1))/(m*m);

  [V, K, H] = rat_krylov(A, eye(N), b, xi);
  n1 = n1 + norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n2 = n2 + norm(V'*V - eye(m+1))/m;    
  
  V = rat_krylov(A, b, K, H);
  n3 = n3 + norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n4 = n4 + norm(V'*V - eye(m+1))/(m*m);
  
  xi = 2:5:25;
  xi = util_cplxpair([xi+5i xi-5i]);
  m = length(xi);
  param.real = 1;
  param.refinement = 1;
  
  [V, K, H] = rat_krylov(A, eye(N), b, xi, param);
  n5 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n6 = norm(V'*V - eye(m+1))/(m);
  
  V = rat_krylov(struct('A', A, 'B', eye(N)), b, K, H, param);
  n7 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));  
  n8 = norm(V'*V - eye(m+1))/(m*m);

  AA = gallery('tridiag', N);
  [V, K, H] = rat_krylov(AA, eye(N), b, K, H, param);
  n7 = n7 + norm(AA*V*K - V*H)/((4*norm(K) + norm(H))*norm(V));
    
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
  n0 = n0 + norm(param.inner_product(V, V) - eye(m+1))/(m*m);
    
  check = [n1 n2 n3 n4 n5 n6 n7 n8*1e-2 n9 n0*1e-2] < tol;
  
end
