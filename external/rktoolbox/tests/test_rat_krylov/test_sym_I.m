function check = test_sym_I
  
  N = 50;
  A = gallery('tridiag', N);
  nrm = 4;
  tol = 1e-14;
  
  b = ones(N, 1);
  b(1)   = 0;
  b(4)   = -1;
  b(end) = 0;
  
  xi = [1+1i, 1-1i, 5, inf, 0 ,0, inf, -1];
  m = length(xi);
  
  % In the following we check the backward error of the produced
  % RAD and the orthonormality of the constructed basis. The first
  % four calls are all mathematically equivalent.  
  
  [V, K, H] = rat_krylov(A, b, xi);
  n1 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n2 = norm(V'*V - eye(m+1));
      
  param.orth = 'CGS';
  param.real = 1;
  param.refinement = 1;
  [V, K, H] = rat_krylov(A, b, xi, param);
  n3 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n4 = norm(V'*V - eye(m+1));

  [V, K, H] = rat_krylov(A, speye(N), b, xi, 'real');
  n5 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n6 = norm(V'*V - eye(m+1));
 
  [V, K, H] = rat_krylov(struct('A', A), b, xi, 'real');
  n7 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n8 = norm(V'*V - eye(m+1));

  % A slightly modified example.
  
  E = diag(rand(N, 1));
  param.inner_product = @(x, y) y'*((10*eye(N) + E)*x);
  AB.multiply = @(rho, eta, x)  rho*A*x-eta*x;
  AB.solve    = @(nu,  mu,  x)  (nu*A-mu*speye(N))\x;
  AB.isreal   = 1;  
  
  [V, K, H] = rat_krylov(AB, b, xi, param);
  n9 = norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n0 = norm(param.inner_product(V, V) - eye(m+1));

  param.real = 0;
  param.refinement = 0;
  [V, K, H] = rat_krylov(AB, b, xi, param);
  n9 = n9 + norm(A*V*K - V*H)/(nrm*norm(K) + norm(H));
  n0 = n0 + norm(param.inner_product(V, V) - eye(m+1));

  check = [n1 n2 n3 n4 n5 n6 n7 n8 n9/2 n0/2] < tol;
end