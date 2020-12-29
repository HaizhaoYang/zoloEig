function check = test_householder_real()
  N   = 80;
  tol = 1e-15;

  A = gallery('tridiag', N);
  b = ones(N, 1); 
  x = A*b;

  % Case x(1) = 0.
  x(1) = 0;
  nrm = norm(x);
  [w, sigma, gamma] = util_householder(x, 1);
  Px = x-((w'*x)*sigma)*w;
  n1 = norm(imag(sigma))/nrm;
  n2 = norm(Px(2:end))/nrm;
  n3 = (Px(1)-gamma*norm(x))/nrm;

  [w, sigma, gamma] = util_householder(x, 2);
  Px = x-((w'*x)*sigma)*w;
  n4 = norm(imag(sigma))/nrm;
  n5 = norm(Px(2:end))/nrm;
  n6 = (Px(1)-gamma*norm(x))/nrm;

  [w, sigma, gamma] = util_householder(x, 3);  
  Px = x-((w'*x)*conj(sigma))*w;
  n7 = norm(Px(2:end))/nrm;
  n8 = (Px(1)-gamma*norm(x))/nrm;

  [w, sigma, gamma] = util_householder(x, 4);
  Px = x-((w'*x)*conj(sigma))*w;
  n9 = norm(Px(2:end))/nrm;
  n0 = (Px(1)-gamma*norm(x))/nrm;

  check1 = [n1 n2 n3 n4 n5 n6 n7 n7 n9 n0] < tol;

  % Case x(1) ~= 0.
  x(1) = 0.2;
  nrm = norm(x);

  [w, sigma, gamma] = util_householder(x, 1);
  Px = x-((w'*x)*sigma)*w;
  n1 = norm(imag(sigma))/nrm;
  n2 = norm(Px(2:end))/nrm;
  n3 = (Px(1)-gamma*norm(x))/nrm;

  [w, sigma, gamma] = util_householder(x, 2);
  Px = x-((w'*x)*sigma)*w;
  n4 = norm(imag(sigma))/nrm;
  n5 = norm(Px(2:end))/nrm;
  n6 = (Px(1)-gamma*norm(x))/nrm;

  [w, sigma, gamma] = util_householder(x, 3);
  
  Px = x-((w'*x)*conj(sigma))*w;
  n7 = norm(Px(2:end))/nrm;
  n8 = (Px(1)-gamma*norm(x))/nrm;

  [w, sigma, gamma] = util_householder(x, 4);
  Px = x-((w'*x)*conj(sigma))*w;
  n9 = norm(Px(2:end))/nrm;
  n0 = (Px(1)-gamma*norm(x))/nrm;
    
  check2 = [n1 n2 n3 n4 n5 n6 n7 n7 n9 n0] < tol;
  
  check = [check1 check2];
  
end
