function check = test_nlep()
  
  % Linearise a scalar rational function and find its roots.  
  Xi    = [-4, -5, -3, -4];
  Sigma = [1, 3];

  Nmax = 8;
  n    = 1;
  tol  = 1e-15;    
  
  A = @(z) (3*z.^2 - 2*z - 3)./(z+4)./(z+5);

  AB = util_linearise_nlep(A, Sigma, Xi, tol, Nmax, 1);
 
  D     = AB.D;
  xi    = AB.xi;
  sigma = AB.sigma;
  beta  = AB.beta;
  N     = AB.N;
    
  % Build linearisation explicitly, and
  % compare with pencil representation.
  [AN2, BN2] = AB.get_matrices();
  ee = eig(full(AN2), full(BN2));
  [~, ind] = sort(abs(ee)); ee = ee(ind);
  
  v = randn(N*n, 1); 
  rkxi = 1.5*ones(1, 2);  

  V1 = rat_krylov(AB,       v, rkxi);
  V2 = rat_krylov(AN2, BN2, v, rkxi);

  check = [norm(V1-V2)<1e-14, N == 3 , max(A(ee(1:2))) < 1e-14];
  
end