function [vf, V, K, H] = util_markovfunmv(A, b, m, fm, Xi)   
% UTIL_MARKOVFUNMV    Black-box rational Arnoldi method for
%                     Markov matrix functions based on [1].
%
% [vf, V, K, H] = util_markovfunmv(A, b, m, fm, Xi) computes an RAD of order m
% with poles chosen from Xi in order to obtain a good approximation vf to
% fm(A)*b from the corresponding rational Krylov space. Here fm has to be a
% function handle to a matrix function.
%   
% The algorithm is mathematically equivalent to the one given is [1], but the
% implementation is simpler as it is built on top the RKToolbox. 
%
% Example: Approximate sqrtm(A)*b by order 15 rational Arnoldi approximant
%
%  A = gallery('tridiag', 100); b = randn(100, 1); 
%  m = 15; fm = @(X) sqrtm(X);
%  [vf, V, K, H] = util_markovfunmv(A, b, m, fm);
%  norm(sqrtm(full(A))*b - vf)/norm(sqrtm(full(A))*b) % rel. error
%
% References:
%
% [1] S. Guettel and L. Knizhnerman. A black-box rational Arnoldi variant 
%     for Cauchy--Stieltjes matrix functions, BIT, 53 (2013), pp. 595--616.
   
  if nargin < 5,
      Xi = -logspace(-8,8,1601);
  end

  V = b; 
  K = zeros(1, 0); 
  H = zeros(1, 0); 
  
  xi = inf(1, 1);
  
  j = 1;
  while j <= m    
    [V, K, H] = rat_krylov(A, V, K, H, xi);
  
    [Q, ~] = qr(K);    
    r = rkfun(K, H, Q(:, end));
    [~, index] = min(abs(r(Xi)));
    xi = Xi(index);
   
    Am = K\H; 
    vf = V*(K*fm(Am)*(K\(norm(b)*eye(j+1, 1))));
    j = j+1;
  end  
end
