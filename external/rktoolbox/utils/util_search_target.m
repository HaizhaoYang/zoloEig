function [V, K, H, S, T] = util_search_target(A, b, xi, k)
% UTIL_SEARCH_TARGET    Constructs search and target spaces for RKFIT.
%
% The function construct the search and target spaces used by
% RKFIT, see [1, Algorithm 3.1]. This is only a limited version,
% and it is not used by RKFIT itself.
%
% [V, K, H, S, T] = util_search_target(A, b, xi, k)
%
% On exit we have A*V*K = V*H. S is an orthonormal basis of the
% rational Krylov space representing the search space, having first
% column collinear to b. T is an orthonormal basis of the polynomial
% Krylov space representing the target space.
%
% References:
%
% [1] M. Berljafa and S. G{\"u}ttel. The RKFIT algorithm for
%     nonlinear rational approximation, MIMS EPrint 2015.38,
%     Manchester Institute for Mathematical Sciences, The
%     University of Manchester, UK, 2015.

  m = length(xi);

  if k > 0, xi = [xi(:); inf*ones(k, 1)]'; end

  [V, K, H] = rat_krylov(A, b, xi);
  [~, ~, Q] = util_hh2th(K, H);

  S = V(:, 1:m+1);
  T = V*Q(1:m+1+k, :)';

end
