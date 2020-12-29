function ratfun = reduce(obj,tol,A,b)
%REDUCE   Reduce the type of an RKFUN. -- not ready yet! --

  if nargin < 2, tol = 1e-15; end

  % Construct matrix A for which we have AVK = VH, with V = I.
  K = obj.K;
  H = obj.H;
  y = null(K');
  % x = null(H');                    % one option
  alph = mean(eig(K\H)); x = alph*y; % another option
  A = H/K + x*y';

  % The vector b is the first unit vector.
  b = eye(size(A,1), 1);
  V = rat_krylov(A, b, K, H);

  % Construct F explictly.
  F = feval(obj, A, eye(size(A)));

  maxit = 5;
  xi = poles(obj);

  [xi, ratfun, misfit, out] = ...
      rkfit(F, A, b, xi, maxit, tol);
  keyboard
end
