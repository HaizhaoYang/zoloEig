function [V, W, T] = util_bilanczos(A, p, q)
% UTIL_BILANCZOS    Bi-orthogonal Lanczos.
%
% Returns bi-orthogonal V and W such that W'*A*V is tridiagonal.
% This is a very basic implementation without reorthogonalization
% or look-ahead. Use with caution.
%
% Example of usage: [V, W, T] = util_bilanczos(A, p, q);
%
% The implementation follows the algorithm on page 78 of [Anne Greenbaum,
% "Iterative Methods for Solving Linear Systems", SIAM, 1997].

  ip = @(x,y) y'*x;
  n  = size(A,1);
  V  = p/sqrt(ip(p,p));
  W  = q/(ip(q,V));
  for j = 1:n
    Av = A*V(:,j);
    Aw = A'*W(:,j);
    alph(j,1) = ip(Av,W(:,j));
    if j == n, break, end
    vv = Av - alph(j)*V(:,j);
    ww = Aw - conj(alph(j))*W(:,j);
    if j > 1
      vv = vv - bet(j-1)*V(:,j-1);
      ww = ww - gam(j-1)*W(:,j-1);
    end
    gam(j,1) = sqrt(ip(vv,vv)); V(:,j+1) = vv/gam(j);
    bet(j,1) = ip(V(:,j+1),ww); W(:,j+1) = ww/conj(bet(j));
  end

  T = A^0;
  T(1:n+1:end) = alph;
  T(2:n+1:end) = gam;
  T(n+1:n+1:end) = bet;

end