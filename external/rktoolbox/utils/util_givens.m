function [s, c] = util_givens(x, y)
% UTIL_GIVENS    Computes sin and cos for a plane rotation.
%
% [s, c] = util_givens(x, y) for scalars x and y computes the sin
% and cos for the plane rotation G given by
%               |c  -s|
%               |s'  c|,
% such that (G*[x; y]) = [h; 0].

  if abs(x) < eps(1)
    c = 0; s = 1;
    return
  end

  t = -y/x;
  c = 1/sqrt(1+abs(t)^2);
  s = c*conj(t);
end
