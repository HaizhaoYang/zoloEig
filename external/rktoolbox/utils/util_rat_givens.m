function [s, c] = util_rat_givens(h, k, xi)
% UTIL_RAT_GIVENS    Computes sin and cos for a plane rotation that
%                    changes a pole.
%
% [s, c] = util_rat_givens(h, k, xi) for vectors h and k of
% length 2 and xi a complex number or infinity, computes the
% sin and cos for the plane rotation G given by
%               |c  -s|
%               |s'  c|,
% for which the ratio (G*h)(2)/(G*k)(2) is set to xi.

  if isinf(xi)
    tmp = k; k = h; h = tmp;
    xi  = 0;
  end
  
  if abs(xi - h(1)/k(1)) < eps(1)*abs(xi)
    s = 1;
    c = 0;
    return
  end
  
  t = ((xi*k(2)-h(2))/(h(1)-xi*k(1)));
  c = 1/sqrt(1+abs(t)^2);
  s = c*conj(t);
end
