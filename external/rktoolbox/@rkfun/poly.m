function [p, q, pq] = poly(obj)
%POLY    Convert RKFUN into a quotient of two polynomials.
%
% [p,q,pq] = poly(obj)
% Returns two vectors p and q of monomial coefficients such that 
% the RKFUN is represented as r(z) = polyval(p,z)/polyval(q,z). 
%
% The function handle pq can be used for scalar evaluation of p/q.
%
% The conversio of an RKFUN to monomial basis may be ill-conditioned 
% and multiple precision arithmetic may be required. 

  p  = poly(roots(obj));
  q  = poly(poles(obj));  
  nr = 0.594583697409052;
  nf = polyval(p,nr)/polyval(q,nr)/feval(obj,nr);
  p  = p/nf;
  pq = @(z) polyval(p,z)./polyval(q,z);
  
end