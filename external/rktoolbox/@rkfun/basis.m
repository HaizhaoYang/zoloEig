function b = basis(obj)
%BASIS    Rational basis functions of an RKFUN.
%
% Returns a cell array of RKFUNs that represent obj, i.e., 
% obj = obj.coeffs(1)*b{1} + obj.coeffs(2)*b{2} + ...
% Note that the RKFUNs returned by basis do not necessarily
% need to be linearly independent.

  for j = 1:length(obj.coeffs)
    Kj = obj.K(1:j, 1:j-1);
    Hj = obj.H(1:j, 1:j-1);
    b{j} = rkfun(Kj, Hj, [zeros(j-1, 1); 1]);
  end

end
