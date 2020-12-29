function flag = isreal(obj)
%ISREAL    Return true iff an RKFUN is real-valued.

  flag = isreal(obj.K) & isreal(obj.H) & isreal(obj.coeffs);

end
