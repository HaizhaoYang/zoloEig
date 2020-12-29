function obj = vpa(obj, d)
%VPA    Convert RKFUN into variable precision format. 
  
  if nargin < 2, d = digits; end
  
  obj.K = vpa(obj.K, d);
  obj.H = vpa(obj.H, d);
  obj.coeffs = vpa(obj.coeffs, d);
  
end