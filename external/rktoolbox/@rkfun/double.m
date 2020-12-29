function obj = double(obj)
%DOUBLE    Convert RKFUN into double precision (undo vpa or mp). 
  
  obj.K = double(obj.K);
  obj.H = double(obj.H);
  obj.coeffs = double(obj.coeffs);
  
end