function obj = mp(obj,d)
%MP    Convert RKFUN into Advanpix Multiple Precision format.
%
% This requires the ADVANPIX Multiprecision Computing Toolbox to be
% installed.
  
  if nargin < 2, d = mp.Digits; end
  
  if isa(obj.coeffs, 'sym')
    warning(['RKFUN:MP: Conversion from VPA to MP is not supported.' ...
             ' Going via double!']);
    obj = double(obj);
  end
  
  obj.K = mp(obj.K, d);
  obj.H = mp(obj.H, d);
  obj.coeffs = mp(obj.coeffs, d);
  
end