function c = mtimes(a, b)
%MTIMES    Scalar multiplication. 
   
  if isnumeric(a) && isscalar(a)
    c = b; c.coeffs = a*c.coeffs; return
  elseif isnumeric(b) && isscalar(b)
    c = a; c.coeffs = b*c.coeffs; return
  end
  
  if isa(a,'rkfun') && isa(b,'rkfun')
     error(['MTIMES: Use .* for point-wise multiplication ' ...
            'of RKFUNs.']);     
  else   
     error('MTIMES: This operation is not supported.');
  end
    
end