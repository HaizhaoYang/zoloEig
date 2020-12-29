function c = mrdivide(a, b)
%MRDIVIDE    Scalar division. 

  if isnumeric(b) && isscalar(b)
    c = a;
    c.coeffs = c.coeffs/b;
    return
  end

  error(['MRDIVIDE: This form of division is not implemented. ' ...
         'Try using ''./''.'])
  
end

