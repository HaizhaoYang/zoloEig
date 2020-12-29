function c = uminus(a)
%UMINUS    Unary minus. 
  
  c = a;
  c.coeffs = -c.coeffs;
  
end