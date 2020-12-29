function c = power(a, n)
%POWER    Integer exponentiation of an RKFUN. 
  
  assert(isfinite(n) && round(n) == n, ...
         'POWER: Exponent needs to be an integer.')
  
  c = gallery('constant', 1);
  
  if n < 0,  a = 1./a; n = -n; end
  if n == 0, return;           end
  
  while n > 1
    if mod(n, 2) == 0, n = n/2;
    else               c = c.*a; n = (n-1)/2; end
    a = a.*a;
  end
  c = c.*a;
  
end
