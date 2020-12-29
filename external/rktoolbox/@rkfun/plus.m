function c = plus(a, b)
%PLUS    Scalar addition. 
  
  if isnumeric(a) && isscalar(a)
    c = b; c.coeffs(1) = a + c.coeffs(1); return
  elseif isnumeric(b) && isscalar(b)
    c = a; c.coeffs(1) = c.coeffs(1) + b; return
  end

  % Prepare pencil of the sum c = a+b. (mp safe)
  K = []; K(size(a.K,1)+size(b.K,2), size(a.K,2)+size(b.K,2)) = 0;    
  H = K;

  K(1:size(a.K, 1),         1:size(a.K, 2))   = a.K;
  K([1 size(a.K, 1)+1:end], size(a.K, 1):end) = b.K;
  H(1:size(a.K, 1),         1:size(a.K, 2))   = a.H;
  H([1 size(a.K, 1)+1:end], size(a.K, 1):end) = b.H;     
    
  % Prepare coeffs vector of the sum c = a+b. (mp safe)
  coeffs = []; 
  coeffs(size(a.K, 1) + size(b.K, 2), 1) = 0; 

  coeffs(1, 1)                      = a.coeffs(1)+b.coeffs(1); 
  coeffs(2:length(a.coeffs), 1)     = a.coeffs(2:end); 
  coeffs(length(a.coeffs)+1:end, 1) = b.coeffs(2:end);
  
  t1 = type(a); 
  t2 = type(b);
  k  = max(t1(1)+t2(2), t2(1)+t1(2)) - (t1(2) + t2(2)); 
  c  = rkfun(K, H, coeffs, k);
  
end
