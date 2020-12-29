function c = times(a, b)
%TIMES    Multiplication of two RKFUNs.

  % Hand scalar multiply to mtimes.
  if     isnumeric(a),          c = mtimes(a, b);        return
  elseif length(a.coeffs) == 1, c = mtimes(a.coeffs, b); return
  elseif isnumeric(b),          c = mtimes(a, b);        return
  elseif length(b.coeffs) == 1, c = mtimes(a,b.coeffs);  return
  end

  % Prepare pencil of the product c = a*b. (mp safe)
  K = []; K(size(a.K,1)+size(b.K,2), size(a.K,2)+size(b.K,2)) = 0;    
  H = K;

  % We now compute
  %     KK = [eye(length(a.coeffs), length(a.coeffs)-1) a.coeffs]\a.K;
  %     HH = [eye(length(a.coeffs), length(a.coeffs)-1) a.coeffs]\a.H;
  % more efficiently, exploiting the structure.
  
  coeffs = a.coeffs;
  scl1   = coeffs(end);
  coeffs = -coeffs/scl1;

  KK = a.K;
  KK(1:end-1, end) = KK(1:end-1, end) + KK(end)*coeffs(1:end-1);  
  KK(end) = KK(end)/scl1;

  HH = a.H;
  if size(HH, 2) ~= 1 && abs(HH(end, end-1)) > eps(1)    
    HH(1:end-1, end-1) = ...
        HH(1:end-1, end-1) + HH(end, end-1)*coeffs(1:end-1);
    HH(end, end-1) = HH(end, end-1)/scl1;
  end
  HH(1:end-1, end) = HH(1:end-1, end) + HH(end)*coeffs(1:end-1);  
  HH(end) = HH(end)/scl1;
  
  K(1:size(a.K, 1), 1:size(a.K, 2))     = KK;
  K(size(a.K, 1):end, size(a.K, 1):end) = b.K;
  H(1:size(a.K, 1), 1:size(a.K, 2))     = HH;
  H(size(a.K, 1):end, size(a.K, 1):end) = b.H;

  % Prepare coeffs vector of the product c = a*b. (mp safe)
  coeffs = []; 
  coeffs(size(a.K, 1) + size(b.K, 2), 1) = 0; 
  coeffs(length(a.coeffs):end)           = b.coeffs;  
  
  t = type(a) + type(b);

  k = t(1) - t(2);
  c = rkfun(K, H, coeffs, k);

end
