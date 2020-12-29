function c = rdivide(a, b)
%RDIVIDE    Division of two RKFUNs. 
  
  if isnumeric(b) && isscalar(b), c = mrdivide(a, b);
  elseif length(b.coeffs) == 1,   c = mrdivide(a, b.coeffs);
  else
    if isreal(b), flag = 'real';
    else          flag = 'complex'; end
    
    [KT, HT, QT, ZT] = ...
        move_poles_impl(b.K, b.H, b.coeffs, flag);
    
    [t1, t2] = type(b); k = t2 - t1; 
    invb = rkfun(KT, HT, QT(:, 1)/norm(b.coeffs), k);

    % work out scaling
    %point = [roots(invb).' poles(invb).'];
    %point = max(point(isfinite(point))) + 1.0;
    point = -0.468615581100624;
    scale = b.feval(point)*invb.feval(point);
    
    invb = rkfun(invb.K, invb.H, invb.coeffs/scale, k);
    c = a.*invb;
  end

end
