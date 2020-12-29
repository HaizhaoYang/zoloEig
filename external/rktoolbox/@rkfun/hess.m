function retobj = hess(obj)
%HESS    Convert RKFUN pencil into (strict) upper-Hessenberg form.

  retobj = obj;
  [~, ~, Q, Z] = qz(obj.H(2:end, :), obj.K(2:end, :)); 
  Q = blkdiag(1, Q);
  retobj.H = Q*obj.H*Z;
  retobj.K = Q*obj.K*Z;
  retobj.coeffs = Q*obj.coeffs;
  
  retobj.H = triu(retobj.H, -1);
  retobj.K = triu(retobj.K, -1);
end
