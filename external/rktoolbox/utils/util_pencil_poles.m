function xi = util_pencil_poles(K, H)
% UTIL_PENCIL_POLES    Computes the poles of the pencil (H, K).
%
% xi = util_pencil_poles(K, H) returns the poles of the
% (m+1)-by-m pencil (H, K), i.e, it returns the eigenvalues of the
% lower m-by-m part.
  
  m = size(H, 2); xi = zeros(0,1);
  j = 1;
  while j <= m
    if j+2 <= m+1
      if double(H(j+2, j)) == 0
        xi(1,j) = H(j+1, j)/K(j+1, j);
      else
	[HH, KK] = qz(H(j+1:j+2, j:j+1), K(j+1:j+2, j:j+1));
        cxi = diag(HH)./diag(KK);
	if not(or(isinf(cxi(1)), isreal(cxi(1))))
	  cxi(1) = (cxi(1)+conj(cxi(2)))/2; 
	  cxi(2) = conj(cxi(1)); 
	end
        xi(1, j)   = cxi(1);
        xi(1, j+1) = cxi(2);
        j = j+1;
      end
    else
      xi(1,j) = H(j+1, j)/K(j+1, j);
    end
    j = j+1;
  end

end
