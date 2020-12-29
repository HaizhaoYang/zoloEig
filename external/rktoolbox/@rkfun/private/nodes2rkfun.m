function r = nodes2rkfun(rts, pls, scl)
%NODES2RKFUN    Construct an RKFUN from its roots and poles.

  if nargin < 3, scl = 1;  end
  if nargin < 2, pls = []; end

  rts = rts(isfinite(rts));
  pls = pls(isfinite(pls));   
  
  k = length(rts) - length(pls);
  m = max(length(rts), length(pls));
  
  if m == 0
    r = rkfun(zeros(1, 0), zeros(1, 0), scl, k);
    return
  end
  
  [rts, f1] = util_cplxpair(rts);
  [pls, f2] = util_cplxpair(pls);
  
  rts(length(rts)+1:m) = inf;
  pls(length(pls)+1:m) = inf;  
  
  K = []; H = []; c = [];
  
  if not(f1 && f2)
    % May produce a complex-valued pencil.
    for j = 1:m 
      [K1, H1, c1] = nodes2rkfun11(rts(j), pls(j));     
      [K,  H,  c]  = merge(K, H, c, K1, H1, c1);
    end            
  else      
    % Produces a real-valued pencil.
    j = 1;
    inR = @(x)  isreal(x) && isfinite(x);
    inC = @(x) ~isreal(x);                % in C\R !
    while j <= m
      if isreal(rts(j)) && isreal(pls(j))
        [K1, H1, c1] = nodes2rkfun11(rts(j), pls(j));
        j = j-1;
        
      elseif inC(rts(j)) && inC(pls(j))
        [K1, H1, c1] = ...
            nodes2rkfun22(rts(j), rts(j+1), pls(j), pls(j+1));
        
      elseif inC(rts(j)) && inR(pls(j)) && inR(pls(j+1))        
        [K1, H1, c1] = ...
            nodes2rkfun22(rts(j), rts(j+1), pls(j), pls(j+1));
        
      elseif inC(rts(j)) && inR(pls(j)) && isinf(pls(j+1))      
        [K1, H1, c1] = nodes2rkfun21(rts(j), pls(j));   

      elseif inC(rts(j)) && isinf(pls(j))
        [K1, H1, c1] = nodes2rkfun20(rts(j));

      elseif inR(rts(j)) && isinf(rts(j+1)) && inC(pls(j))
        [K1, H1, c1] = nodes2rkfun12(rts(j), pls(j));

      elseif inR(rts(j)) && inR(rts(j+1)) && inC(pls(j))
        [K1, H1, c1] = ...
            nodes2rkfun22(rts(j), rts(j+1), pls(j), pls(j+1));

      elseif isinf(rts(j)) && isinf(rts(j+1)) && inC(pls(j))
        [K1, H1, c1] = nodes2rkfun02(pls(j));

      else
        error('NODES2RKFUN: ... may have been attacked by a bug!')      
      end % if
      
      j = j+2;
      [K,  H,  c]  = merge(K, H, c, K1, H1, c1);
    end % while
  end
  
  r = rkfun(K, H, c*scl, k);
end


function [K, H, c] = merge(K1, H1, c1, K2, H2, c2)
% MERGE    Concatenates two pencil representation in one, equivalent to
%          multiplying rational functions.
%
% [K, H, c] = merge(K1, H1, c1, K2, H2, c2) expands (K1, H1) by
% (K2, H2) in order to perform multiplication of the rational
% functions represented by (K1, H1, c1) and (K2, H2, c2). The
% coefficient vectors cj are assumed to be scalar multiples
% of the size(Kj, 1)-th canonical vector.
% 
% See also rkfun/times.
  
  assert(norm(c1(1:end-1), inf) == 0);
  assert(norm(c2(1:end-1), inf) == 0);
  
  if length(K1) == 0
    K = K2; H = H2; c = c2;
    return
  end
  
  K = []; K(size(K1,1)+size(K2,2), size(K1,2)+size(K2,2)) = 0;  
  H = K;
  
  K(1:size(K1, 1), 1:size(K1, 2))     = K1;
  K(size(K1, 1):end, size(K1, 1):end) = K2;
  H(1:size(K1, 1), 1:size(K1, 2))     = H1;
  H(size(K1, 1):end, size(K1, 1):end) = H2;
  
  c = []; 
  c(size(K, 1), 1) = 0; 
  if   length(c2) ~= 0, c(length(c1):end) = c2*c1(end);
  else                  c = c1;                         end
 
end


function [K, H, c] = nodes2rkfun11(zeta, xi)
%NODES2RKFUN11    Construct pencil representation of a rational
%                 function of exact type (0, 1), (1, 0) or (1, 1).
%
% [K, H, c] = nodes2rkfun11(zeta, xi) constructs the pencil
% representation for (x-zeta)/(x-xi) if both zeta and xi are
% finite. Otherwise, for (x-zeta) if xi infinite, or for 1/(x-xi) if
% zeta is infinite.

  assert(isfinite(zeta) || isfinite(xi), ...
         'NODES2RKFUN11: type (0, 0) not supported.')
  
  if isinf(xi),       K = [1;  0];  H = [zeta;   1];
  elseif isinf(zeta), K = [0;  1];  H = [1;     xi];    
  else,               K = [-1; 1];  H = [-zeta; xi]; end
  
  c = [0; 1];
  
end


function [K, H, c] = nodes2rkfun02(xi)
%NODES2RKFUN02    Construct pencil representation of a rational
%                 function of exact type (0, 2) with a pair of
%                 complex-conjugate poles.
%
% [K, H, c] = nodes2rkfun02(xi) constructs the pencil
% representation for 1/((x-xi)(x-conj(xi))), where it is assumed
% that imag(xi) is nonzero.
  
  assert(imag(xi) ~= 0, ...
         'NODES2RKFUN02: xi cannot be a real scalar.')
  
  K = [0 0; eye(2)];
  H = [1 0; real(xi) imag(xi); -imag(xi) real(xi)];
  c = [0; 0; 1/imag(xi)];
  
end


function [K, H, c] = nodes2rkfun12(zeta, xi)
%NODES2RKFUN12    Construct pencil representation of a rational
%                 function of exact type (1, 2) with a pair of
%                 complex-conjugate poles.
%
% [K, H, c] = nodes2rkfun12(zeta, xi) constructs the pencil
% representation for (x-zeta)/((x-xi)(x-conj(xi))), where it is assumed
% that imag(xi) is nonzero, and zeta is real.
  
  assert(imag(xi) ~= 0, ...
         'NODES2RKFUN12: xi cannot be a real scalar.')
  assert(isfinite(zeta) && isreal(zeta), ...
         'NODES2RKFUN12: zeta must be a real number.')
  
  K = [0 0; eye(2)];
  H = [1 0; real(xi) imag(xi); -imag(xi) real(xi)];
  
  if abs(zeta - real(xi)) < abs(zeta)*eps(1)
    P = [0 1; 1 0];
    H = blkdiag(1, P)*H*P;
    c = [0; 0; 1];
    return
  end
  
  a1 = 0;
  a2 = 1;
  a3 = (real(xi)-zeta)/imag(xi);

  if abs(a2) > abs(a3)      
    [a2, a3] = deal(a3, a2);    
    P = [0 1; 1 0];
    H = blkdiag(1, P)*H*P;
  end
  
  a1 = -a1/a3;    
  a2 = -a2/a3;    

  K(1:2, :) = K(1:2, :) + [a1 a1; a2 a2]*diag(K(3, :));
  H(1:2, :) = H(1:2, :) + [a1 a1; a2 a2]*diag(H(3, :));
  
  c = [0; 0; a3];
  
end


function [K, H, c] = nodes2rkfun20(zeta)
%NODES2RKFUN20    Construct pencil representation of a rational
%                 function of exact type (2, 0) with a pair of
%                 complex-conjugate roots.
%
% [K, H, c] = nodes2rkfun20(zeta) constructs the pencil
% representation for (x-zeta)(x-conj(zeta)), where it is assumed
% that imag(zeta) is nonzero.
  
  assert(imag(zeta) ~= 0, ...
         'NODES2RKFUN20: zeta cannot be a real scalar.')
  assert(isfinite(zeta),  ...
         'NODES2RKFUN20: zeta must be finite.')  
  
  K = [eye(2); 0 0];
  H = [0 0; eye(2)];

  a1 = abs(zeta)^2;
  a2 = -2*real(zeta);
  a3 = 1;

  a1 = -a1/a3;    
  a2 = -a2/a3;
 
  K(1:2, :) = K(1:2, :) + [a1 a1; a2 a2]*diag(K(3, :));
  H(1:2, :) = H(1:2, :) + [a1 a1; a2 a2]*diag(H(3, :));
    
  c = [0 0 a3].';
  
end


function [K, H, c] = nodes2rkfun21(zeta, xi)
%NODES2RKFUN21    Construct pencil representation of a rational
%                 function of exact type (2, 1) with a pair of
%                 complex-conjugate roots, and a real pole.
%
% [K, H, c] = nodes2rkfun21(zeta, xi) constructs the pencil
% representation for (x-zeta)(x-conj(zeta))/(x-xi), where it is
% assumed that imag(zeta) is nonzero and xi is real.
  
  assert(imag(zeta) ~= 0, ...
         'NODES2RKFUN21: zeta cannot be a real scalar.')
  assert(isfinite(zeta) && isreal(xi) && isfinite(xi), ...
         'NODES2RKFUN21: zeta must be finite and xi real.')
    
  K = [1 0; 0 0; 0 1];
  H = [0 1; 1 0; 0 xi];

  a1 = xi - 2*real(zeta);
  a2 = 1;
  a3 = abs(zeta)^2 + a1*xi;
  
  if abs(a2) > abs(a3)         
    [a2, a3] = deal(a3, a2);    
    P = [0 1; 1 0];
    K = blkdiag(1, P)*K*P;
    H = blkdiag(1, P)*H*P;
  end
  
  a1 = -a1/a3;    
  a2 = -a2/a3;

  K(1:2, :) = K(1:2, :) + [a1 a1; a2 a2]*diag(K(3, :));
  H(1:2, :) = H(1:2, :) + [a1 a1; a2 a2]*diag(H(3, :));  
  
  c = [0 0 a3].';
  
end


function [K, H, c] = nodes2rkfun22(zeta1, zeta2, xi1, xi2)
%NODES2RKFUN22    Construct pencil representation of a rational
%                 function of exact type (2, 2) with a pair of
%                 complex-conjugate poles or roots.
%
% [K, H, c] = nodes2rkfun22(zeta1, zeta2, xi1, xi2) constructs the
% pencil representation for
% ((x-zeta1)(x-zeta2))/((x-xi)(x-conj(xi))), for the following
% cases:
%   - xi1 = conj(xi2), and zeta1 = conj(zeta2);
%   - xi1 = conj(xi2), and both zeta1 and zeta2 are real scalars;
%   - xi1 and xi2 are both real scalars, and zeta1 = conj(zeta2).
  
  assert(xi2 == conj(xi1) || and(isreal(xi1), isreal(xi2)), ...  
         ['NODES2RKFUN22: poles do not allow for a real-valued' ...
          ' representation.'])
  
  assert(zeta2 == conj(zeta1) || and(isreal(zeta1), isreal(zeta2)), ...
         ['NODES2RKFUN22: roots do not allow for a real-valued' ...
          ' representation.'])
  
  assert(~isreal(xi1) || ~isreal(xi2) || ~isreal(zeta1) || ~isreal(zeta2), ...
         'NODES2RKFUN22: case not supported.')
  
  if xi2 == conj(xi1) && imag(xi2) ~= 0
    % Complex-conjugate poles.
    xi = xi1;
    
    K = [0 0; eye(2)];
    H = [1 0; real(xi) imag(xi); -imag(xi) real(xi)];
    
    p = real(zeta1+zeta2);
    q = real(zeta1*zeta2);
    
    a1 = 1;
    a2 = 2*real(xi) - p;
    a3 = (q - abs(xi)^2 + a2*real(xi))/imag(xi);        
    
    if abs(a2) > abs(a3)  
      [a2, a3] = deal(a3, a2);  
      P = [0 1; 1 0];
      K = blkdiag(1, P)*K*P;
      H = blkdiag(1, P)*H*P;
    end
  else
    % Complex-conjugate roots.    
    zeta = zeta1;
    
    % Have to see when it is better to change them...
    if abs(abs(zeta)^2 + xi2^2 - 2*real(zeta)*xi2) > ...
          abs(abs(zeta)^2 + xi1^2 - 2*real(zeta)*xi1)
      [xi1, xi2] = deal(xi2, xi1);
    end
    
    K = [0 0; eye(2)];
    H = [1 0; xi1 1; 0 xi2];
    
    a1 = 1;
    a2 = xi1+xi2 - 2*real(zeta);
    a3 = abs(zeta)^2 + xi2^2 - 2*real(zeta)*xi2;
  end

  a1 = -a1/a3;    
  a2 = -a2/a3;
  
  K(1:2, :) = K(1:2, :) + [a1 a1; a2 a2]*diag(K(3, :));
  H(1:2, :) = H(1:2, :) + [a1 a1; a2 a2]*diag(H(3, :));
  
  c = [0 0 a3].';
  
end
