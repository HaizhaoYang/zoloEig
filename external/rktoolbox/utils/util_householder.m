function [w, sigma, gamma] = util_householder(x, variant)
% UTIL_HOUSEHOLDER    Computes the representation for a Householder
%                     reflector.
%
% For an n-vector x
%     [w, sigma, gamma] = util_householder(x);
% returns an n-vector w and scalars sigma and gamma, representing
% a Householder reflector
%     P = eye(n) - sigma*w*w'; 
% such that
%     P'*x = gamma*norm(x)*eye(n, 1).
% gamma is of unit module.  
%  
% The implementation of all the four variants from [1] are
% supported, and one can be chosen by the optional variable
% variant. 
% 
% variant 1 [Wilkinson approach]
%     The scalar sigma is guaranteed to be real, and hence P is
%     both unitary  and Hermitian. gamma may by complex.
% variant 2 (default)
%     The scalar sigma is guaranteed to be real and moreover in the
%     interval [.5, 1], hence P is both unitary and
%     Hermitian. gamma may by complex.  
% variant 3 [Hammarling, Du Croz approach]
%     The scalar sigma may be complex, hence P may be non
%     Hermitian. gamma equals to -1.
% variant 4 [LAPACK approach]
%     The scalar sigma may be complex, hence P may be non
%     Hermitian. gamma equals to -1.     
% 
% References
%
% [1] R.B. Lehoucq. The computation of elementary unitary matrices,
% LAPACK Working Note 72, 1995.  
  
  if nargin == 1, variant = 2; end
  switch variant
   case 1
    [w, sigma, gamma] = householder01(x);
   case 2
    [w, sigma, gamma] = householder02(x);
   case 3
    [w, sigma, gamma] = householder03(x);
   case 4
    [w, sigma, gamma] = householder04(x);
   otherwise
    error('Incorrect variant in util_householder.')
  end
  
end



function [w, sigma, gamma] = householder01(x)
% HOUSEHOLDER01    Construct an elementary unitary matrix.
% 
% For an n-vector x
%     [w, sigma, gamma] = householder01(x);
% constructs an n-vector w, and scalars sigma and gamma so that
%     P'*x = gamma*norm(x)*eye(n, 1),
% where 
%     P = eye(n) - sigma*w*w';
% P is unitary and Hermitian. The scalar sigma is real, and gamma
% is unimodular. The implementation is based on [1, Section 2.1].
% 
% References
%
% [1] R.B. Lehoucq. The computation of elementary unitary matrices,
% LAPACK Working Note 72, 1995.  

  scale = sum(abs(real(x)))+sum(abs(imag(x)));
  x = x/scale;
  
  norm_x = norm(x);
  rho    = abs(x(1));
  
  w = x;
  
  if isreal(x) 
    if x(1) >=0 , gamma = -1; else gamma = 1; end    
  else
    theta = angle(x(1));
    gamma = -exp(theta*1i);
  end
  
  w(1) = w(1) - gamma*norm_x;   
  sigma = 1/(norm_x*(rho+norm_x));

end



function [w, sigma, gamma] = householder02(x)
% HOUSEHOLDER02    Construct an elementary unitary matrix.
% 
% For an n-vector x
%     [w, sigma, gamma] = householder02(x);
% constructs an n-vector w, and scalars sigma and gamma so that
%     P'*x = gamma*norm(x)*eye(n, 1),
% where 
%     P = eye(n) - sigma*w*w';
% P is unitary and Hermitian. The scalar sigma is real, and gamma
% is unimodular. The implementation is based on [1, Section 2.2].
% 
% References
%
% [1] R.B. Lehoucq. The computation of elementary unitary matrices,
% LAPACK Working Note 72, 1995.  
  
  norm_x = norm(x);
  rho    = abs(x(1));

  if isreal(x) 
    if x(1) >=0 , gamma = -1; else gamma = 1; end
    w = -gamma*x/norm_x; 
  else
    theta = angle(x(1));
    gamma = -exp(theta*1i);
    w = x*exp(-theta*1i)/norm_x;
  end
  
  w(1) = w(1) + 1;  
  sigma = norm_x/(rho+norm_x);
  
end



function [w, sigma, gamma] = householder03(x)
% HOUSEHOLDER03    Construct an elementary unitary matrix.
% 
% For an n-vector x
%     [w, sigma, gamma] = householder03(x);
% constructs an n-vector w, and scalars sigma and gamma so that
%     P'*x = gamma*norm(x)*eye(n, 1),
% where 
%     P = eye(n) - sigma*w*w';
% P is unitary, but not necessarily Hermitian. It is symmetric in
% the real case. The scalar sigma may be complex, and gamma is
% unimodular. The implementation is based on [1, Section 2.3].
% 
% References
%
% [1] R.B. Lehoucq. The computation of elementary unitary matrices,
% LAPACK Working Note 72, 1995.
  
  norm_x = norm(x);
  gamma  = -1;
  
  nu = sign(real(x(1)))*norm_x;
  if nu == 0, nu = norm_x; end 
  kappa  = (abs(real(x(1)))+norm_x)/norm_x;
  
  w = x; 
  w(1) = w(1) + nu;
  w = (sqrt(kappa)/(x(1)+nu))*w;  
  
  sigma = (x(1)+nu)/(nu*kappa);  
end



function [w, sigma, gamma] = householder04(x)
% HOUSEHOLDER04    Construct an elementary unitary matrix.
% 
% For an n-vector x
%     [w, sigma, gamma] = householder04(x);
% constructs an n-vector w, and scalars sigma and gamma so that
%     P'*x = gamma*norm(x)*eye(n, 1),
% where 
%     P = eye(n) - sigma*w*w';
% P is unitary, but not necessarily Hermitian. It is symmetric in
% the real case. The scalar sigma may be complex, and gamma is
% unimodular. The implementation is based on [1, Section 2.4].
% 
% References
%
% [1] R.B. Lehoucq. The computation of elementary unitary matrices,
% LAPACK Working Note 72, 1995.
     
  norm_x = norm(x);
  gamma  = -1;  
  
  nu = sign(real(x(1)))*norm_x;  
  if nu == 0, nu = norm_x;, end  
  
  w = x;
  w(1) = w(1) + nu;
  w = w/(x(1)+nu);
  
  sigma = (x(1)+nu)/nu;
  
end



function [v, beta] = util_householder_old(x)
% [v, beta] = util_householder_old(x) return vectors v and a scalar
% beta, representing a Householder reflector
% P = eye(length(x)) - beta*v*v'; such that
% P*x = norm(x, 'fro')*eye(length(x), 1).

  if norm(x, inf) == 0
    error('UTIL_HOUSEHOLDER: Zero vector encountered.');
  end

  if isreal(x)
    sigma(1) = x(2:end)'*x(2:end);
    v = [1; x(2:end)];

    if     sigma(1) == 0 && x(1) > 0, beta = 0; 
    elseif sigma(1) == 0 && x(1) < 0, beta = -2; 
    else
      mu = sqrt(sigma(1) + x(1)^2);
      if x(1) < 0, v(1) = x(1) - mu;
      else         v(1) = -sigma(1)/(x(1) + mu); end
      beta = 2*v(1)^2/(sigma(1) + v(1)^2);
      v = v/v(1);
    end

  else
    
    x = x/norm(x, 'inf');

    th = angle(x(1));
    v  = x;
    c  = exp(1i*th)*norm(x);

    if abs(v(1) + c) < abs(v(1) - c)
      v(1) = v(1) - c; 
    else
      v(1) = v(1) + c; 
    end
    beta = 2/(v'*v);
  end

end % function
