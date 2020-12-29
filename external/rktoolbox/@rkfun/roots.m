function rts = roots(obj, varargin)
%ROOTS    Compute the roots of an RKFUN. 
%
% rts = roots(ratfun) 
% Returns all roots of ratfun.
%
% rts = roots(ratfun, 'real') 
% Returns all real roots of ratfun.

  method = 0; 
  realonly = 0;
  for j = 1:length(varargin)
    if isnumeric(varargin{j}),      method = varargin{j}; end
    if strcmp(varargin{j}, 'real'), realonly = 1;         end
  end
  
  K = obj.K;
  H = obj.H;
  m = size(H, 2);
  coeffs = obj.coeffs;
  nrm = norm(coeffs);
  
  if nrm == 0, [t1,t2] = type(obj); rts = zeros(t1+1, 1); return, end

  coeffs = coeffs/nrm;  
  
  switch method

   case 0 
    %%% Using the trailing m-by-m part. %%%
    [u, sigma] = util_householder(coeffs);    
    if norm(u) == 0, rts = []; return, end
    K = K - sigma*u*(u'*K);
    H = H - sigma*u*(u'*H);
    if isa(H, 'sym'), rts = eig(H(2:m+1, :)/K(2:m+1, :)); 
    else              rts = eig(H(2:m+1, :),K(2:m+1, :)); end    
    % Only output the k smallest roots, where k is the numerator
    % degree.          
    t = type(obj);
    [~, ind] = sort(abs(rts), 'ascend'); 
    rts = rts(ind);
    rts = rts(1:t(1));    

   case 1 
    %%% Using leading m-by-m part. %%%
    
    % We define the m-by-m identity in this way for for vpa/mp.    
    M = eye(m+1, m+1)^0;
    M(:, end) = coeffs;
    K = M\K;
    H = M\H; 
    if isa(H, 'sym'), rts = eig(H(1:m,:)/K(1:m,:)); 
    else              rts = eig(H(1:m,:),K(1:m,:)); end    
    % Only output the k smallest roots, where k is the numerator
    % degree.      
    t = type(obj);
    [~, ind] = sort(abs(rts), 'ascend');
    rts = rts(ind);
    rts = rts(1:t(1));
    
   case 2 
    %%% Transforming to polynomial RAD first. %%%
    
    [K, H, Q] = util_hh2th(K, H);
    t = type(obj); 
    t = t(1);    
    coeffs = Q*coeffs;
    coeffs = coeffs(1:t+1);
    K = K(1:t+1, 1:t); % K = triu(K(1:t+1, 1:t));
    H = H(1:t+1, 1:t); % H = triu(H(1:t+1, 1:t), -1);    
    % Now compute the roots.
    rts = roots(rkfun(K, H, coeffs), 0);
            
  end % switch method
  
  %{
  % Newton correction.
  iter = 1; 
  val  = feval(obj, rts); 
  N    = 2*length(rts);
  while iter <= 5 && norm(val, inf) > 1e-16
    lam([1:2:N 2:2:N], 1) = [rts;rts];
    b = []; b(N, 1) = 0; b(2:2:end) = 1;
    A = spdiags([lam, b], 0:1, N, N);
    
    der = feval(obj, A, b); 
    der = der(1:2:end);
    val = feval(obj, rts); 
    rts = rts - val./der;
    iter = iter + 1;
  end
  %}
  
  % Filter real roots if required.
  if realonly
    ind = find(abs(imag(rts)) < 1e-15);
    rts = sort(real(rts(ind)));
  end
  
end % function

