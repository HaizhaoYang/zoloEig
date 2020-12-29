function obj = gallery(varargin)

  obj = 0;
  switch lower(varargin{1})
   case 'constant'
    p = varargin{2};
    obj = rkfun(zeros(1, 0), zeros(1, 0), p(1));     
     
   case 'cheby'
    if nargin < 2, n = 5; else n = varargin{2}; end
    e = .5*ones(n,1); e(1) = 1;     
    H = full(spdiags([e,0*e,e],-1:1,n+1,n));    
    K = eye(n+1,n); 
    c = zeros(n+1,1); 
    c(n+1) = 1;
    obj = rkfun(K, H, c, n);
    
   case 'cayley'
    obj = nodes2rkfun(1,-1,-1);
    
   case 'moebius'
    p = varargin{2};
    obj = rkfun([-p(1) ; p(3)], [p(2) ; -p(4)], [0;1]);

   case 'sqrt'
    % Zolotarev's sqrt approximation on [1, b] of type (d, d-1).         
    if nargin < 2, d = 8;   else d = varargin{2}; end
    if nargin < 3, b = 1e2; else b = varargin{3}; end
    d  = max(d, 1); 
    b  = max(b, 1);
    obj = 1./rkfun('invsqrt', d, b);
        
   case 'invsqrt'
    % Zolotarev's invsqrt approximation on [1, b] of type (d-1, d).         
    % using quotient form
    if nargin < 2, d = 8;   else d = varargin{2}; end
    if nargin < 3, b = 1e2; else b = varargin{3}; end
    d  = max(d, 1); 
    b  = max(b, 1);
    
    % use MATLAB's elliptic functions:
    %K = ellipke(1-1/b);
    %[sn, cn, dn] = ellipj((0:2*d)*K/(2*d), 1-1/b);
    
    % alternatively use functions from Driscoll's SC Toolbox
    L = .5*log(b)/pi; 
    [~, K] = util_ellipkkp(L);
    [sn, cn, dn] = util_ellipjc(1i*(0:2*d)*K/(2*d), L);
    cn = 1./cn; dn = dn.*cn; sn = -1i*sn.*cn;
        
    c = sn.^2./cn.^2;
    zero = -c(3:2:end-1); 
    pole = -c(2:2:end);
    extr = dn.^-2;
    
    % Determine scaling factor D.
    r = 1;
    for j = 1:d,
      if j < d,
        r = r.*(extr-zero(j))./(extr-pole(j));
      else
        r = r./(extr-pole(j));
      end
    end
    r = r.*sqrt(extr);
    D = 2/(min(r) + max(r));
    obj = nodes2rkfun(zero, pole, D);
    
   case 'invsqrt_pf'
    % Zolotarev's invsqrt approximation on [1, b] of type (d-1, d).         
    % using partial fraction form
    if nargin < 2, d = 8;   else d = varargin{2}; end
    if nargin < 3, b = 1e2; else b = varargin{3}; end
    d  = max(d, 1); 
    b  = max(b, 1);
    k2 = 1/b; % elliptic functions parameter k^2
    Kp = ellipke(1-k2);
    t = 1i*(.5:d)*Kp/d;
    [sn,cn,dn] = ellipj(imag(t),1-k2);
    cn = 1./cn; dn = dn.*cn; 
    pls = (1i*sn.*cn).^2;
    res = 2*Kp*cn.*dn/(pi*d);
    % get scaling right by evaluating error at extrema
    [~,~,ddn] = ellipj((0:.5:d)*Kp/d,1-k2);
    ext = 1./ddn.^2;
    s = 0*ext; 
    for j = 1:d,
        s = s + res(j)./(ext-pls(j));
    end
    minS = min(s.*sqrt(ext));
    maxS = max(s.*sqrt(ext));
    alph = 2/(minS + maxS);
    res = alph*res;
    % build rkfun
    obj = 0;
    for j = 1:d,
        obj = obj + rkfun('moebius',[0,res(j),1,-pls(j)]);
    end
    obj.k = -1;

   case 'invsqrt2h'
    % balanced Zolotarev approximant for 1/sqrt(x+(h*x/2)^2) 
    % of degree (m-1,m) on [a1,b1]U[a2,b2], where a1<b1<0<a2<b2.
    % The construction is described in 
    % V. Druskin, D. Guettel, L. Knizhnerman, Near-optimal perfectly 
    % matched layers for indefinite Helmholtz problems, SIAM Review, 2016.
    a1 = varargin{2}; b1 = varargin{3}; a2 = varargin{4}; b2 = varargin{5};
    if nargin < 6, m = 8;   else m = varargin{6}; end
    if nargin < 7, h = 0;   else h = varargin{7}; end
    
    % balance approximation degrees on [a1,b1] with [a2,b2]
    delta1 = sqrt(b1/a1); mu1 = ((1-sqrt(delta1))/(1+sqrt(delta1)))^2; 
    mu1p = sqrt(1-mu1^2); rho1 = exp(-pi/4*ellipke(mu1p^2)/ellipke(mu1^2));
    delta2 = sqrt(a2/b2); mu2 = ((1-sqrt(delta2))/(1+sqrt(delta2)))^2; 
    mu2p = sqrt(1-mu2^2); rho2 = exp(-pi/4*ellipke(mu2p^2)/ellipke(mu2^2));
    m1 = round(2*m*log(rho2)/(log(rho1) + log(rho2)));
    m2 = round(2*m*log(rho1)/(log(rho1) + log(rho2)));
    
    if (m1 + m2)/2 ~= m,
        error('degree mismatch');
    end
    
    % compute interpolation points
    j = (2*m1 - 2*(1:m1) + 1)/(2*m1);
    L = .5*log(a1/b1)/pi; 
    [~, K] = util_ellipkkp(L);
    [sn, cn, dn] = util_ellipjc(1i*j*K, L);
    cn = 1./cn; dn = dn.*cn; 
    rts1 = a1*dn.^2;
    
    j = (2*m2 - 2*(1:m2) + 1)/(2*m2);
    L = .5*log(b2/a2)/pi; 
    [~, K] = util_ellipkkp(L);
    [sn, cn, dn] = util_ellipjc(1i*j*K, L);
    cn = 1./cn; dn = dn.*cn; 
    rts2 = b2*dn.^2;
    
    % compute 'interpolant' by LS fit of exact type
    x = [rts1,rts2];
    A = diag(x);
    f = 1./sqrt(x+(h*x/2).^2);
    F = diag(f);
    b = ones(size(A,1),1);
    xi = inf(1, m);
    [xi, obj, misfit] = rkfit(F, A, b, xi, ...
                                    struct('tol', 5e-15, ...
                                    'maxit', 10, ...
                                    'k', -1, ...
                                    'reduction', 0));
                       
   case 'sqrt0h'
    % balanced Remez approximant for sqrt(x+(h*x/2)^2) 
    % of degree (m,m-1) on [a1,b2], where a1<=0<=b2.
    a1 = varargin{2}; b2 = varargin{3};
    if nargin < 4, m = 8;   else m = varargin{4}; end
    if nargin < 5, h = 0;   else h = varargin{5}; end
    obj = ipsqrt(a1,b2,m,h);       
    
   case 'sqrt2h'
    % balanced Zolotarev approximant for sqrt(x+(h*x/2)^2) 
    % of degree (m,m-1) on [a1,b1]U[a2,b2], where a1<b1<0<a2<b2.
    % The construction is described in 
    % V. Druskin, D. Guettel, L. Knizhnerman, Near-optimal perfectly 
    % matched layers for indefinite Helmholtz problems, SIAM Review, 2016.
    a1 = varargin{2}; b1 = varargin{3}; a2 = varargin{4}; b2 = varargin{5};
    if nargin < 6, m = 8;   else m = varargin{6}; end
    if nargin < 7, h = 0;   else h = varargin{7}; end
    
    % balance approximation degrees on [a1,b1] with [a2,b2]
    delta1 = sqrt(b1/a1); mu1 = ((1-sqrt(delta1))/(1+sqrt(delta1)))^2; 
    mu1p = sqrt(1-mu1^2); rho1 = exp(-pi/4*ellipke(mu1p^2)/ellipke(mu1^2));
    delta2 = sqrt(a2/b2); mu2 = ((1-sqrt(delta2))/(1+sqrt(delta2)))^2; 
    mu2p = sqrt(1-mu2^2); rho2 = exp(-pi/4*ellipke(mu2p^2)/ellipke(mu2^2));
    m1 = round(2*m*log(rho2)/(log(rho1) + log(rho2)));
    m2 = round(2*m*log(rho1)/(log(rho1) + log(rho2)));
    
    if (m1 + m2)/2 ~= m,
        error('degree mismatch');
    end
    
    % compute interpolation points
    j = (2*m1 - 2*(1:m1) + 1)/(2*m1);
    L = .5*log(a1/b1)/pi; 
    [~, K] = util_ellipkkp(L);
    [sn, cn, dn] = util_ellipjc(1i*j*K, L);
    cn = 1./cn; dn = dn.*cn; 
    rts1 = a1*dn.^2;
    
    j = (2*m2 - 2*(1:m2) + 1)/(2*m2);
    L = .5*log(b2/a2)/pi; 
    [~, K] = util_ellipkkp(L);
    [sn, cn, dn] = util_ellipjc(1i*j*K, L);
    cn = 1./cn; dn = dn.*cn; 
    rts2 = b2*dn.^2;
    
    % compute 'interpolant' by LS fit of exact type
    x = [rts1,rts2];
    A = diag(x);
    f = sqrt(x+(h*x/2).^2);
    F = diag(f);
    b = ones(size(A,1),1);
    xi = inf(1, m-1);
    [xi, obj, misfit] = rkfit(F, A, b, xi, ...
                                    struct('tol', 5e-15, ...
                                    'maxit', 10, ...
                                    'k', 1, ...
                                    'reduction', 0));
                                
   case 'sign'
    % Zolotarev's sign approximation on +/-[1, b] 
    % of type (2*d-1, 2*d).    
    if nargin < 2, d = 8;   else d = varargin{2}; end
    if nargin < 3, b = 1e4; else b = varargin{3}; end
    d  = max(d, 1); 
    b  = max(b, 1);

    % use MATLAB's elliptic functions
    %K = ellipke(1-1/b^2);
    %[sn, cn, dn] = ellipj((0:2*d)*K/(2*d), 1-1/b^2);
    
    % alternatively use functions from Driscoll's SC Toolbox
    L = log(b)/pi; 
    [~, K] = util_ellipkkp(L);
    [sn, cn, dn] = util_ellipjc(1i*(0:2*d)*K/(2*d), L);
    cn = 1./cn; dn = dn.*cn; sn = -1i*sn.*cn;
    
    c = sn.^2./cn.^2;
    zero = -c(3:2:end-1); 
    pole = -c(2:2:end);
    extr = dn.^-2;
    
    % Determine scaling factor D.
    r = 1;
    for j = 1:d,
      if j < d,
        r = r.*(extr-zero(j))./(extr-pole(j));
      else
        r = r./(extr-pole(j));
      end
    end
    r = r.*sqrt(extr);
    D = 2/(min(r) + max(r));
    
    % now D*r is relative approximation to 1/sqrt(x) on [0,b^2]
    % now use sign(x) = x/sqrt(x^2)
    rts(1:2:2*(d-1)) = +sqrt(zero);
    rts(2:2:2*(d-1)) = -sqrt(zero);
    pls(1:2:2*d) = +sqrt(pole);
    pls(2:2:2*d) = -sqrt(pole);
    rts(2*d-1) = 0;
    obj = nodes2rkfun(rts, pls, D);
    
    case 'step'
     if nargin < 2, d = 8;   else d = varargin{2}; end
     if nargin < 3, b = 1e4; else b = varargin{3}; end
     d = max(d, 1); 
     b = max(b, 1);
     % Zolotarev's sign approximation on [1, b] of 
     % type (2*d-1, 2*d).
     r = gallery('sign', d, b);
     % Compose with Moebius transform.
     m = -sqrt(b)*gallery('moebius',[1, 1, 1, -1]); 
     obj = 0.5*feval(r, m) + 0.5;
			
   otherwise 
    obj = sym2rkfun(varargin{1});    
  end
        
end % function gallery
