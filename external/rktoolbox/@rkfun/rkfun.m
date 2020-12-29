%RKFUN    RKFUN constructor.
%
% Usages:
% 1) obj = rkfun(K, H, coeffs, k) constructs an RKFUN from the
%    recursion matrices (H,K), and the coefficient vector coeffs.
%    The optional parameter k specifies sub- or superdiagonal
%    rational functions.
% 2) obj = rkfun(str, params) returns an RKFUN from the gallery. 
%    To see a list of the gallery functions type help RKFUN/gallery.
%    Example: r = rkfun('step'); ezplot(r, [-2, 2])
% 3) obj = rkfun(str) constructs an RKFUN by evaluating the string 
%    str symbolically. Example: rkfun('(2*x-3)/(x+4)')


classdef rkfun
properties
K
H
coeffs
k
iscolumn
prefact
end % properties
    
methods(Access = public, Static = true)
% Constructor. (in methods)
function obj = rkfun(varargin)
  
  % rkfun representing r(x) = x.
  if nargin == 0
    obj = nodes2rkfun(0);
    return
  end
  
  % Constant rkfun.
  if nargin == 1    
    if isnumeric(varargin{1}) && isscalar(varargin{1})
      obj = nodes2rkfun([] ,[], varargin{1});
      return
    end
  end
  
  % Gallery.
  if ischar(varargin{1})
    obj = rkfun.gallery(varargin{:});
    return
  end
  
  % Construct rkfun from pencil and coeffs.
  obj.K = varargin{1};
  obj.H = varargin{2};
  obj.coeffs = varargin{3};
  if nargin >= 4
    obj.k = varargin{4};
  else
    obj.k = 0;
  end
  if nargin >= 5
    obj.iscolumn = varargin{5};
  else
    obj.iscolumn = 1;
  end
  
end % function rkfun
        
function obj = gallery(varargin)
% GALLERY    Collection of rational functions.
%
% obj = rkfun.gallery(funname, param1, param2, ...) takes 
% funname, a case-insensitive string that is the name of 
% a rational function family, and the family's input 
% parameters. 
%
% See the listing below for available function families.
%
% constant   Constant function of value param1.
% cheby      Chebyshev polynomial (first kind) of degree param1.
% cayley     Cayley transformation (1-z)/(1+z).
% moebius    Moebius transformation (az+b)/(cz+d) with 
%            param1 = [a,b,c,d].
% sqrt       Zolotarev sqrt approximation of degree param1 on 
%            the positive interval [1,param2].
% invsqrt    Zolotarev invsqrt approximation of degree param1 on 
%            the positive interval [1,param2].
% sqrt0h     balanced Remez approximation to sqrt(x+(h*x/2)^2)
%            of degree param3 on [param1,param2], 
%            where param1 <= 0 <= param2 and h = param4.
% sqrt2h     balanced Zolotarev approximation to sqrt(x+(h*x/2)^2)
%            of degree param5 on [param1,param2]U[param3,param4], 
%            where param1 < param2 < 0 < param3 < param4 and h = param6.
% invsqrt2h  balanced Zolotarev approximation to 1/sqrt(x+(h*x/2)^2)
%            of degree param5 on [param1,param2]U[param3,param4], 
%            where param1 < param2 < 0 < param3 < param4 and h = param6.
% sign       Zolotarev sign approximation of degree 2*param1 on
%            the union of [1,param2] and [-param2,-1].
% step       Unit step function approximation for [-1,1] of 
%            degree 2*param1 with steepness param2.  

  obj = gallery(varargin{:});

end
                
function obj = nodes2rkfun(varargin)
% NODES2RKFUN    Create an RKFUN by providing roots and poles.
%
% obj = rkfun.nodes2rkfun(rts, pls, scl) constructs an RKFUN with
% roots rts and poles pls. If no scaling paramter scl is provided,
% the RKFUN is such that both the numerator and denominator are
% monic polynomials. Otherwise the scalar multiple of this monic
% function is returned.    

  obj = nodes2rkfun(varargin{:});
  
end

end % methods    
end % classdef rkfun
