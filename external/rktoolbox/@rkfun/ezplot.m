function varargout = ezplot(obj, varargin)
%EZPLOT   Easy-to-use function plotter. 
%
% Usages: - ezplot(rkfun)
%         - ezplot(rkfun, interval)
%         - ezplot(rkfun, opts)
%         - ezplot(rkfun, interval, opts)
% This function plots the real part of an rkfun over an x-range 
% determined by interval, or, if interval is not given, 
% by Ritz values associated with the pencil (K,H). The function
% accepts the plot parameters of the standard PLOT function and
% can return a handle to the plot if required.
% Example: hdl = ezplot(rkfun('cheby',3),[-1.2,1.2],'r','LineWidth',2)

  obj = double(obj);
  interval = []; opts = [];
  if nargin >= 2,
     if isnumeric(varargin{1}),
         interval = varargin{1};
         if nargin >= 3, opts = varargin(2:end); end
     else
         opts = varargin;
     end
  end
  
  if isempty(interval),
    ee = eig(obj.H(1:end-1,1:end),obj.K(1:end-1,1:end));
    x1  = min(real(ee(isfinite(ee))));
    x2  = max(real(ee(isfinite(ee))));
    d   = max(x2-x1,1);
    x1p = x1 - 0.2*d;
    x2p = x2 + 0.2*d;
  else
    x1 = interval(1); x1p = x1;
    x2 = interval(2); x2p = x2;
  end
    
  xx = linspace(x1p, x2p, 813);
  yy = double(feval(obj, xx));
  if isempty(opts), 
      hdl = plot(xx,yy); 
  else
      hdl = plot(xx,yy,opts{:}); 
  end
  xlim([x1p,x2p])
  yy = real(yy(xx >= x1p & xx <= x2p));
  y1 = min(yy);
  y2 = max(yy);
  d  = max(y2-y1,1e-12);
  ylim([y1 - 0.1*d,y2 + 0.1*d]);
  
  if nargout >= 1,
      varargout = {hdl}; 
  else
      varargout = {};
  end
  
end

