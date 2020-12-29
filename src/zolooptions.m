function [nc,r,reltol,itmax,verbose,itsolacc] = ...
    zolooptions(varargin)

nc       = 200;
r        = [];
reltol   = 1e-8;
itmax = 10;
verbose  = 0;
itsolacc = min(1e-13,reltol/100);

if nargin == 1 && ~isnumeric(varargin{1})
    opt = varargin{1};
    if isfield(opt,'nc')
        nc = opt.nc;
    end
    if isfield(opt,'r')
        r = opt.r;
    end
    if isfield(opt,'reltol')
        reltol = opt.reltol;
    end
    if isfield(opt,'itmax')
        itmax = opt.itmax;
    end
    if isfield(opt,'verbose')
        verbose = opt.verbose;
    end
    if isfield(opt,'itsolacc')
        itsolacc = opt.itsolacc;
    end
end

if nargin >= 1 && isnumeric(varargin{1})
    if nargin >= 1
        nc = varargin{1};
    end
    if nargin >= 2
        r = varargin{2};
    end
    if nargin >= 3
        reltol = varargin{3};
    end
    if nargin >= 4
        itmax = varargin{4};
    end
    if nargin >= 5
        verbose = varargin{5};
    end
    if nargin >= 6
        itsolacc = varargin{6};
    end
end

end