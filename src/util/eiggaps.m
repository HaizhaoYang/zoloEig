function varargout = eiggaps(varargin)
% find buffer range [a1,a2] and [b1,b2], where a1<a2<b1<b2

A = varargin{1};

if isa(A,'function_handle')
    B = varargin{2};
    if isa(B,'function_handle')
        n = varargin{3};
        offset = 3;
        opts = struct('Tol',1e-13,'MaxIt',5000); % TODO: add preconditioner
        eiggapfun = @(a)real(jdqzGap(n,A,B,a,opts));
    else
        n = varargin{2};
        offset = 2;
        opts.issym = false;
        opts.isreal = false;
        eiggapfun = @(a)real(eigs(A,n,1,a,opts));
    end
else
    B = varargin{2};
    if all(size(A) == size(B))
        offset = 2;
        eiggapfun = @(a)real(eigs(A,B,1,a));
    else
        offset = 1;
        eiggapfun = @(a)eigs(A,1,a);
    end
end

varlen = nargin - offset;
varargout = cell(1,2*varlen);

emax = eiggapfun('LR');
emin = eiggapfun('SR');

for it = 1:varlen
    
    a = varargin{offset+it};
    
    if a < emin
        
        em = a - (emax-emin);
        ep = a;
        
    elseif a > emax
        
        em = a;
        ep = a + (emax-emin);
        
    else
        
        ep = [];
        em = [];
        
        while isempty(em) || isempty(ep)
            eval = eiggapfun(a);
            if eval <= a
                em = eval;
            end
            if eval >= a
                ep = eval;
            end
            a = 2*a-eval;
        end
        
    end
    
    varargout{2*it-1} = em;
    varargout{2*it  } = ep;
end

end

function a = jdqzGap(n,Afun,Bfun,SIGMA,opts)
[~,a] = jdqzFun(n,Afun,Bfun,1,SIGMA,opts);
end