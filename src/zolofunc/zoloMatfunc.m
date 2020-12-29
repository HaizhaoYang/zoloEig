function func = zoloMatfunc(l,r,scale)

if nargin == 1
    r = zolopara(l);
end

if nargin < 3
    scale = 0;
end

if length(r) <= 1

    [c,a] = zolocoef(l,r,scale);
    
    func = @(A) 0;
    for it = 1:r
        func = @(A) func(A) + (A*A+c(2*it-1)*eye(size(A)))\(a(it)*A);
    end
    
else
    
    func0 = zoloMatfunc(l,r(1:end-1),1);
    
    lnext = func0(l);
    
    [c,a] = zolocoef(lnext,r(end),scale);
    func = @(A) 0;
    for it = 1:r(end)
        func = @(A) func(A) ...
            + (func0(A)*func0(A)+c(2*it-1)*eye(size(A)))\(a(it)*func0(A));
    end
    
end

end