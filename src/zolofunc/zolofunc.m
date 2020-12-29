function func = zolofunc(l,r,scale)

if nargin == 1
    r = zolopara(l);
end

if nargin < 3
    scale = 0;
end

if length(r) <= 1

    [c,a] = zolocoef(l,r,scale);
    
    func = @(x) 0;
    for it = 1:r
        func = @(x) func(x) + (a(it)*x)./(x.^2+c(2*it-1));
    end
    
else
    func0 = zolofunc(l,r(1:end-1),1);
    
    lnext = func0(l);
    
    [c,a] = zolocoef(lnext,r(end),scale);
    func = @(x) 0;
    for it = 1:r(end)
        func = @(x) func(x) + (a(it)*func0(x))./(func0(x).^2+c(2*it-1));
    end
    
end

end