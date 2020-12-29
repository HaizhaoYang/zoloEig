function [c,a] = zolocoef(l,r,scale)

if nargin < 3
    scale = 1;
end

c = nan(1,2*r-1);
a = nan(1,r);
b = nan(1,r-1);

kp = l;
alpha = acos(kp);
K = mellipke(alpha); % modified elliptic functions

for it = 1:2*r-1
    [sn,cn,~] = mellipj(it*K/(2*r),alpha);
    c(it) = (sn/cn*kp)^2;
end

for iti=1:r-1
    enu = 1;
    for itj = 1:r-1
        enu = enu*(c(2*iti-1)-c(2*itj));
    end
    den = 1;
    for itj = 1:r-1
        if iti~=itj
            den = den*(c(2*iti-1)-c(2*itj-1));
        end
    end
    b(iti) = -enu/den;
end

a(1:r-1) = b./(c(2*r-1)-c(2*(1:r-1)-1));
a(r) = 1-sum(a(1:r-1));

dn = nan(1,2);
[~,~,dn(1)] = mellipj(1*K/(2*r),alpha);
[~,~,dn(2)] = mellipj(2*K/(2*r),alpha);
ex = dn;

ey = ex./(ex.^2+c(2*r-1));
for itj = 1:r-1
    ey = ey.*(ex.^2+c(2*itj))./(ex.^2+c(2*itj-1));
end

if scale
    a = a./max(ey);
else
    a = a./mean(ey);
end

end