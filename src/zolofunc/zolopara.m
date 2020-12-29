function r = zolopara(l,tol)
% choose the Zolotarev degree according to the size of the eigengap
if nargin < 2
    tol = 5e-15;
end

for r1 = 1:3
    r = r1;
    if 1-zoloeval(l,l,r) < tol/10
        return;
    end
end

for r1 = 2:8
    for r2 = r1:r1% TODO: There is a bug when r1 ~= r2
        r = [r1 r2];
        l1 = zoloeval(l,l,r1,1);
        if 1-zoloeval(l1,l1,r2,0) < tol/10
            return;
        end
    end
end

r = [8 8];
warning('The gap is too small, r=[8,8] is used anyway');

end