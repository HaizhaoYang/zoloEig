function [Y,niter] = applyZoloDenseMat(cc,Ainvzolo,X,Bfun)
if nargin <4
    BX = X;
else
    BX = Bfun(X);
end

niter = 0;
r = length(Ainvzolo);

Y = cc*X;
for it=1:r
    Ainvfun = Ainvzolo{it};
    Y = Y + 2*real(Ainvfun(BX));
end

end