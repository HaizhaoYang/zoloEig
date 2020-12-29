function [Y,niter] = applyZoloMat(A,B,cc,prefact,X)

niter = 0;
r = length(prefact);

if isreal(A) && isreal(B)
    BX = B*X;
    Y = cc*X;
    for it=1:r
        MF = prefact{it};
        Y = Y + 2*real(MF\BX);
    end
else
    BX = B*X;
    Y = cc*X;
    for it=1:r
        MF = prefact{it};
        Y = Y + MF\BX + (BX'/MF)';
    end
end

end