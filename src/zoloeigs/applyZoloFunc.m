function [Y,niter] = applyZoloFunc(R,c,a,X,itsolacc)

r = length(a);
Y = zeros(size(X));

if isreal(X)

    [SRX,~,~,niter] = gmres_shift(R,X,sqrt(c(2*(1:r)-1))*1i,100,itsolacc,100);
    SRX = real(SRX);
    for it=1:r
        Y = Y + a(it)*SRX(:,:,it);
    end
    Y = real(Y);
    
else
    
    [SRX,~,~,niter] = gmres_shift(R,X, ...
        [sqrt(c(2*(1:r)-1))*1i -sqrt(c(2*(1:r)-1))*1i],100,itsolacc,100);
    for it=1:r
        Y = Y + a(it)/2*SRX(:,:,it) + a(it)/2*SRX(:,:,it+r);
    end
    
end

end