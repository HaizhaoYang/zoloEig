function y = applyA2D(x, v, h, ecut)

if nargin < 4
    ecut = Inf;
end

n = size(v,1);
kk=[(0:(n+mod(n,2))/2-1),(-(n-mod(n,2))/2:-1)];
kk=(2*pi*kk).^2;
kk = bsxfun(@plus,kk,kk');
kk(kk>ecut)=ecut;

xten = reshape(x, n, n, []);

xten = ifft2(bsxfun(@times,fft2(xten),kk));

y = reshape( (1/(n*h)^2)*xten + ...
    bsxfun( @times, v, reshape(x,n,n,[]) ), n^2, [] );

end