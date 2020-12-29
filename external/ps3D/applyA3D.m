function f = applyA3D(u,v,h)

n=size(v,1);
kk=[(0:(n+mod(n,2))/2-1),(-(n-mod(n,2))/2:-1)];
kk=(2*pi*kk).^2;
[kk1,kk2,kk3]=ndgrid(kk);
kk=kk1+kk2+kk3;

u=reshape(u,n,n,n,[]);
f=(1/(n*h)^2)*ifft3(bsxfun(@times,fft3(u),kk))+bsxfun(@times,u,v);
f=reshape(f,n^3,[]);

end