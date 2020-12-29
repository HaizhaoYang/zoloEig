n = 8 * 4;
boxsize = 8;

v = myinterpft2D(rand(8, 8), n);
vmin = min(v(:));
vmax = max(v(:));
v = (2 * pi * n / 16) ^ 2 * (2 / (vmax - vmin) * (v - vmin) - 1);

nlocal = n;

h = 1 / n;
allvmode = 2;
reach = 2;
spread = 2;

ecut=Inf;
ecut=0.5*2*(2*pi*n/2)^2;

N = n ^ 2;
Afun=@(u)applyA2D(u,v,h,ecut);
Ainvfun=@(sigma)genAinvFunLocal2D(v, reach, spread, allvmode, boxsize, h, n, nlocal, sigma, ecut);

Ainv0=Ainvfun(0);
f=randn(N,1);
u=Ainv0(f);