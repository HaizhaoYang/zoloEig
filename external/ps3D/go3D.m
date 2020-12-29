n = 4 * 4;
boxsize = 8;

v = myinterpft3D(rand(4, 4, 4), n);
vmin = min(v(:));
vmax = max(v(:));
v = (pi * n / 2) ^ 2 * (2 / (vmax - vmin) * (v - vmin) - 1);

nlocal = n;

h = 1 / n;
allvmode = 1;
reach = 2;
spread = 2;

N = n ^ 3;
Afun=@(u)applyA3D(u,v,h);
Ainvfun=@(sigma)genAinvFunLocal3D(v, reach, spread, allvmode, boxsize, h, n, nlocal, sigma);