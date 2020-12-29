n = 1000;

Q = randn(n);
[Q,~] = qr(Q);

ev = linspace(0.13,10,n)';
A = Q*diag(ev)*Q';

a = 4;
nc = 10;
b = ev(find(ev>a,1)+nc)*0.4+0.6*ev(find(ev>a,1)+nc-1);

global ev_ref
ev_ref = ev(ev>a & ev<b);

tic;
opt = [];
opt.reltol = 1e-13;
opt.verbose = 1;
ev2 = zoloeigs(A,a,b,opt);
toc;

err = norm(sort(ev2) - ev_ref)/norm(ev_ref)