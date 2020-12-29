n = 50;
G = numgrid('D', n); 
A = (n+1)^2*delsq(G); 
N = size(A, 1);
B = speye(N);
[V, D]   = eigs(A, 11, 'SM');
[D, ind] = sort(diag(D)); V = V(:, ind);

figure(1)
for j = 1:9
  subplot(3, 3, j)
  v = V(:, j);
  U = NaN*G;
  U(G>0) = full(v(G(G>0)));
  contourf(U);
  prism, axis square ij, axis off
  title(['\lambda = ' num2str(D(j))]);
end

%%
fprintf('test sparse zoloeigs\n');
opt = [];
opt.reltol = 1e-13;
opt.verbose = true;
opt.nc = 8;
eigvals = zologeigs(A,B,0,200,opt);

norm(sort(eigvals) - D(1:8))/norm(D(1:8))

%%
[Uoutf,eigvalsf,relerrsf] = rkfeast(A,0,200,8,16,opt);
norm(sort(eigvalsf) - D(1:8))/norm(D(1:8))

%%
fprintf('test dense/function handle zoloeigs\n');
opt = [];
opt.reltol = 1e-13;
opt.verbose = true;
opt.nc = 8;
n = size(A,1);
Afun = @(x) A*x; 
Bfun = @(x) B*x;
GenAinvfun = @(sigma) Ainvfunc(A,B,sigma);
a = 0; b = 200; 
eigvals = zologeigsdense(Afun,GenAinvfun,Bfun,n,a,b,opt);
norm(sort(eigvals) - D(1:8))/norm(D(1:8))

