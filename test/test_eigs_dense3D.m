maxNumCompThreads(1);

run('../EigSol_startup.m');

siz = 8:4:64;

nmat = length(siz);

trail = 10;

nc = 10;

timeigs = zeros(trail,nmat);
timzolo = zeros(trail,nmat);
errzolo = zeros(trail,nmat);
gcond = zeros(1,nmat);

for it = 1:nmat
    
    n = siz(it);
    
    fprintf('=====================================================\n');
    fprintf('Matrix is of size %4d by %4d.\n', ...
        n^3, n^3 );
    fprintf('-----------------------------------------------------\n');
    
    boxsize = 1;
    
    v = myinterpft2D(rand(8, 8), n);
    vmin = min(v(:));
    vmax = max(v(:));
    v = (pi * n / 2) ^ 2 * (2 / (vmax - vmin) * (v - vmin) - 1);
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
    
    for ittrail = 1:trail
        
        opteigs = [];
        opteigs.isreal = true;
        opteigs.issym = true;
        tic;
        ev = eigs(Afun,N,nc+1,'SA',opteigs);
        timeigs(ittrail,it) = toc;
        ev_ref = ev(1:nc);
        a = [ev(1)-10 ev(1)];
        b = [ev(nc) ev(nc+1)];
        
        tic;
        opt = [];
        opt.verbose = 1;
        opt.nc = nc;
        opt.reltol = 1e-5;
        ev2 = zoloeigsdense(Afun,Ainvfun,N,a,b,opt);
        timzolo(ittrail,it) = toc;
        
        errzolo(ittrail,it) = norm(sort(ev2) - ev_ref)/norm(ev_ref);
        
        fprintf('Zolo eig finishes in %.2e sec.\n', timzolo(ittrail,it) );
        fprintf('The rel err is %.2e.\n', errzolo(ittrail,it) );
        
        fprintf('-----------------------------------------------------\n');
    end
end

save('ZoloDense3D');

%%

siz = siz.^3;
figure(1)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
errorbar(siz,mean(timeigs), ...
    mean(timeigs)-min(timeigs),max(timeigs)-mean(timeigs), ...
    'Linewidth',2);
errorbar(siz,mean(timzolo), ...
    mean(timzolo)-min(timzolo),max(timzolo)-mean(timzolo), ...
    'Linewidth',2);
plot(siz,min(min(timzolo))/min(siz)*siz,':','Linewidth',2);
plot(siz,min(min(timzolo))/min(siz)^2*siz.^2,':','Linewidth',2);
xlabel('Size');
ylabel('Time');
legend('Eigs time','ZoloEig time','linear ref','quadratic ref');
saveas(gcf,'TimZoloDense3D.fig','fig');

figure(2)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
errorbar(siz,mean(errzolo), ...
    mean(errzolo)-min(errzolo),max(errzolo)-mean(errzolo), ...
    'Linewidth',2);
xlabel('Size');
ylabel('Relative Error');
legend('ZoloEig RelErr');
saveas(gcf,'ErrZoloDense3D.fig','fig');
