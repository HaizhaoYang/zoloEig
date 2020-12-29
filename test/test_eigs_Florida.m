maxNumCompThreads(1);

UF = UFget;           % get index of the UF Sparse Matrix Collection
siz = UF.nrows;
nz = UF.nnz;

nmat = length(UF.Name);

nc = 100;
opt.itmax = 1;
opt.verbose = 0;

if exist('Eigs_Florida_Data.mat','file')
    load Eigs_Florida_Data;
else
    sit = 1;
    timeig = zeros(1,nmat);
    timfeast = zeros(1,nmat);
    timzolo = zeros(1,nmat);
    errfeast = zeros(1,nmat);
    errzolo = zeros(1,nmat);
    idxeig = false(1,nmat);
    gcond = zeros(1,nmat);
end

for it = sit:nmat
    
    if UF.numerical_symmetry(it)~=1 || ~UF.isReal(it)
        continue;
    end
    
    if siz(it) > 6000 || siz(it) < 200
        continue;
    end
    
    sit = it;
    save('Eigs_Florida_Data','sit','timeig','timfeast','timzolo', ...
        'errfeast','errzolo','idxeig','gcond');
    
    Prob = UFget(it);
    
    A = Prob.A;
    
    fprintf('=========================================================\n');
    fprintf('Matrix name: %s/%s\n',UF.Group{it},UF.Name{it});
    fprintf('Matrix index: %3d out of %3d\n',it,nmat);
    fprintf('Matrix is of size %4d with %8d non-zeros.\n', ...
        siz(it), nz(it) );
    fprintf('---------------------------------------------------------\n');
    
    tic;
    ev = eig(A);
    timeig(it) = toc;
    ev = sort(ev);
    
    a0 = min(ev);
    b0 = max(ev);
    
    fprintf('Matlab eig finishes in %.2e sec.\n', timeig(it) );
    
    ittry = 1;
    gcond(it) = 1e15;
    while ittry < 100 && ( gcond(it) > 1e8 || isnan(gcond(it)) )
        mid = randi(siz(it)-nc-1);
        gcond(it) = abs(b0-a0) ...
            /min(abs(ev(mid+1)-ev(mid)),abs(ev(mid+nc+1)-ev(mid+nc)));
        ittry = ittry+1;
    end
    
    if ittry >= 100
        continue;
    end
    
    fprintf('The general condition number is %.2e.\n', gcond(it) );
    
    fprintf('---------------------------------------------------------\n');
    
    idxeig(it) = true;
    a = ev(mid)*0.4+0.6*ev(mid+1);
    b = ev(mid+nc)*0.4+0.6*ev(mid+nc+1);
    
    ev_ref = ev(ev>a & ev<b);
    
    %--------------------------------------------------------------------
    % FEAST eig
    tic;
    opt = [];
    opt.itmax = 1;
    evslic = rkfeast(A,a,b,nc,32,opt);
    timfeast(it) = toc;
    
    errfeast(it) = norm(sort(evslic) - ev_ref)/norm(ev_ref);
    
    fprintf('FEAST eig finishes in %.2e sec.\n', timfeast(it) );
    fprintf('The rel err is %.2e.\n', errfeast(it) );
    
    fprintf('---------------------------------------------------------\n');
    
    %--------------------------------------------------------------------
    % Zolo eig
    tic;
    opt = [];
    opt.r = [4 4];
    opt.nc = nc;
    opt.verbose = true;
    opt.itsolacc = 1e-6;
    opt.itmax = 1;
    aint = [ev(mid) ev(mid+1)];
    bint = [ev(mid+nc) ev(mid+nc+1)];
    evslic = zoloeigs(A,aint,bint,opt);
    timzolo(it) = toc;
    
    errzolo(it) = norm(sort(evslic) - ev_ref)/norm(ev_ref);
    
    fprintf('Zolo eig finishes in %.2e sec.\n', timzolo(it) );
    fprintf('The rel err is %.2e.\n', errzolo(it) );
    
    fprintf('---------------------------------------------------------\n');
end

%%

figure(1)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
loglog(siz(idxeig),timeig(idxeig),'*');
loglog(siz(idxeig),timfeast(idxeig),'*');
loglog(siz(idxeig),timzolo(idxeig),'*');
sizline = min(siz(idxeig)):.1:max(siz(idxeig));
plot(sizline,min(timzolo(idxeig))/min(siz(idxeig))*sizline,':','Linewidth',2);
plot(sizline,min(timzolo(idxeig))/min(siz(idxeig))^2*sizline.^2,':','Linewidth',2);
plot(sizline,min(timzolo(idxeig))/min(siz(idxeig))^3*sizline.^3,':','Linewidth',2);
xlabel('Size');
ylabel('Time');
legend('Eig time','FEAST time','ZoloEig time','linear ref','quadratic ref','cubic ref');

figure(2)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
loglog(nz(idxeig),timeig(idxeig),'*');
loglog(nz(idxeig),timfeast(idxeig),'*');
loglog(nz(idxeig),timzolo(idxeig),'*');
loglog(nz(idxeig),min(timzolo(idxeig))/min(nz(idxeig))*nz(idxeig),':','Linewidth',2);
loglog(nz(idxeig),min(timzolo(idxeig))/min(nz(idxeig))^2*nz(idxeig).^2,':','Linewidth',2);
loglog(nz(idxeig),min(timzolo(idxeig))/min(nz(idxeig))^3*nz(idxeig).^3,':','Linewidth',2);
xlabel('Non-Zeros');
ylabel('Time');
legend('Eig time','FEAST time','ZoloEig time','linear ref','quadratic ref','cubic ref');

figure(3)
set(gca,'YScale','log','FontSize',16);
 hold all;
suberrzolo = errzolo(idxeig);
suberrfeast = errfeast(idxeig);
[~,idx] = sort(suberrzolo);
semilogy(suberrfeast(idx),'*');
semilogy(suberrzolo(idx),'*');
ylabel('Relative Error');
legend('FEAST RelErr','ZoloEig RelErr');

figure(4)
set(gca,'YScale','log','FontSize',16);
hold all;
suberrzolo = errzolo(idxeig);
suberrfeast = errfeast(idxeig);
[~,idx] = sort(suberrfeast);
semilogy(suberrfeast(idx),'*');
semilogy(suberrzolo(idx),'*');
ylabel('Relative Error');
legend('FEAST RelErr','ZoloEig RelErr');
