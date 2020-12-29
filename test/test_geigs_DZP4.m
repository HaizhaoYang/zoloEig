if exist('GEigs_DZP4_Data.mat','file')
    load GEigs_DZP4_Data;
else
    sit = 2;
    nmat = 7;
    siz = zeros(1,nmat);
    timfeast = zeros(1,nmat);
    timzolo = zeros(1,nmat);
    errsfeast = cell(1,nmat);
    errszolo = cell(1,nmat);
    eiggap = zeros(1,nmat);
end

for natom = sit:7
    
    load(['../data/DZP4/' num2str(natom) '/H_DZP4_' num2str(natom) ...
        '.dat']);
    H_DZP4 = eval(['H_DZP4_' num2str(natom)]);
    A = sparse(H_DZP4(:,1),H_DZP4(:,2),H_DZP4(:,3));
    B = sparse(H_DZP4(:,1),H_DZP4(:,2),H_DZP4(:,4));
    
    A = (A+A')/2;
    B = (B+B')/2;
    
    n = size(A,1);
    ncols = 93;%16*natom^3;
    
    [V, D]   = eigs(A,B,ncols+2,'SA');
    [D, ind] = sort(diag(D)); V = V(:, ind);
    
    lmin = D(1)-10;
    lmax = (1.1*D(ncols)+D(ncols+1))/2.1;
    
    opt = [];
    opt.verbose = true;
    opt.nc = ncols;
    opt.r = [3 3];
    opt.itsolacc = 1e-9;
    opt.reltol = 1e-8;
    aint = [D(1)-10 D(1)];
    bint = [D(ncols) D(ncols+1)];
    tic;
    [Uout,eigvals,relerrs] = zologeigs(A,B,aint,bint,opt);
    zolotim = toc;
    
    opt = [];
    opt.verbose = true;
    opt.reltol = 1e-8;
    tic;
    [Uoutf,eigvalsf,relerrsf] = rkfeastg(A,B,lmin,lmax,ncols,16,opt);
    feasttim = toc;
    
    %%
    figure;
    hold all;
    
    plot(relerrs,'-^');
    plot(relerrsf,'-o');
    
    legend({['ZoloEig,   ' num2str(zolotim) ' sec'], ...
        ['FEAST,   ' num2str(feasttim) ' sec']});
    
    xlabel('k');
    ylabel('Relative Error');
    
    set(gca,'YScale','log');
    
    saveas(gca, ['GEigs_DZP4_' num2str(natom) '.eps'], 'epsc' );
    
    %%
    fprintf(['%d & %d & %d & %.2e & (%d,%d) & %.2e & %d ' ...
        '& %.2e & %d & %.2e & %d & %.2e\n'], natom, n, ncols, ...
        D(ncols+1)-D(ncols), 3, 3, relerrs(end), length(relerrs), ...
        zolotim, 16, relerrsf(end), length(relerrsf), feasttim);
    
    sit = natom+1;
    siz(natom) = n;
    timfeast(natom) = feasttim;
    timzolo(natom) = zolotim;
    eiggap(natom) = D(ncols+1)-D(ncols);
    errsfeast{natom} = relerrsf;
    errszolo{natom} = relerrs;
    
    save('GEigs_DZP4_Data','sit','siz','timfeast','timzolo', ...
        'errsfeast','errszolo','eiggap');
    
end

%%
ncols = 93;
for natom = 1:6

    relerrs = errszolo{natom};
    relerrsf = errsfeast{natom};
    
    figure;
    set(gca,'YScale','log','FontSize',16);
    hold all;
    
    plot(relerrs,'-^','Linewidth',2);
    plot(relerrsf,'-o','Linewidth',2);
    
    legend({['ZoloEig,   ' num2str(timzolo(natom)) ' sec'], ...
        ['FEAST,   ' num2str(timfeast(natom)) ' sec']});
    
    xlabel('k');
    ylabel('Relative Error');
    
    fprintf(['%d & %d & %d & %.2e & (%d,%d) & %.2e & %d ' ...
        '& %.2e & %d & %.2e & %d & %.2e\n'], natom, siz(natom), ncols, ...
        eiggap(natom), 3, 3, relerrs(end), length(relerrs), ...
        timzolo(natom), 16, relerrsf(end), length(relerrsf), ...
        timfeast(natom));
    
end