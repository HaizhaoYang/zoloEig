
%niter = 5:3:20;
niter = 2.^(8:8);
siz = nan(1,length(niter));
timSym = nan(1,length(niter));
timMF  = nan(1,length(niter));
timSol = nan(1,length(niter));
timSolMatlab = nan(1,length(niter));
for it = 1:length(niter)
    n = niter(it);
    A = getHfd2D(n,4);
    siz(it) = size(A,1);
    
    tic;
    SMF = SymbolMF(A);
    timSym(it) = toc;
    
    tic;
    AMF = Multifrontal(A,SMF);
    timMF(it) = toc
    
    X = randn(size(A,1),10);
    Y = A*X;
    tic;
    XX = AMF\Y;
    timSol(it) = toc;
    
    tic;
    XXMatlab = A\Y;
    timSolMatlab(it) = toc
    
end

figure(1)
loglog(siz,timSym,'*');
hold all;
loglog(siz,min(min(timSym))/min(siz)*siz);
loglog(siz,min(min(timSym))/min(siz)^2*siz.^2);
loglog(siz,min(min(timSym))/min(siz)^3*siz.^3);
xlabel('Size');
ylabel('Time');
legend('Sym time','linear ref','quadratic ref','cubic ref');

figure(2)
loglog(siz,timMF,'*');
hold all;
loglog(siz,min(min(timMF))/min(siz)*siz);
loglog(siz,min(min(timMF))/min(siz)^2*siz.^2);
loglog(siz,min(min(timMF))/min(siz)^3*siz.^3);
xlabel('Size');
ylabel('Time');
legend('MF time','linear ref','quadratic ref','cubic ref');

figure(3)
loglog(siz,timSol,'*');
hold all;
loglog(siz,min(min(timSol))/min(siz)*siz);
loglog(siz,min(min(timSol))/min(siz)^2*siz.^2);
loglog(siz,min(min(timSol))/min(siz)^3*siz.^3);
xlabel('Size');
ylabel('Time');
legend('Sol time','linear ref','quadratic ref','cubic ref');

figure(4)
loglog(siz,timSolMatlab,'*');
hold all;
loglog(siz,min(min(timSolMatlab))/min(siz)*siz);
loglog(siz,min(min(timSolMatlab))/min(siz)^2*siz.^2);
loglog(siz,min(min(timSolMatlab))/min(siz)^3*siz.^3);
xlabel('Size');
ylabel('Time');
legend('Matlab Sol time','linear ref','quadratic ref','cubic ref');

tim = timSym + timMF + timSol;
figure(5)
loglog(siz,tim);
hold all;
loglog(siz,timSolMatlab);
xlabel('Size');
ylabel('Time');
legend('MuFiM time','Matlab Sol time');