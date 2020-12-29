maxNumCompThreads(1);

siznpc = 2.^(4);
K = 4*ones(1,length(siznpc));

siz = (K.*siznpc).^2;

nmat = length(siz);

trail = 1;

nz = zeros(1,nmat);
timeig = zeros(trail,nmat);
timfeast = zeros(trail,nmat);
timzolo = zeros(trail,nmat);
errfeast = zeros(trail,nmat);
errzolo = zeros(trail,nmat);
gcond = zeros(1,nmat);

for it = 1:nmat
    
    A = getHfd2D(siznpc(it),K(it));
    nz(it) = nnz(A);
    
    fprintf('=====================================================\n');
    fprintf('Matrix is of size %4d with %8d non-zeros.\n', ...
        siz(it), nz(it) );
    fprintf('-----------------------------------------------------\n');
    
    for ittrail = 1:trail
        
        tic;
        ev = eig(A);
        timeig(ittrail,it) = toc;
        [ev,~] = sort(ev);
        
        fprintf('Matlab eigs finishes in %.2e sec.\n', timeig(ittrail,it));
        
        fprintf('-----------------------------------------------------\n');

        %-----------------------------------------------------------------
        % Zolo eig
        tic;
        opt = [];
        opt.verbose = 1;
        evslic = zoloeig(A,opt);
        evslic = sort(evslic);
        timzolo(ittrail,it) = toc;
        
        errzolo(ittrail,it) = norm(real(evslic)-ev)/norm(ev);
        
        fprintf('Zolo eig finishes in %.2e sec.\n', timzolo(ittrail,it) );
        fprintf('The rel err is %.2e.\n', errzolo(ittrail,it) );
        
        fprintf('-----------------------------------------------------\n');
        
    end
end

save('Eigs_FD_Data','timeig','timzolo','siz', 'errzolo');

%%

figure(1)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
errorbar(siz,mean(timeig), ...
    mean(timeig)-min(timeig),max(timeig)-mean(timeig), ...
    'Linewidth',2);
errorbar(siz,mean(timzolo), ...
    mean(timzolo)-min(timzolo),max(timzolo)-mean(timzolo), ...
    'Linewidth',2);
plot(siz,min(mean(timzolo))/min(siz)*siz,':','Linewidth',2);
plot(siz,min(mean(timzolo))/min(siz)^2*siz.^2,':','Linewidth',2);
xlabel('Size');
ylabel('Time');
legend('Eigs time','ZoloEig time','linear ref','quadratic ref');

figure(2)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
errorbar(siz,mean(errzolo), ...
    mean(errzolo)-min(errzolo),max(errzolo)-mean(errzolo), ...
    'Linewidth',2);
xlabel('Size');
ylabel('Relative Error');
legend('ZoloEig RelErr');
