maxNumCompThreads(1);

siznpc = 2.^(3:4);
K = 4*ones(1,length(siznpc));

siz = (K.*siznpc).^2;

nmat = length(siz);

trail = 10;

nc = 96;

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
        [v,ev] = eigs(A,[],128);
        timeig(ittrail,it) = toc;
        ev = diag(ev);
        [ev,idx] = sort(ev);
        v = v(:,idx);
        
        fprintf('Matlab eigs finishes in %.2e sec.\n', timeig(ittrail,it));
        
        fprintf('-----------------------------------------------------\n');
        
        a = ev(1)-10;
        b = ev(nc)*0.4+ev(nc+1)*0.6;
        
        idx = ev>a & ev<b;
        ev_ref = ev(idx);
        v_ref = v(:,idx);
        
        
        %         %-----------------------------------------------------------------
        %         % FEAST eig
        %         opt.quarType = 'GL';
        %         tic;
        %         evslic = rkfeast(A,a,b,nc,16);
        %         timfeast(ittrail,it) = toc;
        %
        %         errfeast(ittrail,it) = norm(sort(evslic) - ev_ref)/norm(ev_ref);
        %
        %         fprintf('FEAST eig finishes in %.2e sec.\n', timfeast(ittrail,it) );
        %         fprintf('The rel err is %.2e.\n', errfeast(ittrail,it) );
        %
        %         fprintf('-----------------------------------------------------\n');
        
        %-----------------------------------------------------------------
        % Zolo eig
        tic;
        opt = [];
        opt.verbose = 1;
        opt.r = [4 4];
        opt.nc = nc;
        aint = [a ev(1)];
        bint = [ev(nc) ev(nc+1)];
        evslic = zoloeigs(A,aint,bint,opt);
        timzolo(ittrail,it) = toc;
        
        errzolo(ittrail,it) = norm(sort(real(evslic))-ev_ref)/norm(ev_ref);
        
        fprintf('Zolo eig finishes in %.2e sec.\n', timzolo(ittrail,it) );
        fprintf('The rel err is %.2e.\n', errzolo(ittrail,it) );
        
        fprintf('-----------------------------------------------------\n');
        
    end
end

save('Eigs_FD_Data','timeig','timfeast','timzolo','siz', ...
    'errfeast','errzolo');

%%

figure(1)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
errorbar(siz,mean(timeig), ...
    mean(timeig)-min(timeig),max(timeig)-mean(timeig), ...
    'Linewidth',2);
%errorbar(siz,mean(timfeast), ...
%    mean(timfeast)-min(timfeast),max(timfeast)-mean(timfeast), ...
%    'Linewidth',3);
errorbar(siz,mean(timzolo), ...
    mean(timzolo)-min(timzolo),max(timzolo)-mean(timzolo), ...
    'Linewidth',2);
plot(siz,min(min(timzolo))/min(siz)*siz,':','Linewidth',2);
plot(siz,min(min(timzolo))/min(siz)^2*siz.^2,':','Linewidth',2);
xlabel('Size');
ylabel('Time');
legend('Eigs time','ZoloEig time','linear ref','quadratic ref');

figure(2)
set(gca,'XScale','log','YScale','log','FontSize',16);
hold all;
%errorbar(siz,mean(errfeast), ...
%    mean(errfeast)-min(errfeast),max(errfeast)-mean(errfeast), ...
%    'Linewidth',2);
errorbar(siz,mean(errzolo), ...
    mean(errzolo)-min(errzolo),max(errzolo)-mean(errzolo), ...
    'Linewidth',2);
xlabel('Size');
ylabel('Relative Error');
legend('ZoloEig RelErr');
