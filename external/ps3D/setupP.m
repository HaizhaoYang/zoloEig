function C = setupP(N,npw,P)

    GRD = reshape(1:N^3,[N,N,N]);
    REM = ones(N,N,N);
    Lcur = P;
    
    maxlen = N^3/4;
    clear C;
    C(1:maxlen) = struct('rdidx', 0, 'skidx', 0,'Arrinv', 0, 'X', 0, 'Y', 0, 'S', 0);
    cntc = 0;
    
    NL = round(log2(N/npw))+1;
    for ell=1:NL
        W = 2^(ell-1)*npw;
        nck = N/W;
        RDS = zeros(N,N,N);
        
        cnta = 0;
        cntc_start = cntc;
        tic;
        for K=1:nck
            for J=1:nck
                for I=1:nck
                    %---------
                    if(    ell<NL)
                        ia = (I-1)*W+1;        ib = I*W+1;                is = ia:ib;
                        if(ib==N+1)
                            ib=1; is(end)=1;
                        end
                        ja = (J-1)*W+1;        jb = J*W+1;                js = ja:jb;
                        if(jb==N+1)
                            jb=1; js(end)=1;
                        end
                        ka = (K-1)*W+1;        kb = K*W+1;                ks = ka:kb;
                        if(kb==N+1)
                            kb=1; ks(end)=1;
                        end
                        
                        [ii,jj,kk] = ndgrid(is,js,ks);
                        curGRD = GRD(is,js,ks);
                        curREM = REM(is,js,ks);
                        curLCL = zeros(size(curGRD));
                        curLCL(find(curREM==1)) = 1:numel(find(curREM==1)); %local ordering
                        GALL = curGRD(find(curREM==1));
                        inn = (ii~=ia & ii~=ib & jj~=ja & jj~=jb & kk~=ka & kk~=kb);
                        rdidx = curGRD(find(curREM==1 & inn==1));
                        skidx = curGRD(find(curREM==1 & inn==0));
                    elseif(ell==NL)
                        rdidx = find(REM==1);
                        skidx = zeros(0,1);
                    end
                    Arr = full(Lcur(rdidx,rdidx));                %if(ell>1) Arr=full(Arr); end;
                    Asr = full(Lcur(skidx,rdidx));
                    Ars = full(Lcur(rdidx,skidx));
                    Arrinv = inv((Arr));
                    
                    X = Asr*Arrinv;
                    Y = Arrinv*Ars;
                    S = -X*Ars;
                    C(cntc+1) = struct('rdidx', rdidx, 'skidx', skidx, ...
                                       'Arrinv', Arrinv, 'X', X, 'Y', Y, 'S', S);
                    cntc = cntc+1;
                    
                    cnta = cnta + numel(S);
                    RDS(rdidx)=1;
                end
            end
        end
        toc;
        
        Ia = zeros(1,cnta);
        Ja = zeros(1,cnta);
        Sa = zeros(1,cnta);
        cnta = 0;
        for g=(cntc_start+1):cntc
            [It,Jt] = ndgrid(C(g).skidx);
            S = C(g).S; C(g).S=[];
            It = It(:);
            Jt = Jt(:);
            St = S(:);
            num = numel(S);
            Ia(cnta+[1:num]) = It;
            Ja(cnta+[1:num]) = Jt;
            Sa(cnta+[1:num]) = St;
            cnta = cnta + num;            clear It Jt St S;
        end
        %update the matrix
        tmp = find(RDS==1);        %OLDREM = REM;
        REM(tmp) = 0;
        Sch = sparse(Ia,Ja,Sa,N^3,N^3); %block diagonal
        clear Ia Ja Sa;
        
        [Ib,Jb,Sb] = find(Lcur); clear Lcur;
        act1 = find(REM(Ib)==1&REM(Jb)==1); %remains
                                            %act2 = find(REM(Ib)==0&REM(Jb)==0); %deleted
        act = [act1];
        Ib = Ib(act);
        Jb = Jb(act);
        Sb = Sb(act);
        clear act;
        tmp = sparse(Ib,Jb,Sb,N^3,N^3); %the ones restricted to the REM ones
        clear Ib Jb Sb;
        
        Lcur = tmp + Sch; clear tmp Pat Sch;
    end
    
    C = C(1:cntc);
    
end

