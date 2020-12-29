function MF = Factorization(MF,A)
% FACTORIZATION Multifrontal factorization
%   MF = FACTORIZATION(MF,A) factorizes the matrix A as a data-sparse
%   multiplication of lower or upper trangular matrices. If A is a
%   numerically symmetric matrix, A is factorized as
%       L_1 L_2 ... L_k D L_k^T ... L_2^T L_1^T,
%   where L_i is lower trangular matrix. If A is a pattern symmetric
%   matrix, A is factorized as
%       A = L_1 L_2 ... L_k U_k ... U_2 U_1,
%   where L_i is lower trangular matrix and U_j is upper trangular matrix.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

P = zeros(size(A,1),1);
it = 0;

if MF.symm == 1
    MF.Ltree = repmat(struct( 'Mat',[], 'Matinv',[], 'AMatinv',[] ), ...
        MF.Nnode, 1);
    MF.Dtree = repmat(struct( 'Mat',[], 'Matinv',[] ), MF.Nnode, 1);
    MF.idxtree = repmat(struct( 'idx',[], 'actidx',[] ), MF.Nnode, 1);
    SymmFactorizationRecursion(MF.symboltree);
elseif MF.symm == 2
    MF.Ltree = repmat(struct( 'Mat',[], 'Matinv',[], 'AMatinv',[] ), ...
        MF.Nnode, 1);
    MF.Utree = repmat(struct( 'Mat',[], 'Matinv',[], 'AMatinv',[] ), ...
        MF.Nnode, 1);
    MF.idxtree = repmat(struct( 'idx',[], 'actidx',[] ), MF.Nnode, 1);
    PatSymmFactorizationRecursion(MF.symboltree);
end

%====================================================================
    function [extidx,Aupdate] = ...
            SymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [lidx,lA] = SymmFactorizationRecursion(Stree.ltree);
            
            [ridx,rA] = SymmFactorizationRecursion(Stree.rtree);
            
            [cidx,cA] = MergeUpdate(lidx,lA,ridx,rA);
        else
            cidx = [];
            cA = [];
        end
        
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        
        [Aidx,Aact,Aupdate] = ExtUpdate(idx,actidx,cidx,cA);
        
        [L,D] = ldl(full(Aidx));
        Linv = inv(L);
        Dinv = inv(D);
        
        ALDinv = Aact*Linv'*Dinv;
        
        extidx = actidx;
        Aupdate = Aupdate - ALDinv*Linv*Aact';
        
        it = it + 1;
        MF.idxtree(it).idx    = idx;
        MF.idxtree(it).actidx = actidx;
        MF.Ltree(it).Mat      = L;
        MF.Ltree(it).Matinv   = Linv;
        MF.Ltree(it).AMatinv  = ALDinv;
        MF.Dtree(it).Mat      = D;
        MF.Dtree(it).Matinv   = Dinv;
        
    end

%====================================================================
    function [extidx,Aupdate] = ...
            PatSymmFactorizationRecursion(Stree)
        
        if strcmpi(Stree.type,'node')
            [lidx,lA] = PatSymmFactorizationRecursion(Stree.ltree);
            
            [ridx,rA] = PatSymmFactorizationRecursion(Stree.rtree);
            
            [cidx,cA] = MergeUpdate(lidx,lA,ridx,rA);
        else
            cidx = [];
            cA = [];
        end
        
        idx = Stree.idx;
        actidx = Stree.actidx;
        
        [Aidx,Aactidx,Aidxact,Aupdate] = PatExtUpdate(idx,actidx,cidx,cA);
        
        [L,U] = lu(Aidx);
        Linv = inv(L);
        Uinv = inv(U);
        
        AUinv = Aactidx*Uinv;
        ALinv = (Linv*Aidxact)';
        
        extidx = actidx;
        Aupdate = Aupdate - AUinv*ALinv';
        
        it = it + 1;
        MF.idxtree(it).idx    = idx;
        MF.idxtree(it).actidx = actidx;
        MF.Ltree(it).Mat      = L;
        MF.Ltree(it).Matinv   = Linv;
        MF.Ltree(it).AMatinv  = AUinv;
        MF.Utree(it).Mat      = U';
        MF.Utree(it).Matinv   = Uinv';
        MF.Utree(it).AMatinv  = ALinv;
        
    end

%====================================================================
    function [idx,A] = MergeUpdate(lidx,lA,ridx,rA)
        
        llen = length(lidx);
        rlen = length(ridx);
        [idx,~,IC] = unique([lidx,ridx],'sorted');
        lI = IC(1:llen);
        rI = IC(llen+1:llen+rlen);
        A = zeros(length(idx));
        A(lI,lI) = A(lI,lI) + lA;
        A(rI,rI) = A(rI,rI) + rA;
        
    end

    function [Aidx,Aact,Aupdate] = ExtUpdate(idx,actidx,cidx,cA)
        
        ilen = length(idx);
        alen = length(actidx);
        [I,J,V] = find(A(:,idx));
        
        P(idx) = 1:ilen;
        Aidx = zeros(ilen,ilen);
        iI = ismembc(I,idx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aidx(P(Isub) + (Jsub - 1)*ilen) = Vsub;
        
        P(actidx) = 1:alen;
        Aact = zeros(alen,ilen);
        iI = ismembc(I,actidx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aact(P(Isub) + (Jsub - 1)*alen) = Vsub;
        
        Aupdate = zeros(length(actidx));
        
        iicI = ismembc(idx,cidx);
        cicI = ismembc(cidx,idx);
        aacI = ismembc(actidx,cidx);
        cacI = ismembc(cidx,actidx);
        Aidx(iicI,iicI) = Aidx(iicI,iicI) + cA(cicI,cicI);
        Aact(aacI,iicI) = Aact(aacI,iicI) + cA(cacI,cicI);
        Aupdate(aacI,aacI) = cA(cacI,cacI);
        
    end

    function [Aidx,Aactidx,Aidxact,Aupdate] = ...
            PatExtUpdate(idx,actidx,cidx,cA)
        
        ilen = length(idx);
        alen = length(actidx);
        [I,J,V] = find(A(:,idx));
        
        P(idx) = 1:ilen;
        Aidx = zeros(ilen,ilen);
        iI = ismembc(I,idx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aidx(P(Isub) + (Jsub - 1)*ilen) = Vsub;
        
        P(actidx) = 1:alen;
        Aactidx = zeros(alen,ilen);
        iI = ismembc(I,actidx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aactidx(P(Isub) + (Jsub - 1)*alen) = Vsub;
        
        P(idx) = 1:ilen;
        [I,J,V] = find(A(:,actidx));
        Aidxact = zeros(ilen,alen);
        iI = ismembc(I,idx);
        Isub = I(iI);
        Jsub = J(iI);
        Vsub = V(iI);
        Aidxact(P(Isub) + (Jsub - 1)*ilen) = Vsub;
        
        Aupdate = zeros(length(actidx));
        
        iicI = ismembc(idx,cidx);
        cicI = ismembc(cidx,idx);
        aacI = ismembc(actidx,cidx);
        cacI = ismembc(cidx,actidx);
        Aidx(iicI,iicI) = Aidx(iicI,iicI) + cA(cicI,cicI);
        Aactidx(aacI,iicI) = Aactidx(aacI,iicI) + cA(cacI,cicI);
        Aidxact(iicI,aacI) = Aidxact(iicI,aacI) + cA(cicI,cacI);
        Aupdate(aacI,aacI) = cA(cacI,cacI);
        
    end

end