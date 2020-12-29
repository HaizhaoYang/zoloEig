function C = mtimes(A,B)
% MTIMES C = A B
%   C = A*B multiplies matrix A and B together when either of them is
%   factorized.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

if isa(A,'Multifrontal')
    
    MF = A;
    C = B;
    
    if MF.symm == 1
        LeftMulSymmUp  (MF.Nnode,MF.idxtree,MF.Ltree);
        LeftMulSymmDiag(MF.Nnode,MF.idxtree,MF.Dtree);
        LeftMulSymmDown(MF.Nnode,MF.idxtree,MF.Ltree);
    elseif MF.symm == 2
        LeftMulSymmUp  (MF.Nnode,MF.idxtree,MF.Utree);
        LeftMulSymmDown(MF.Nnode,MF.idxtree,MF.Ltree);
    end
    
elseif isa(B,'Multifrontal')
    
    MF = B;
    C = A;
    
    if MF.symm == 1
        RightMulSymmUp  (MF.Nnode,MF.idxtree,MF.Ltree);
        RightMulSymmDiag(MF.Nnode,MF.idxtree,MF.Dtree);
        RightMulSymmDown(MF.Nnode,MF.idxtree,MF.Ltree);
    elseif MF.symm == 2
        RightMulSymmUp  (MF.Nnode,MF.idxtree,MF.Ltree);
        RightMulSymmDown(MF.Nnode,MF.idxtree,MF.Utree);
    end
    
end

%=====================================================================
    function LeftMulSymmUp(Nnode,idxtree,Utree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(idx,:);
            Cidx = Utree(it).Mat'*Cidx;
            Cidx = Cidx + Utree(it).AMatinv'*C(actidx,:);
            C(idx,:) = Cidx;
        end
        
    end

    function LeftMulSymmDiag(Nnode,idxtree,Dtree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            C(idx,:) = Dtree(it).Mat*C(idx,:);
        end
        
    end

    function LeftMulSymmDown(Nnode,idxtree,Ltree)
        
        for it = Nnode:-1:1
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(idx,:);
            C(actidx,:) = C(actidx,:) + Ltree(it).AMatinv*Cidx;
            Cidx = Ltree(it).Mat*Cidx;
            C(idx,:) = Cidx;
        end
        
    end

%=====================================================================
    function RightMulSymmUp(Nnode,idxtree,Ltree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(:,idx);
            Cidx = Cidx*Ltree(it).Mat;
            Cidx = Cidx + C(:,actidx)*Ltree(it).AMatinv;
            C(:,idx) = Cidx;
        end
        
    end

    function RightMulSymmDiag(Nnode,idxtree,Dtree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            C(:,idx) = C(:,idx)*Dtree(it).Mat;
        end
        
    end

    function RightMulSymmDown(Nnode,idxtree,Utree)
        
        for it = Nnode:-1:1
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(:,idx);
            C(:,actidx) = C(:,actidx) + Cidx*Utree(it).AMatinv';
            Cidx = Cidx*Utree(it).Mat';
            C(:,idx) = Cidx;
        end
        
    end

end