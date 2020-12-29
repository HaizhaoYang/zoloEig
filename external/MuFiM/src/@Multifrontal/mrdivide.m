function C = mrdivide(A,B)
% MRDIVIDE C = A/B
%   C = A/B matrix A divided by matrix B which is factorized.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

if isa(B,'Multifrontal')
    
    MF = B;
    C = A;
    
    if MF.symm == 1
        RightDivSymmUp  (MF.Nnode,MF.idxtree,MF.Ltree);
        RightDivSymmDiag(MF.Nnode,MF.idxtree,MF.Dtree);
        RightDivSymmDown(MF.Nnode,MF.idxtree,MF.Ltree);
    elseif MF.symm == 2
        RightDivSymmUp  (MF.Nnode,MF.idxtree,MF.Utree);
        RightDivSymmDown(MF.Nnode,MF.idxtree,MF.Ltree);
    end
    
else
    error('Multifrontal as a numeriter has not been implemented yet');
end

%=====================================================================
    function RightDivSymmUp(Nnode,idxtree,Utree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(:,idx);
            Cidx = Cidx*Utree(it).Matinv';
            C(:,idx) = Cidx;
            C(:,actidx) = C(:,actidx) - Cidx*Utree(it).AMatinv';
        end
        
    end

    function RightDivSymmDiag(Nnode,idxtree,Dtree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            C(:,idx) = C(:,idx)*Dtree(it).Matinv;
        end
        
    end

    function RightDivSymmDown(Nnode,idxtree,Ltree)
        
        for it = Nnode:-1:1
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(:,idx);
            Cidx = Cidx - C(:,actidx)*Ltree(it).AMatinv;
            Cidx = Cidx*Ltree(it).Matinv;
            C(:,idx) = Cidx;
        end
        
    end

end