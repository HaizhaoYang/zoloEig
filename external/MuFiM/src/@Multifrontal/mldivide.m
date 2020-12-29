function C = mldivide(A,B)
% MRDIVIDE C = A\B
%   C = A\B matrix B is left divided by matrix A which is factorized.
%
%   See also SYMBOLMF, MULTIFRONTAL.

%   Copyright 2016 Yingzhou Li, Stanford University

if isa(A,'Multifrontal')
    
    MF = A;
    C = B;
    
    if MF.symm == 1
        LeftDivSymmUp  (MF.Nnode,MF.idxtree,MF.Ltree);
        LeftDivSymmDiag(MF.Nnode,MF.idxtree,MF.Dtree);
        LeftDivSymmDown(MF.Nnode,MF.idxtree,MF.Ltree);
    elseif MF.symm == 2
        LeftDivSymmUp  (MF.Nnode,MF.idxtree,MF.Ltree);
        LeftDivSymmDown(MF.Nnode,MF.idxtree,MF.Utree);
    end
    
else
    error('Multifrontal as a numeriter has not been implemented yet');
end

%=====================================================================
    function LeftDivSymmUp(Nnode,idxtree,Ltree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(idx,:);
            Cidx = Ltree(it).Matinv*Cidx;
            C(idx,:) = Cidx;
            C(actidx,:) = C(actidx,:) - Ltree(it).AMatinv*Cidx;
        end
        
    end

    function LeftDivSymmDiag(Nnode,idxtree,Dtree)
        
        for it = 1:Nnode
            idx = idxtree(it).idx;
            C(idx,:) = Dtree(it).Matinv*C(idx,:);
        end
        
    end

    function LeftDivSymmDown(Nnode,idxtree,Utree)
        
        for it = Nnode:-1:1
            idx = idxtree(it).idx;
            actidx = idxtree(it).actidx;
            Cidx = C(idx,:);
            Cidx = Cidx - Utree(it).AMatinv'*C(actidx,:);
            C(idx,:) = Utree(it).Matinv'*Cidx;
        end
        
    end

end