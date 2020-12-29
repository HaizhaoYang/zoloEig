classdef Multifrontal < SymbolMF
% MULTIFRONTAL Class for Multifrontal factorization
%   MF = MULTIFRONTAL(A) first generates the symbolic factorization for
%   general sparse matrices A and then factorizes the matrix A as a
%   data-sparse multiplication of lower or upper trangular matrices. If A
%   is a numerically symmetric matrix, A is factorized as
%       L_1 L_2 ... L_k D L_k^T ... L_2^T L_1^T,
%   where L_i is lower trangular matrix. If A is a pattern symmetric
%   matrix, A is factorized as
%       A = L_1 L_2 ... L_k U_k ... U_2 U_1,
%   where L_i is lower trangular matrix and U_j is upper trangular matrix.
%
%   MF = MULTIFRONTAL(A,SMF) factorizes the matrix A based on the symbolic
%   factorization given in SMF.
%
%   MF = MULTIFRONTAL(A,cutoff) generates the symbolic factorization with
%   the cutoff and then factorizes the matrix A based on the symbolic
%   factorization given in SMF.
%
%   The class is a derived class of SymbolMF and also includes the
%   following extra properties.
%       Name         Explaination  
%       Ltree        Tree stores the lower trangular matrices
%       Utree        Tree stores the upper trangular matrices
%       Dtree        Tree stores the diagonal matrices
%
%   See also SYMBOLMF, BITREEPARTITION.

%   Copyright 2016 Yingzhou Li, Stanford University

    properties
        Ltree
        Utree
        Dtree
        idxtree
    end
    methods
        function MF = Multifrontal(A,SMF)
            
            if nargin < 2
                SMF = SymbolMF(A);
            elseif nargin == 2
                if ~isa(SMF,'SymbolMF')
                    SMF = SymbolMF(A,SMF);
                end
            end
            
            MF.M          = SMF.M;
            MF.N          = SMF.N;
            MF.symm       = SMF.symm;
            MF.symboltree = SMF.symboltree;
            MF.Nnode      = SMF.Nnode;
            
            MF = Factorization(MF,A);
            
        end
    end
end