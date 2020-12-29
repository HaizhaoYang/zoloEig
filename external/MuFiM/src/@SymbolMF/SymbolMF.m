classdef SymbolMF
% SYMBOLMF Class for Multifrontal symbolic factorization
%   SMF = SYMBOLMF(A) generates the symbolic factorization for general
%   sparse matrices A. The graph of A is recursively bipartitioned via
%   Metis.
%
%   The class includes the following properties.
%       Name         Explaination  
%       M,N          Size of the matrix
%       symm         Symmetric flag. symm = 0 means non-symmetric,
%                                    symm = 1 means numerically symmetric
%                                    symm = 2 means pattern symmetric
%       symboltree   Symbol tree for the graph.
%
%   See also MULTIFRONTAL, BITREEPARTITION.

%   Copyright 2016 Yingzhou Li, Stanford University

    properties
        M
        N
        symm
        symboltree
        Nnode
    end
    methods
        function SMF = SymbolMF(A)
            
            if nargin < 1
                return;
            end
            
            [SMF.M,SMF.N] = size(A);
            
            if ishermitian(A)
                % numerically symmetric, factorized as LDLT
                SMF.symm = 1;
            elseif issymmetric(spones(A))
                % symbolic symmetric, factorized as LU
                SMF.symm = 2;
            else
                SMF.symm = 2;
                A = (A+A')/2;
            end
            
            if SMF.symm > 0 && issparse(A)
                [SMF.symboltree, SMF.Nnode] = BiTreePartition(A);
            else
                SMF.symboltree.type = 'leaf';
                SMF.symboltree.idx = 1:SMF.M;
                SMF.symboltree.actidx = [];
                SMF.Nnode = 1;
            end
        end
    end
end