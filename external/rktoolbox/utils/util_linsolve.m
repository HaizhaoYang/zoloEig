function y = util_linsolve(M,v)
% UTIL_LINSOLVE    Solves linear system of equations.
%
% y = util_linsolve(M,v) solves the system M*v = y, where the LU
% factors of the coefficient matrix M are reused when the same
% input matrix is identified by its (1,1) element.
%
% To delete the LU factors call UTIL_LINSOLVE with empty argument
% list.

    persistent sol indm

    verbose = 1;

    if nargin < 2,
        sol = []; indm = [];
        if verbose, disp('util_linsolve reset'); end;
        return
    end

    indm(1) = NaN; % dummy element
    m = M(1,1);
    ind = find(indm == m,1,'first');

    if issparse(M),
       if isempty(ind),
           ind = length(indm) + 1;
           if verbose, disp(['util_linsolve: sparse factorization nr ' num2str(ind-1)]); end;
           [sol(ind).L,sol(ind).U,sol(ind).P,sol(ind).Q,sol(ind).R] = lu(M);
           indm(ind) = m;
       end
       y = sol(ind).Q * (sol(ind).U \ (sol(ind).L \ (sol(ind).P * (sol(ind).R \ v))));
    else
       if isempty(ind),
           ind = length(indm) + 1;
           if verbose, disp(['util_linsolve: full factorization nr ' num2str(ind-1)]); end;
           [sol(ind).L,sol(ind).U,sol(ind).p] = lu(M,'vector');
           indm(ind) = m;
       end
       y = sol(ind).U \ (sol(ind).L \ v(sol(ind).p,:));
    end

end
