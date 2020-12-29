function AB = util_linearize_nlep(A, Sigma, Xi, tol, Nmax, cyclic)
% UTIL_LINEARIZE_NLEP    Constructs a pencil for the linearization
%                        of a nonlinear eigenvalue problem.
% This is the same as UTIL_LINEARISE_NLEP (with an S)!

% Call the English function (has better manners):
AB = util_linearise_nlep(A, Sigma, Xi, tol, Nmax, cyclic);

