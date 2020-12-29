function [x, success] = util_cplxpair(x, tol)
% UTIL_CPLXPAIR    Sorts a vector keeping together
%                  exact complex-conjugate pairs.  
%
% [x, success] = util_cplxpair(x, tol) is used to sort the array x 
% in complex conjugate pairs if possible. Matlab's function cplxpair
% does the same but it will break with an error if no pairing
% exists. If the pairing is not possible the array x remains
% unaltered, and the flag success is set to false, otherwise it set
% to true. If the pairing exists, x will be ordered by cplxpair.  
%
% See also cplxpair.
  
  if nargin == 1, tol = 100*eps; end

  try
    x = cplxpair(x, tol);   
    success = true;
  catch
    success = false;
  end
  
end
