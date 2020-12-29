function f = util_dtnfun(h, c, x)
% UTIL_DTNFUN    Discrete Dirichlet-to-Neumann function.
%
% Returns a function handle to the discrete Dirichlet-to-Neumann function 
% for step size h and with variable-coefficient c. The first entry in c 
% is the coefficient for the infinite extension, and then c(2),c(3),... 
% are variable coefficients towards the x = 0 interface.
% If c = 0, then f(x) = sqrt(x+(hx/2)^2).
% The handle f corresponds to the function defined by the continued fraction
%
% h(z+c(1))                          1
% --------- + -----------------------------------------------------------------
%     2                                        1
%             h + -------------------------------------------------------------
%                                                    1
%                 h(z+c(2)) + -------------------------------------------------
%                                                           1  
%                             h + ... + ---------------------------------------
%                                                                1
%                                       h(z+c(end)) + -------------------------
%                                                                  1
%                                                     h + ---------------------
%                                                         hz/2+sqrt(z+(hz/2)^2)
%                                                      
% This corresponds to formula (2.3) in 
%
%   V. Druskin, S. Guettel, L. Knizhnerman. Compressing variable-coefficient 
%   exterior Helmholtz problems via RKFIT, MIMS EPrint 2016.53 
%   (<http://eprints.ma.man.ac.uk/2511/>), Manchester Institute for 
%   Mathematical Sciences, The University of Manchester, UK, 2016.

if nargin == 2,
    f = @(x) sqrt(((x+c(1)) + (h*(x+c(1))/2).^2));
    for j = 2:length(c),
        %f = @(x) f(x) + h*(x+c(j-1))/2;
        %f = @(x) h + 1./f(x);
        %f = @(x) h*(x+c(j))/2 + 1./f(x);
        f = @(x) h*(x+c(j))/2 + 1./(h + 1./(f(x) + h*(x+c(j-1))/2));
    end
end
if nargin == 3,
    f = sqrt(((x+c(1)) + (h*(x+c(1))/2).^2));
    for j = 2:length(c),
        f = h*(x+c(j))/2 + 1./(h + 1./(f + h*(x+c(j-1))/2));
    end    
end