function xt = util_log2lin(x,r,wz) 
% UTIL_LOG2LIN    Maps points to doubly logarithmic scale.
%
% util_loglin(x,r,wz) maps the points in x so that the outer points 
% of r = [neg-left,neg-right,pos-left,pos-right] correspond to [0,1].
% The indefinite interval [neg-right,pos-left] is mapped to a linear
% region of relative width wz.

[sm,sn] = size(x);
x = x(:).';
if nargin < 3,
    wz = 0.1; % relative width of linear region about zero
end
wn = log10(r(1)/r(2));
wp = log10(r(4)/r(3));
s = (1-wz)/(wn+wp);
wn = wn*s; wp = wp*s; % wn+wz+wp=1, proportions along [0,1] axis

xn = ( log10(-r(1)) - log10(-x(x<=r(2))) ) / log10(r(1)/r(2)) * wn;
xz = (x(x>r(2) & x<r(3)) - r(2))/(r(3)-r(2))*wz + wn;
xp = ( log10(x(x>=r(3))) - log10(r(3)) ) / log10(r(4)/r(3))*wp + wn + wz;
xt = [ xn, xz, xp ];
xt = reshape(xt,sm,sn);