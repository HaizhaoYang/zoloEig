function xi = util_ieee_poles(lower, upper, scale, m)
% UTIL_IEEE_POLES    Computes potentially good poles to initialize
%                    RKFIT with.
%
% xi = util_ieee_poles(lower, upper, scale, m) returns m (assumed
% to be even) poles of the form -z*scale + 1i*z and -z*scale - 1i*z,
% where the z are log-spaced on [10^lower, 10^upper].

  ixi = logspace(lower, upper, m/2);
  xi  = [];
  for j = 1:length(ixi)
    rxi = -ixi(j)*scale;
    xi  = [xi (rxi-1i*ixi(j)) (rxi+1i*ixi(j))];
  end

end % function