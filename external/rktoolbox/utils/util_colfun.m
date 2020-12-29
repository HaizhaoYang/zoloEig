function W = util_colfun( fun, V )
% UTIL_COLFUN   Apply a function column-wise.
%
%    W = UTIL_COLFUN(FUN, V) applies the function FUN to every column of V.

W = 0*V;
for j = 1:size(V,2),
   W(:,j) = fun(V(:,j)); 
end

end

