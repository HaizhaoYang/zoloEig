function s = size(obj, dim)
%SIZE    Return the size of an RKFUN.
  
  if obj.iscolumn, s = [inf, numel(obj)];
  else             s = [numel(obj), inf]; end

  if nargin == 2,  s = s(dim); end

end