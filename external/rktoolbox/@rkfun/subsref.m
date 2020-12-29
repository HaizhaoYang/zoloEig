function w = subsref(obj, varargin)
%SUBSREF   Evaluate an RKFUN (calls feval).
  
if nargin == 2 && length(varargin{1}) == 1 && strcmp(varargin{1}.type, '.')
  if strcmp(varargin{1}.subs, 'K'), w = obj.K; end
  if strcmp(varargin{1}.subs, 'H'), w = obj.H; end
  if strcmp(varargin{1}.subs, 'coeffs'), w = obj.coeffs; end
  if strcmp(varargin{1}.subs, 'k'), w = obj.k; end
  return
end

if nargin == 2 && length(varargin{1}) == 2   
  tmp = varargin{1};
  if and(strcmp(tmp(1).type, '.'), ...
         and(strcmp(tmp(1).subs, 'coeffs'), ...
             strcmp(tmp(2).type, '()')))
    w = cell2mat(tmp(2).subs);
    tmp = obj.coeffs;
    w = tmp(w);
    return
  end        
end


if length(varargin{1}.subs) == 1
  w = feval(obj, varargin{1}.subs{1});
  return
end

if length(varargin{1}.subs) == 2
  w = feval(obj, varargin{1}.subs{1}, varargin{1}.subs{2});
  return
end

end
