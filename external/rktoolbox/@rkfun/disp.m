function disp(obj)
%DISP   Display information about an RKFUN.

  fprintf('\tRKFUN object of type (%d, %d).\n', type(obj));

  if isreal(obj.K) && isreal(obj.H), f = 'Real';
  else                               f = 'Complex'; end
  s = [f '-valued Hessenberg pencil (H, K) of size '];
  fprintf('\t%s%d-by-%d.\n', s, size(obj.H, 1), size(obj.H, 2));

  c = obj.coeffs;
  if isa(c, 'sym')
    fprintf('\tVariable precision arithmetic (VPA) activated.\n');
  elseif isa(c, 'mp')
    fprintf(['\tMultiple precision arithmetic (ADVANPIX)' ...
             ' activated.\n']);
  end

  c = double(c);
  L = min(length(c), 5);
  if isreal(obj.coeffs)
    fprintf('\tcoeffs = [');
    if L > 1, fprintf('%1.3f, ', c(1:L-1)); end
    fprintf('%1.3f',   c(L));
  else
    fprintf('\t|coeffs| = [');
    if L > 1, fprintf('%1.3f, ', abs(c(1:L-1))); end
    fprintf('%1.3f',   abs(c(L)));
  end
  L = length(c) - L;
  if L <= 0, fprintf(']\n');
  else       fprintf(', ...]\n', L); end

end % function
