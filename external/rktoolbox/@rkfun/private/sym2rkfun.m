function r = sym2rkfun(str)
%SYM2RKFUN    Construct an RKFUN from rational symbolic expression.

  % The conversion is performed by first finding roots and poles.  
  variablename = symvar(str);
  
  if length(variablename) == 1
    syms(char(variablename{1}));
  else
    error(['SYM2RKFUN: Failed to evaluate input string. Only one' ...
           ' independent variable allowed.']);    
  end

  try   ratf = eval(str);
  catch error(['SYM2RKFUN: Failed to evaluate input string. Must' ...
               ' be a rational function.']); 
  end 
    
  orig_state = warning;
  warning('off','all');
  
  try   rts = double(solve(simplify(ratf)));
  catch rts = [];  end
  
  try   pls = double(solve(simplify(1/ratf)));
  catch pls = [];  end  
  
  warning(orig_state);

  p = 0.601487492996570;
  
  try   valp = double(subs(ratf, sym(p)));
  catch error('SYM2RKFUN: Conversion failed.'); end
  
  r = nodes2rkfun(rts, pls);
  r = (valp/feval(r, p))*r;

  % Point-wise check if conversion was successful.
  p = -1.106081920075201;
  valp = double(subs(ratf, sym(p)));
  if abs(feval(r,p) - valp) > 1e-14*abs(valp)
    error('SYM2RKFUN: Conversion failed.');
  end
  
end