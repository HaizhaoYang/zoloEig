function check = test_givens()

  tol  = 1e-13;
  
  x{1} = [1; 2];
  x{2} = [0; 2];
  x{3} = [1; 0];
  x{4} = [1i; 0];
  x{5} = [1e-9; 2i+3];
  x{6} = [3+1i; 2i-3];

  err1 = 0;
  for j = 1:length(x)
    xx = x{j};
    [s, c] = util_givens(xx(1), xx(2));
    G = [c -s; conj(s) c];
    
    err1 = err1 + abs((c.^2 + abs(s).^2)-1);
    err1 = err1 + abs([0 1]*(G*xx));   
  end
  
  y{1} = [-1; 2]; 
  y{2} = [3; pi];
  y{3} = [1i; 2i+3];
  y{4} = [3+1i; 2i-3];
  y{5} = [1; 0];
  y{6} = [1i; 0];
  z = [1e-2, 1i, 5.3, -1.2, 1+2i, 5]; 
  
  err2 = 0;
  for j = 1:length(x)
    xx = x{j};
    yy = y{j};
    [s, c] = util_rat_givens(xx, yy, z(j));
    G = [c -s; conj(s) c];
    
    err2 = err2 + abs((c.^2 + abs(s).^2)-1);
    err2 = err2 + abs(([0 1]*(G*xx))/([0 1]*(G*yy))-z(j))/abs(z(j));
    
    [s, c] = util_rat_givens(yy, xx, z(j));
    G = [c -s; conj(s) c];
    
    err2 = err2 + abs((c.^2 + abs(s).^2)-1);
    err2 = err2 + abs(([0 1]*(G*yy))/([0 1]*(G*xx))-z(j))/abs(z(j));
  end

  check = [err1  err2] < tol;  
  
end