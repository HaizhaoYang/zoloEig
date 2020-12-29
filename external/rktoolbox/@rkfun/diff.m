function ratfun = diff(obj)
%DIFF   Differentiate an RKFUN.
  
  % Get Ritz values and poles
  ee = eig(obj.H(1:end-1,1:end),obj.K(1:end-1,1:end));
  ee = ee(isfinite(ee));
  xi = poles(obj);
  
  % Build 3th order jordan matrix A and vector b.
  if(isreal(obj))
    [A, b] = internal_build_jordan_form(ee, 4, 'real');
  else
    [A, b] = internal_build_jordan_form(ee, 4);
  end
    
  % Vector of 3rd, 2nd and 1st order derivatives.
  v  = feval(obj, A, b);
  
  real_n    = sum(b(4:4:end));
  imag_n    = length(b)/4-real_n;
  real_end  = 4*real_n;
  
  % Rescaled derivatives.
  d3 = 3*v(1:4:real_end); % The 3rd derivative is 6*v(1:4:real_end)!  
  d2 = 2*v(2:4:real_end);
  d1 =   v(3:4:real_end);
  
  if imag_n
    index = real_end+1:8:length(v);
    index = kron(index, [1 1]) + repmat([0, 1], 1, length(index));    
    d3 = [d3; 3*v(index  )];
    d2 = [d2; 2*v(index+2)];
    d1 = [d1;   v(index+4)];
  end

  % Build 2nd order jordan matrix A, F, and vector b.
  AA = sparse([]);
  FF = sparse([]);
  bb = [];
  for j = 1:real_n
    jb = (j-1)*4+1;
    je = jb+2;

    FF = blkdiag(FF, ...
                 sparse([ ...
                     d1(j), d2(j), d3(j); ...
                     0,     d1(j), d2(j); ...
                     0,     0,     d1(j)]));
    
    AA = blkdiag(AA, A(jb:je, jb:je));
    bb = [bb; 0; 0; 1];
  end
  
  if imag_n
    O2 = zeros(2);
    for j = 1:2:imag_n
      jb = real_end + (ceil(j/2)-1)*8+1;
      je = jb+5;
      
      D1 = [d1(real_n+j), -d1(real_n+j+1); d1(real_n+j+1), d1(real_n+j)];      
      D2 = [d2(real_n+j), -d2(real_n+j+1); d2(real_n+j+1), d2(real_n+j)];
      D3 = [d3(real_n+j), -d3(real_n+j+1); d3(real_n+j+1), d3(real_n+j)];
      
      FF = blkdiag(FF, ...
                   sparse([ ...
                       D1, D2, D3; ...
                       O2, D1, D2; ...
                       O2, O2, D1]));
      
      AA = blkdiag(AA, A(jb:je, jb:je));
      bb = [bb; 0; 0; 0; 0; 1; 0];
    end    
  end
  
  t = type(obj);
  param.k     = t(1)-t(2)-1;
  param.maxit = 0;
  param.tol   = 0;
  param.real  = isreal(obj);
  param.reduction = 0;
  
  [~, ratfun] = rkfit(FF, AA, bb, [xi; xi], param);

end


function [A, b] = internal_build_jordan_form(lambda, order, flag)
% INTERNAL_BUILD_JORDAN_FORM    Constructs a Jordan matrix from given
%                               eigenvalues and algebraic multiplicity. 
% 
% Examples:
%
% [Ar, br] = internal_build_jordan_form([8, 4+3i, 4-3i], 2);
% Ar =
%  8  1            
%     8            
%        4+3i  1      
%              4+3i      
%                    4-3i  1
%                          4-3i  
% br' =  
%  0  1  0  1  0  1
%
% [Ar, br] = internal_build_jordan_form([8, 4+3i, 4-3i], 2, 'real');
% Ar =
%  8  1        
%     8        
%        4 -3  1  
%        3  4     1
%              4 -3
%              3  4
% br' =  
%  0  1  0  0  s  0
% 
% (Here s = 1, although it 'should' be sqrt(2).)
  
  if nargin == 2, flag = 'complex'; end
  
  A = sparse([]); b = [];
  
  if strcmp(flag, 'complex')
    for j = 1:length(lambda)
      A = blkdiag(A, ...
                  sparse(lambda(j)*eye(order) + ...
                         diag(ones(order-1, 1), 1)));
      b = [b; zeros(order-1, 1); 1];
    end
    return
  end
  
  % else: flag == 'real'
  
  lambda_r = lambda(imag(lambda) == 0);
  lambda_i = lambda(imag(lambda) ~= 0);
  
  [A, b] = internal_build_jordan_form(lambda_r, order);

  [lambda_i, success] = util_cplxpair(lambda_i);
  
  assert(success, ['internal_build_jordan_form: Real Jordan' ...
                   'form cannot be build.']);

  s = 1;
  lambda_i = lambda_i(1:2:end);
  
  for j = 1:length(lambda_i)
    a = lambda_i(j);
    a = [real(a) imag(a); -imag(a) real(a)];
    
    % A = blkdiag(A, ...
    %            sparse(kron(eye(order), a) + ...
    %                   kron(diag(ones(order-1, 1), 1), eye(2))));
    
    kron_A = zeros(order*2, order*2);
    for kron_j = 1:order      
      jb = (kron_j-1)*2+1;
      je = jb+1;
      kron_A(jb:je, jb:je) = a;
    end
    for kron_j = 1:order-1
      jb = (kron_j-1)*2+1;
      je = jb+1;
      kron_A(jb:je, (jb:je)+2) = eye(2);
    end
    A = blkdiag(A, sparse(kron_A));
    
    b = [b; zeros(2*order-2, 1); s; 0];
  end  
end
