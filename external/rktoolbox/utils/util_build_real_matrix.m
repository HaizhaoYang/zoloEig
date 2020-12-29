function A = util_build_real_matrix(data)
% UTIL_BUILD_REAL_MATRIX    Constructs a block-diagonal matrix from
%                           1D array.
% 
% A = util_build_real_matrix(data) takes the complex 1D array of
% length n and constructs a real block diagonal matrix of size
% 2*n-by-2*n in sparse format with eigenvalues [data conj(data)].
%
% For instance
% 
% util_build_real_matrix([1+2i 4i 3 5-6i])
% 
% constructs a sparse matrix with the following nonzero structure:  
%
%    1  -2     
%    2   1     
%                -4     
%             4         
%                      3         
%                          3     
%                              5   6
%                             -6   5

  assert(size(data, 1)==1 || size(data, 2) == 1, ...
         'The input array needs to be a row or column vector.')    
  if size(data, 1) > 1, data = data.'; end  
  
  n = length(data);

  % Real part on the diagonal.
  ii = 1:2*n;
  jj = ii;
  ss = kron(real(data), [1 1]);
    
  % Imaginry on the off-diagonal;
  ii = [ii 2*(1:n) 2*(1:n)-1];
  jj = [jj ii(3*n+1:end) ii(2*n+1:3*n)];
  ss = [ss imag(data)]; 
  ss = [ss -ss(2*n+1:end)];
  
  ii = ii(ss ~= 0);
  jj = jj(ss ~= 0);
  ss = ss(ss ~= 0);
  
  A = sparse(ii, jj, ss, 2*n, 2*n, length(ss));
    
end
