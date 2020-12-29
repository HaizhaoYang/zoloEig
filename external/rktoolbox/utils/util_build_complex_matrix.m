function A = util_build_complex_matrix(data)
% UTIL_BUILD_COMPLEX_MATRIX    Constructs a diagonal matrix from 1D
%                              array.
% 
% A = util_build_complex_matrix(data) takes the complex 1D array of
% length n and constructs a complex diagonal matrix of size
% 2*n-by-2*n in sparse format with eigenvalues [data conj(data)].  

  assert(size(data, 1)==1 || size(data, 2) == 1, ...
         'The input array needs to be a row or column vector.')    
  if size(data, 1) > 1, data = data.'; end  
  
  n = length(data);
  N = 2*n;
  
  A = zeros(N, 1);
  A(1:2:end) = data;
  A(2:2:end) = conj(data);
  A = sparse(1:N, 1:N, A, N, N, N);
  
end
