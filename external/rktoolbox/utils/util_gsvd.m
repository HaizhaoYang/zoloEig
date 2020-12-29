function [U, W, X, C, S, R] = util_gsvd(K, H)
% [U, W, X, C, S, R] = util_gsvd(K, H, mode) computes the GSVD of the
% pencil (K, H) calling MATLAB's GSVD command and return also the R
% factor of the thin QR of [K; H].
  
  [U, W, X, C, S] = gsvd(K, H);  
  [Q, R] = qr([K; H], 0);
    
end
