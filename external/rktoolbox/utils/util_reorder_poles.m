function [KT, HT, QT, ZT] = util_reorder_poles(K, H, perm, flag)
% UTIL_REORDER_POLES    Reorders the poles of (H, K).
%
% [KT, HT, QT, ZT] = util_reorder_poles(K, H, perm, flag) reorders the poles
% of the upper (quasi-)Hessenberg pencil (H, K) of size (n+1)-by-n. That is,
% it produces upper (quasi-)Hessenberg matrices KT and HT and unitary
% matrices QT and ZT such that
%
%      KT = QT*K*ZT, HT = QT*H*ZT,
%
% with the poles (HT, KT) being reordered as specified by perm. The j-th pole
% in (H, K) shall be perm(j)-th pole in (HT, KT). (The vector -perm is passed
% to ordqz as the cluster parameter.) The parameter flag is optional, and can
% be passed as flag = 'real' to obtain upper quasi-Hessenberg matrices KT and
% HT, provided that K, and H are real-valued.     

  if nargin == 3, flag = 'complex'; end
  
  [HT, KT, QT, ZT] = ordqz(H(2:end, :), K(2:end, :), ...
                           eye(size(K, 2)), eye(size(K, 2)), ...			   
			   -perm);
  
  KT = [K(1, :)*ZT; KT];
  HT = [H(1, :)*ZT; HT];
  
  QT = blkdiag(1, QT);
  
end
