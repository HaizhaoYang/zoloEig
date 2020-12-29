function [KT, HT, QT, ZT] = util_recover_rad(K, H, flag)
% UTIL_RECOVER_RAD    Transform (H, K) into an upper (quasi-)Hessenberg
%                     pencil.
%
% [KT, HT, QT, ZT] = util_recover_rad(K, H, flag) for (n+1)-by-n matrices K
% and H, produces upper (quasi-)Hessenberg matrices KT and HT and unitary
% matrices QT and ZT such that
%
%      KT = QT*K*ZT, HT = QT*H*ZT.
%
% The parameter flag is optional, and can be passed as flag =
% 'real' to obtain upper quasi-Hessenberg matrices KT and HT,
% provided that K, and H are real-valued. 
% 
% This algorithm is described in
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.

  if nargin == 2, flag = 'complex'; end

  HT = H; KT = K;
  
  [H, K, QT, ZT] = qz(HT(2:end, :), KT(2:end, :), flag);
  QT = blkdiag(1, QT);
  HT = [HT(1, :)*ZT; H];
  KT = [KT(1, :)*ZT; K];

end
