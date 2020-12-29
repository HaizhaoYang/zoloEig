%% An overview of the example collection
%  Mario Berljafa \and Stefan Guettel
%  
%  October 2016
%  
%  Tags: Overview

%%
% Welcome to this example collection, which intends to demonstrate some of 
% the features of the MATLAB Rational Krylov Toolbox. Simply use the menu
% on the left-hand side to navigate through the collection. Each example is
% available as a MATLAB |m-file| and in PDF format (see the links in the 
% above header). All examples are also 
% included in the |rktoolbox.zip| file available from the 
% <http://rktoolbox.org RKToolbox> website.
%
% RKT_BIGBREAK
%
% New examples will be added over time and contributions are more than 
% welcome. If you would like to add an example to this collection please
% email your MATLAB file to stefan.guettel@manchester.ac.uk. You can use
% any |m-file| of this collection as a template.
%
% RKT_BIGBREAK
%
% Here is a simple example illustrating the fascinating convergence
% behaviour of rational Ritz values [1]. The matrix $A$ is diagonal with
% $100$ equispaced eigenvalues in the interval $[1,100]$. Using the 
% rational Arnoldi method [2,3] implemented in |rat_krylov|, we compute 
% Ritz values associated with rational Krylov spaces of increasing dimension 
% with poles alternating between
% $0$ and $\infty$. We then visualize the distance of each Ritz value of order
% $j=1,\ldots,99$ to its closest eigenvalue:

N = 100; m = 99; 
A = spdiags((1:m+1)', 0, N, N); 
b = ones(N, 1);
xi = zeros(1, m); xi(1:2:end) = inf;

[V, K, H] = rat_krylov(A, b, xi);

Am = H(1:m, 1:m)/K(1:m, 1:m);
R  = ones(N, m);
for j = 1:m
  ritz = eig(Am(1:j,1:j));
  R(round(ritz),j) = abs(ritz - round(ritz));
end
imagesc(R); colormap(hot(100)); colorbar
xlabel('order j'); ylabel('Ritz values');
title('distance of Ritz value to closest eigenvalue')

%%
% [1] B. Beckermann, S. Guettel, and R. Vandebril. _On the convergence of 
%     rational Ritz values,_ SIAM J. Matrix Anal. Appl., 31(4):1740--1774, 2010. 
%
% RKT_BIGBREAK
%
% [2] A. Ruhe. _Rational Krylov: A practical algorithm for large sparse
%     nonsymmetric matrix pencils,_ SIAM J. Sci. Comput., 19(5):1535--1551,
%     1998.
%
% RKT_BIGBREAK
%
% [3] A. Ruhe. _The rational Krylov algorithm for nonsymmetric eigenvalue
%     problems. III: Complex shifts for real matrices,_ BIT,
%     34(1):165--176, 1994.
