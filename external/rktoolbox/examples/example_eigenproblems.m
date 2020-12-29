%% Solving large sparse eigenvalue problems
%  Mario Berljafa \and Stefan Guettel
%
%  June 2015
%
%  Tags: RAT_KRYLOV, eigenvalues

%% Introduction
% The first use of rational Krylov methods was for the solution of large 
% sparse eigenvalue problems $A\mathbf x = \lambda B \mathbf x$, where  
% $A$ and $B$ are $N\times N$ matrices and $(\lambda,\mathbf x)$ are the
% wanted eigenpairs; see [3, 4, 5, 6]. 
% Let $\mathbf{b}$ be an $N\times
% 1$ vector and $m$ a positive integer. 
% The _rational Krylov space_
% $\mathcal{Q}_{m+1}(A,\mathbf{b}, q_m)$ is defined as
% $\mathcal{Q}_{m+1}(A,\mathbf{b}, q_m) =
% q_m(A)^{-1}
% \mathrm{span} \{ \mathbf{b},A\mathbf{b},\ldots,A^m \mathbf{b}\}.$
% Here, $q_m$ is a polynomial of degree at most $m$ having roots 
% $\xi_1,\ldots,\xi_m$, called _poles_ of the rational Krylov space. 
% If $q_m$ is a constant nonzero polynomial, then 
% $\mathcal{Q}_{m+1}(A,\mathbf{b}, q_m)$ is a standard (polynomial)
% Krylov space.
%
% RKT_BIGBREAK
%
% The rational Arnoldi method [4, 5] can be used to compute an orthonormal basis
% $V_{m+1}$ of $\mathcal{Q}_{m+1}$ such that a rational Arnoldi
% decomposition $A V_{m+1}\underline{K_m} = B V_{m+1} \underline{H_m}$ is
% satisfied, where $\underline{H_m}$ and $\underline{K_m}$ are upper
% Hessenberg matrices of size $(m+1)\times m$. 
%
% RKT_BIGBREAK
%
% In the following we initialize two matrices for the pencil $(A,
% B)$ and plot the full spectrum $\Lambda(A, B)$. For simplicity we
% take $B=I$ and choose a rather small size to be able to compute 
% all eigenvalues exactly.

load west0479
A = west0479;   
B = speye(479); 
ee = eig(full(A), full(B)); % Cannot do this for larger matrices!
figure, plot(ee, 'ko', 'linewidth', 2), hold on
title('eigenvalues of (A,B)')
legend('exact eigenvalues')

%%
% We now construct a rational Krylov space with $m=10$ 
% poles all set to zero, i.e., we are interested in the
% generalized eigenvalues of $(A,B)$ closest to zero. The inner
% product for the rational Krylov space can be defined by the user,
% or otherwise is the standard Euclidean one. (Since $B=I$ in this example,
% the two coincide.)

rng(0);
b  = randn(479,1);
m  = 10;
xi = repmat(0, 1, m);
param.inner_product = @(x,y) y'*B*x;
[V, K, H] = rat_krylov(A, B, b, xi, param);

warning off, nrmA = normest(A); nrmB = normest(B); warning on

%%
% We can easily check the validity of the rational Arnoldi decomposition by
% veryfing that the residual norm is close to machine precision:

disp(norm(A*V*K - B*V*H)/(nrmA + nrmB))

%%
% The basis $V_{m+1}$ is close to orthonormal too:

disp(norm(param.inner_product(V, V)-eye(m+1)))

%% Extracting approximate eigenpairs
% A common approach for extracting approximate eigenpairs from a
% search space is by using _Ritz approximations_ or variants thereof.
% Let $C$ be an $N\times N$ matrix, and $X$ an $N\times m$ matrix. 
% The pair $(\vartheta, \mathbf{y}\equiv X\mathbf{z})$ is called a
% 
% * _Ritz pair_ for $C$ with respect to $\mathcal{R}(X)$ if
% $C\mathbf{y}-\vartheta \mathbf{y} \perp \mathcal{R}(X)$;
% * _harmonic Ritz pair_ for $C$ with respect to $\mathcal{R}(X)$ if
% $C\mathbf{y}-\vartheta \mathbf{y} \perp \mathcal{R}(CX)$.
%
% Assume that $B$ is nonsingular. It follows easily (see [1, Lemma
% 2.4] and [2, Theorem 2.1]) from $A V_{m+1}\underline{K_m}
% = B V_{m+1} \underline{H_m}$ and the definition of (harmonic)
% Ritz pairs given above that
%
% * Ritz pairs $(\vartheta, \mathbf{y}\equiv V_{m+1}\underline{K_m} \mathbf{z})$ for
% $B^{-1}A$ with respect to $\mathcal{R}(V_{m+1}\underline{K_m})$
% arise from solutions of the generalized eigenvalue problem 
% $\underline{K_m}^*\underline{H_m}^{\phantom{.}}\mathbf{z} = \vartheta
% \underline{K_m}^*\underline{K_m}^{\phantom{.}}\mathbf{z}$. Since
% $\underline{K_m}^{\phantom{.}}$ is of full rank,
% $\underline{K_m}^*\underline{K_m}^{\phantom{.}}$ is nonsingular,
% and we can equivalently solve the standard eigenvalue problem
% $\underline{K_m}^\dagger\underline{H_m}^{\phantom{.}}\mathbf{z} = \vartheta
% \mathbf{z}$;
% * harmonic Ritz pairs $(\vartheta, \mathbf{y}\equiv V_{m+1}\underline{K_m} \mathbf{z})$ for
% $B^{-1}A$ with respect to $\mathcal{R}(V_{m+1}\underline{K_m})$
% arise from solutions of the generalized eigenvalue problem 
% $\underline{H_m}^*\underline{H_m}^{\phantom{.}}\mathbf{z} = \vartheta
% \underline{H_m}^*\underline{K_m}^{\phantom{.}}\mathbf{z}$. Since
% $\underline{H_m}^{\phantom{.}}$ is of full rank,
% $\underline{H_m}^*\underline{H_m}^{\phantom{.}}$ is nonsingular,
% and we can equivalently solve the standard eigenvalue problem
% $\underline{H_m}^\dagger\underline{K_m}^{\phantom{.}}\mathbf{z} = \lambda
% \mathbf{z}$, and take $\vartheta = \frac{1}{\lambda}$.

[Xr, Dr] = eig(K\H);
[Xh, Dh] = eig(H\K);
ritz = diag(Dr);
hrm_ritz = 1./diag(Dh);

plot(real(ritz), imag(ritz),'bx', 'linewidth', 3)
plot(real(hrm_ritz), imag(hrm_ritz),'m+', 'linewidth', 2)
axis([-0.02, 0.02, -0.1, 0.1])
legend('exact eigenvalues', ...
       'Ritz approximations', ...
       'harmonic Ritz')

%% Accuracy of the approximate eigenpairs
% We can evaluate the accuracy of the (harmonic) Ritz pairs from
% the relative residual norm $\frac{\|A\mathbf{y}-\vartheta B
% \mathbf{y}\|_2}{(\|A\|_2+|\vartheta|\|B\|_2)\|\mathbf{y}\|_2}$. From
% the rational Arnoldi decomposition $A V_{m+1}\underline{K_m} = B
% V_{m+1} \underline{H_m}$ we have $A\mathbf{y}-\vartheta
% B\mathbf{y} = AV_{m+1}\underline{K_m} \mathbf{z} -
% \vartheta BV_{m+1}\underline{K_m} \mathbf{z} =
% BV_{m+1}(\underline{H_m}-\vartheta\underline{K_m})\mathbf{z}$. Hence,
% a cheap estimate of the accuracy of an approximate eigenpair is
% the norm
% $\|(\underline{H_m}-\vartheta\underline{K_m})\mathbf{z}\|_2$. If
% this norm is small compared to $\|B^{-1}A\|_2$, we have computed 
% an eigenpair of a nearby problem. 
% It seems that in this example two eigenpairs have already converged 
% to very high accuracy:

approx_residue = @(X) arrayfun(@(i)norm(X(:, i)), 1:size(X, 2));

approx_res = [approx_residue(H*Xr-K*Xr*Dr).' ...              
              approx_residue(H*Xh-K*Xh*diag(hrm_ritz)).']; 
disp(approx_res)



%% Expanding the rational Arnoldi decomposition
% Let us perform $8$
% further iterations with |rat_krylov|, with $2$ repeated poles being the
% harmonic Ritz eigenvalues expected to converge next. Since the poles
% appear in complex-conjugate pairs, we can turn on the |real| flag
% for |rat_krylov| and end up with a real-valued quasi rational
% Arnoldi decomposition [6].

[~, ind] = sort(approx_res(:, 2));
xi = repmat(hrm_ritz(ind([3, 4]))', 1, 4);
param.real = 1;

[V, K, H] = rat_krylov(A, B, V, K, H, xi, param);

%%
% Let us check the residual norm of the extended rational Arnoldi
% decomposition, and verify that the decomposition has the original
% $10$ poles at zero, and the newly selected $8$ poles.

disp(norm(A*V*K - B*V*H)/(nrmA + nrmB))
disp(util_pencil_poles(K, H).')

%%
% Finally, the (improved) $18$ Ritz pairs are evaluated, both standard
% and harmonic.

[Xr, Dr] = eig(K\H);
[Xh, Dh] = eig(H\K);

ritz = diag(Dr);
hrm_ritz = 1./diag(Dh);

approx_res = [approx_residue(H*Xr-K*Xr*Dr).' ...              
              approx_residue(H*Xh-K*Xh*diag(hrm_ritz)).']; 

disp(approx_res)

%%
%

figure, plot(ee, 'ko', 'linewidth', 2), hold on
plot(real(ritz), imag(ritz),'bx', 'linewidth', 3)
plot(real(hrm_ritz), imag(hrm_ritz),'m+', 'linewidth', 2)
axis([-0.08, 0.08, -0.1, 0.1])
legend('exact eigenvalues', ...
       'Ritz approximations', ...
       'harmonic Ritz')

%%
% Harmonic Ritz pairs are typically better than (standard)
% Ritz pairs for interior eigenvalues, thought this is not yet
% fully understood. Also, for symmetric matrices the two sets of Ritz 
% values interlace each other, and hence their distance is not large  
% as ultimately both sets converge to the same eigenvalues.


%% References
% [1] G. De Samblanx and A. Bultheel, _Using implicitly filtered RKS
%     for generalised eigenvalue problems,_ J. Comput. Appl. Math.,
%     107(2):195--218, 1999.
%
% RKT_BIGBREAK
%
% [2] R. B. Lehoucq and K. Meerbergen. _Using generalized Cayley
%     transformations within an inexact rational Krylov sequence
%     method,_ SIAM J. Matrix Anal. Appl., 20(1):131--148, 1998.
%
% RKT_BIGBREAK
%
% [3] A. Ruhe. _Rational Krylov sequence methods for eigenvalue computation,_
%     Linear Algebra Appl., 58:391--405, 1984.
%
% RKT_BIGBREAK
%
% [4] A. Ruhe. _Rational Krylov algorithms for nonsymmetric eigenvalue 
%     problems,_ Recent Advances in Iterative Methods. Springer New York, 
%     pp. 149--164, 1994.
%
% RKT_BIGBREAK
%
% [5] A. Ruhe. _Rational Krylov: A practical algorithm for large sparse
%     nonsymmetric matrix pencils,_ SIAM J. Sci. Comput.,
%     19(5):1535--1551, 1998.
%
% RKT_BIGBREAK
%
% [6] A. Ruhe. _The rational Krylov algorithm for nonsymmetric eigenvalue
%     problems. III: Complex shifts for real matrices,_ BIT,
%     34(1):165--176, 1994.
