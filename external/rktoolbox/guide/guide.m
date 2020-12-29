%% RKToolbox Guide
% <html><a href="http://www.maths.manchester.ac.uk/~berljafa/">Mario Berljafa</a>
% and <a href="http://guettel.com/">Stefan G&uuml;ttel</a>
% (The University of Manchester, UK)</html>

%% Overview
% Thank you for your interest in the Rational Krylov Toolbox (RKToolbox).
% The RKToolbox is a collection of scientific computing tools based on
% rational Krylov techniques. The development started in 2013 and the 
% current version 2.5 provides
%
% * an implementation of Ruhe's rational Krylov sequence method [6, 7], 
% allowing to control various options, including user-defined 
% inner products, exploitation of complex-conjugate shifts,
% orthogonalization, rerunning [4], and parallelism [5], 
% * algorithms for the implicit and explicit relocation of the poles of 
% a rational Krylov space [3],
% * a collection of utility functions, e.g., for solving nonlinear 
% eigenvalue problems,
% * an implementation of RKFIT [3, 4], a robust algorithm for rational 
% L2 approximation, including automated degree reduction, and
% * the RKFUN class [4] allowing for numerical computations with rational 
% functions, including support for MATLAB Variable Precision Arithmetic 
% and the Advanpix Multiple Precision toolbox [1].
%
% This guide explains the main functionalities of the toolbox. 
% To run the embedded MATLAB codes the RKToolbox needs to be in
% MATLAB's search path. For details about the installation we refer to
% the *Download* section on <http://rktoolbox.org/>.
%
% <html><a href="../latex/guide.pdf">Click here to view the PDF version of this guide.</a></html>

%% Rational Krylov spaces
% A rational Krylov space is a linear vector space of rational functions in a
% matrix times a vector. Let $A$ be a square matrix of size $N\times N$,
% $\mathbf{b}$ an $N\times 1$ starting vector, and let
% $\xi_1,\xi_2,\ldots,\xi_m$ be a  sequence of
% complex or infinite _poles_ all distinct from the eigenvalues of $A$.
% Then the rational Krylov space of order $m+1$ associated with
% $A,\mathbf{b},\xi_j$ is defined as
%
% $$\displaystyle \mathcal{Q}_{m+1}(A,\mathbf{b}, q_m) = q_m(A)^{-1} \mathrm{span} \{ \mathbf{b},A\mathbf{b},\ldots,A^m \mathbf{b}\},$$
%
% where $q_m(z) = \prod_{{j=1,\xi_j\neq \infty}}^m (z - \xi_j)$ is the
% common denominator of the rational functions associated with the
% rational Krylov space. The
% rational Krylov method by Ruhe [6, 7] computes an orthonormal
% basis $V_{m+1}$ of $\mathcal{Q}_{m+1}(A,\mathbf{b},q_m)$. The basis
% matrix $V_{m+1}$ satisfies a rational Arnoldi decomposition of the form
%
% $$\displaystyle  A V_{m+1} \underline{K_m} = V_{m+1} \underline{H_m}, $$
%
% where $(\underline{H_m},\underline{K_m})$ is an (unreduced) upper
% Hessenberg pencil of size $(m+1)\times m$.
%
% RKT_BIGBREAK
%
% Rational Arnoldi decompositions are useful for several purposes.
% For example, the eigenvalues of the upper $m\times m$ part of the pencil
% $(\underline{H_m},\underline{K_m})$ can be excellent approximations to some
% of $A$'s eigenvalues [6, 7]. Other applications include
% matrix function approximation and rational quadrature, model
% order reduction, matrix equations, nonlinear eigenproblems, 
% and rational least squares fitting (see below).

%% Computing rational Krylov bases
% *Relevant functions:* |rat_krylov, util_cplxpair|
%
% RKT_BIGBREAK
%
% Let us compute $V_{m+1}$, $\underline{K_m}$, and
% $\underline{H_m}$ using the |rat_krylov| function, and verify that the
% outputs satisfy the rational Arnoldi decomposition by computing the
% relative residual norm
% $\| A V_{m+1} \underline{K_m} - V_{m+1} \underline{H_m}\|_2 / \| \underline{H_m}\|_2$.
% For $A$ we take the |tridiag| matrix of size
% $100$ from MATLAB's |gallery|, and $\mathbf{b} = [1,0,\ldots,0]^T$. The
% $m=5$ poles $\xi_j$ are, in order, $-1,\infty, -\mathrm{i}, 0, \mathrm{i}$.

N  = 100;                            % matrix size
A  = gallery('tridiag', N);
b  = eye(N, 1);                      % starting vector
xi = [-1, inf, -1i, 0, 1i];          % m = 5 poles
[V, K, H] = rat_krylov(A, b, xi);
resnorm = norm(A*V*K - V*H)/norm(H)  % residual check

%%
% As some of the poles $\xi_j$ in this example are complex, the matrices
% $V_{m+1}$, $\underline{K_m}$, and $\underline{H_m}$ are complex, too:

disp([isreal(V), isreal(K), isreal(H)])

%%
% However, the poles $\xi_j$ can be reordered so that complex-conjugate
% pairs appear next to each other using the function
% |util_cplxpair|. After reordering the poles, we can call the
% function |rat_krylov| with the |'real'| option, thereby
% computing a real-valued rational Arnoldi decomposition [6].

% Group together poles appearing in complex-conjugate pairs.
xi = util_cplxpair(xi);
[V, K, H] = rat_krylov(A, b, xi, 'real');
resnorm = norm(A*V*K - V*H)/norm(H)
disp([isreal(V), isreal(K), isreal(H)])

%%
% Our implementation |rat_krylov| supports many features not shown
% in the basic description above.
%
% * It is possible the use matrix pencils $(A,B)$ instead of a
%   single matrix $A$. This leads to decompositions of the form
%   $A V_{m+1} \underline{K_m} = BV_{m+1} \underline{H_m}$.
% * Both the matrix $A$ and the pencil $(A,B)$ can be passed
%   either explicitly, or implicitly by providing function handles
%   to perform matrix-vector products and to solve shifted
%   linear systems.
% * Non-standard inner products for constructing the orthonormal
%   bases are supported. 
% * One can choose between CGS and MGS with or without
%   reorthogonalization.
% * Iterative refinement for the linear system solves is supported.
%
% For more details type |help rat_krylov|.


%% Moving poles of a rational Krylov space
% *Relevant functions:* |move_poles_expl, move_poles_impl|
%
% RKT_BIGBREAK
%
% There is a direct link between the starting vector $\mathbf{b}$ and the
% poles $\xi_j$ of a rational Krylov space $\mathcal{Q}_{m+1}$.
% A change of the poles $\xi_j$ to $\breve \xi_j$ can be interpreted as
% a change of the starting vector from $\mathbf{b}$ to $\mathbf{\breve b}$,
% and vice versa. Algorithms for moving the poles of a rational Krylov
% space are described in [3] and implemented in the functions
% |move_poles_expl| and |move_poles_impl|.
%
% RKT_BIGBREAK
%
% *Example:* Let us move the $m=5$ poles $-1,\infty, -\mathrm{i},
% 0,$ and $\mathrm{i}$ into $\breve\xi_j = -j$, $j=1,2,\ldots,5$.

N  = 100;
A  = gallery('tridiag', N);
b  = eye(N, 1);
xi = [-1, inf, -1i, 0, 1i];
[V, K, H] = rat_krylov(A, b, xi);
xi_new = -1:-1:-5;
[KT, HT, QT, ZT] = move_poles_expl(K, H, xi_new);

%%
% The poles of a rational Krylov
% space are the eigenvalues of the lower $m\times m$ part
% of the pencil $(\underline{\breve H_m},\underline{\breve K_m})$
% in a rational Arnoldi decomposition
% $A \breve V_{m+1} \underline{\breve K_m} = \breve V_{m+1} \underline{\breve H_m}$
% associated with that space [3]. By transforming a rational Arnoldi
% decomposition we are therefore effectively moving the poles:

VT = V*QT';
resnorm = norm(A*VT*KT - VT*HT)/norm(HT)
moved_poles = util_pencil_poles(HT, KT).'

%% Rational Krylov fitting (RKFIT)
% *Relevant function:* |rkfit|
%
% RKT_BIGBREAK
%
% RKFIT [3, 4] is an iterative
% Krylov-based algorithm for nonlinear rational approximation.
% Given two families of $N\times N$ matrices
% $\{F^{[j]}\}_{j=1}^{\ell}$ and $\{D^{[j]}\}_{j=1}^{\ell}$, an $N\times n$ block of vectors $B$,
% and an $N\times N$ matrix $A$, the algorithm seeks a family of rational
% functions $\lbrace r^{[j]} \rbrace_{j=1}^{\ell}$ of type
% $(m+k, m)$, all sharing a common denominator $q_m$,
% such that the _relative misfit_
%
% $$\displaystyle {\rm misfit} = \sqrt{\frac{{\sum_{j=1}^\ell \| D^{[j]} [ F^{[j]}B  - r^{[j]}(A)B ]  \|_F^2}}{{\sum_{j=1}^\ell \| D^{[j]} F^{[j]} B \|_F^2}}}\to\min$$
%
% is minimal. The matrices $\{D^{[j]}\}_{j=1}^{\ell}$ are optional,
% and if not provided $D^{[j]}=I_N$ is assumed. The algorithm takes
% an initial guess for $q_m$ and iteratively tries to improve it
% by relocating the poles of a rational Krylov space.
%
% RKT_BIGBREAK
%
% We now show on a simple example how to use the |rkfit| function.
% Consider again the tridiagonal matrix $A$ and the vector $\mathbf{b}$
% from above and let $F = A^{1/2}$.

N  = 100;
A  = gallery('tridiag', N);
b  = eye(N, 1);
F  = sqrtm(full(A));
exact = F*b;

%%
% Now let us find a rational
% function $r_m(z)$ of type $(m, m)$ with $m=10$ such that
% $\| F \mathbf{b} - r_m(A)\mathbf{b} \|_2/\|F \mathbf{b}\|_2$ is
% small. The function |rkfit| requires an input vector of $m$
% initial poles and then tries to return an improved set of poles.
% If we had no clue about where to place the initial poles we can easily
% set them all to infinity. In the following we run RKFIT for at most
% $15$ iterations and aim at relative misfit
% $\| F \mathbf{b} - r_m(A)\mathbf{b} \|_2/\|F\mathbf{b}\|_2$
% below $10^{-10}$. We display the error after each iteration.

[xi, ratfun, misfit] = rkfit(F, A, b, ...
                             repmat(inf, 1, 10), ...
                             15, 1e-10, 'real');

disp(misfit)

%%
% The rational function $r_m(A)\mathbf{b}$ of type $(10, 10)$ approximates
% $A^{1/2}\mathbf{b}$ to about $10$ decimal places. A useful output of
% |rkfit| is the RKFUN object |ratfun| representing the rational
% function $r_m$. It can be used, for example, to evaluate $r_m(z)$:
%
% * |ratfun(A,v)| evaluates $r_m(A)\mathbf{v}$ as a matrix function
% times a vector, 
% * |ratfun(A,V)| evaluates $r_m(A)V$ as a matrix function
% times a matrix, e.g., setting $V=I$ as the identity matrix will return
% the full matrix function $r_m(A)$, or 
% * |ratfun(z)| evaluates $r_m(z)$ as a scalar
% function in the complex plane.
%
% RKT_BIGBREAK
%
% Here is a plot of the error $|x^{1/2} - r_m(x)|$ over the
% spectral interval of $A$ (approximately $[0,4]$), together with the
% values at the eigenvalues of $A$:

figure
ee = eig(full(A)).';
xx = sort([logspace(-4.3, 1, 500) , ee]);
loglog(xx,abs(sqrt(xx) - ratfun(xx))); hold on
loglog(ee,abs(sqrt(ee) - ratfun(ee)), 'r.', 'markers', 15)
axis([4e-4, 8, 1e-14, 1e-3]); xlabel('x'); grid on
title('| x^{1/2} - r_m(x) |','interpreter','tex')

%%
% As expected the rational function $r_m(z)$ is a good approximation
% of the square root over $[0,4]$. It is, however, not a uniform approximation
% because we are approximately minimizing the 2-norm error on the eigenvalues of $A$,
% and moreover we are implicitly using a weight function
% given by the components of $\mathbf{b}$ in $A$'s eigenvector basis.
%
% RKT_BIGBREAK

%%
% Additional features of RKFIT are listed below.
%
% * An automated degree reduction procedure [4, Section 4] is implemented; it
%   takes place if a relative misfit below tolerance is achieved,
%   unless deactivated.
% * Nondiagonal rational approximants are supported; can be
%   specified via an additional |param| structure.
% * Utility functions are provided for transforming scalar data
%   appearing in complex-conjugate pairs into real-valued data, as
%   explained in [4, Section 3.5].
%
% For more details type |help rkfit|.
%
% Some of the capabilities of |RKFUN| are shown in the following section.

%% The RKFUN class
% The |rkfun| class is the fundamental data type to represent and work
% with rational functions. It has already been described above how to
% evaluate an |rkfun| object for scalar and matrix arguments by calling
% |ratfun(z)| or |ratfun(A,v)|, respectively. There are more than 20 other
% methods implemented for |rkfun|, and a list of all these can be obtained
% by typing |methods rkfun|. Here we provide a complete list with
% brief descriptions.

%%
%  basis          - Orthonormal rational basis functions of an rkfun.
%  coeffs         - Expansion coefficients of an rkfun.
%  contfrac       - Convert an rkfun into continued fraction form.
%  diff           - Differentiate an rkfun.
%  disp           - Display information about an rkfun.
%  double         - Convert an rkfun into double precision (undo vpa or mp). 
%  ezplot         - Easy-to-use function plotter. 
%  feval          - Evaluate an rkfun at scalar or matrix arguments.
%  hess           - Convert an rkfun pencil to (strict) upper-Hessenberg form.
%  inv            - Invert an rkfun corresponding to a Moebius transform.
%  isreal         - Returns true if an rkfun is real-valued.
%  minus          - Scalar subtraction. 
%  mp             - Convert an rkfun into Advanpix Multiple Precision format.
%  mrdivide       - Scalar division. 
%  mtimes         - Scalar multiplication. 
%  plus           - Scalar addition. 
%  poles          - Return the poles of an rkfun.
%  poly           - Convert an rkfun into a quotient of two polynomials.
%  power          - Integer exponentiation of an rkfun.
%  rdivide        - Division of two rkfuns.
%  residue        - Convert an rkfun into partial fraction form.
%  rkfun          - The rkfun constructor.
%  roots          - Compute the roots of an rkfun. 
%  size           - Returns the size of an rkfun.
%  subsref        - Evaluate an rkfun (calls feval).
%  times          - Multiplication of two rkfuns.
%  type           - Return the type (m+k,m) of an rkfun.
%  uminus         - Unary minus. 
%  uplus          - Unary plus. 
%  vpa            - Convert rkfun into MATLAB's variable precision format. 

%%
% The names of these methods should be self-explanatory. For example,
% |roots(ratfun)| will return the roots of a |ratfun|, and |residue| will
% compute the partial fraction form. Most methods support the use of
% MATLAB's Variable Precision Arithmetic (VPA) or the Advanpix Multiple
% Precision toolbox (MP). So, for example, |contfrac(mp(ratfun))| will
% compute a continued fraction expansion of |ratfun| using multiple
% precision arithmetic. For more details on each of the methods, type 
% |help rkfun.<method name>|. 
%
% The RKFUN gallery provides some predefined rational functions that may be
% useful. A list of the options can be accessed as follows:

help rkfun.gallery

%%
% RKT_BIGBREAK
%
% Another way to create an RKFUN is to make use of MATLAB's symbolic
% engine. For example, |r = rkfun('(x+1)*(x-2)/(x-3)^2')| will create a 
% rational function as expected. Alternatively, one can specify a rational 
% function by its roots and poles (and an optional scaling factor) using 
% the |rkfun.nodes2rkfun| function. For example, 
% |r = rkfun.nodes2rkfun([-1,2],[3,3])| will create the same rational 
% function as above. 

%% References
% [1] Advanpix LLC., _Multiprecision Computing Toolbox for MATLAB,_  
%     version 3.9.4.10481, Tokyo, Japan, 2016. <http://www.advanpix.com/>.
%
% RKT_BIGBREAK
%
% [2] M. Berljafa and S. Guettel. _A Rational Krylov Toolbox for
%     MATLAB,_ MIMS EPrint 2014.56 (<http://eprints.ma.man.ac.uk/2390/>),
%     Manchester Institute for Mathematical Sciences,
%     The University of Manchester, UK, 2014.
%
% RKT_BIGBREAK
%
% [3] M. Berljafa and S. Guettel. _Generalized rational Krylov
%     decompositions with an application to rational approximation,_
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.
%
% RKT_BIGBREAK
%
% [4] M. Berljafa and S. Guettel.
%     _The RKFIT algorithm for nonlinear rational approximation,_ 
%     MIMS EPrint 2015.38 (<http://eprints.ma.man.ac.uk/2309/>), 
%     Manchester Institute for Mathematical Sciences,
%     The University of Manchester, UK, 2015. 
%
% RKT_BIGBREAK
%
% [5] M. Berljafa and S. Guettel,
%     _Parallelization of the rational Arnoldi algorithm,_
%     MIMS EPrint 2016.32 (<http://eprints.ma.man.ac.uk/2503/>),
%     Manchester Institute for Mathematical Sciences,
%     The University of Manchester, UK, 2016.
%
% RKT_BIGBREAK
%
% [6] A. Ruhe. _Rational Krylov: A practical algorithm for large sparse
%     nonsymmetric matrix pencils,_ SIAM J. Sci. Comput., 19(5):1535--1551,
%     1998.
%
% RKT_BIGBREAK
%
% [7] A. Ruhe. _The rational Krylov algorithm for nonsymmetric eigenvalue
%     problems. III: Complex shifts for real matrices,_ BIT,
%     34:165--176, 1994.
