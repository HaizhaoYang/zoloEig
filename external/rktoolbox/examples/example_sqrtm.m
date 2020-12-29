%% Square root of a symmetric matrix
%  Mario Berljafa \and Stefan Guettel
%
%  June 2015
%
%  Tags: RKFIT, square root

%% Introduction
% This is an example of RKFIT [1, 2, 3] being used for
% approximating $\sqrt{A}\mathbf b$, the action of the matrix
% square root onto a vector $\mathbf b$. RKFIT is described in
% [2, 3] and implemented in the RK Toolobox [1]. This code
% reproduces the example from [2, Section 5.2]. We also compare
% RKFIT to the vector fitting code VFIT3. Vector fitting is
% described in [4, 5].

%%
% We first define the matrix $A$, the vector $\mathbf b$, and the matrix $F$ corresponding
% to $\sqrt{A}$. Our aim is then to find a rational function $r(z)$ of type
% $(m,m)$ such that $\| F\mathbf b - r(A)\mathbf b \|_2/\| F\mathbf b\|_2 \rightarrow \min$.

N = 500;
A = gallery('tridiag',N);
b = eye(N,1);
f = @(x) sqrt(x); fm = @(X) sqrtm(full(X));
F = fm(A);
exact = F*b;

[U, D] = eig(full(A)); ee = diag(D).';

%% Running |rkfit|
% In order to run RKFIT we only need to specify the initial poles of $r$.
% In this example we pretend to not know that a reasonable choice
% of initial poles is on the negative real axis (the branch cut of $\sqrt{z}$)
% and choose 16 infinite poles $\xi$.
% As all quantities $F$, $A$, and $\mathbf b$, as well
% as the initial poles are real, it is advised running RKFIT with the
% |'real'| option. This will attempt to produce a rational approximant with
% perfectly complex conjugate (or even real) poles. We perform $10$ RKFIT
% iterations.

m = 16;
[xi, ratfun, misfit] = rkfit(F, A, b, inf(1, m), 10, 'real');

%%
% It turns out that all computed poles are real negative:

disp(reshape(xi, 4, 4))

%%
% Here is a convergence plot of RKFIT, showing the relative misfit
% $\| F\mathbf b - r(A)\mathbf b\|_2/\| F\mathbf b\|_2$ at each iteration.

figure(1)
semilogy(misfit, 'ro--')
xlabel('iteration'), title('relative 2-norm misfit')
legend('RKFIT')

%% Evaluating the rational approximant
% The second output |ratfun| is an object that can be used to
% evaluate the computed rational approximant. This evaluation is implemented
% in two ways. The first option is to evaluate a matrix function $r(A)\mathbf b$ by
% calling |ratfun(A,b)| with two input arguments. For example, here we are
% calculating the absolute misfit or $r$:

disp(norm(F*b - ratfun(A, b)))

%%
% Alternatively, we can evaluate $r(z)$ pointwise by giving only one input
% argument. Let us plot the modulus of the scalar error function
% $\mathrm{err}(x) = f(x) - r(x)$ over the spectral interval of $A$, which is
% approximately $[0,4]$.

figure(2)
xx = sort([logspace(-6,log10(15.728), 5000), ee ]);
errf = f(xx) - ratfun(xx);
loglog(xx,abs(errf), 'm-')
hold on
errf = f(ee) - ratfun(ee);
loglog(ee,abs(errf), 'g+')
xlabel('x'), title('abs(sqrt(x) - r(x))')
legend('RKFIT error','error at evs')

%%
% How often does this error function change sign on $[0,4]$ (or rather on the
% discretisation set |xx|, which was chosen very fine with 5000 points)?
% Equivalently, how many interpolation points are there?

sc = length(find(diff(sign(errf))));
disp(['The error function has ' num2str(sc) ' sign changes on [0,4].'])

%% Some different choices for the initial poles
% Choosing all initial poles equal to infinity seemed to work fine. Of
% course we can achieve convergence in fewer iterations by choosing
% negative initial poles, e.g., logarithmically spaced on the negative real
% axis. Let us rerun RKFIT with the new initial poles and compare
% the convergence.

tmp = -logspace(-8, 8, m);
[xi, ratfun, misfit] = rkfit(F, A, b, tmp, 10, 'real');
figure(1), hold on
semilogy(misfit, 'rs--')
legend('RKFIT (infinite init poles)', ...
       'RKFIT (negative init poles)')

%%
% Now let us try a very bad choice for the initial poles, which is to place
% them right into the spectral interval $[0,4]$ of $A$:

tmp = linspace(0, 4, m);
[xi, ratfun, misfit] = rkfit(F, A, b, tmp, 10, 'real');
semilogy(misfit, 'r*--')
legend('RKFIT (infinite init poles)', ...
       'RKFIT (negative init poles)', ...
       'RKFIT (bad init poles)')

%%
% In this example, RKFIT is not very affected by the choice of the initial poles
% and converges robustly.

%% Comparison with vector fitting
% We now compare RKFIT with the vector fitting code described in
% [4, 5].
%
% RKT_BIGBREAK
%
% The following code makes sure the |vectfit3| implementation of
% VFIT is on MATLAB's search path, and if it is, defines some
% options.

if exist('vectfit3') ~= 2
  warning('Vector fitting VFIT3 not found in the Matlab path.');
  warning('VFIT can be downloaded from:')
  warning('http://www.sintef.no/Projectweb/VECTFIT/Downloads/VFUT3/')
  warning('Skipping comparison with VFIT3.');
  return
end

xi = -logspace(-8, 8, m);
opts.relax  = 1;   % Relaxed non-triviality constraint.
opts.stable = 0;   % Do not enforce stable poles.
opts.asymp  = 2;   % Options are 1, 2, or 3.
opts.skip_pole = 0;
opts.skip_res  = 0;
opts.spy1 = 0;    opts.spy2 = 0;
opts.logx = 0;    opts.logy = 0;
opts.errplot = 0; opts.phaseplot = 0;

%%
% For VFIT3 we need interpolation nodes and weights, and in this
% example these need to be chosen as the eigenvalues of $A$ and the
% components of $\mathbf b$ in the eigenvector basis of $A$, respectively.
%
% RKT_BIGBREAK
%
% *Note:* _The below warnings outputted by VFIT are due to ill-conditioned 
% linear algebra problems being solved, a problem that is circumvented with 
% RKFIT by the use of discrete-orthogonal rational basis functions._

nodes = diag(D); nodes = nodes(:).';
fvals = f(nodes);
weights = U'*b; weights = weights(:).';

for iter = 1:10
  [SER, xi, rmserr, fit] = ...
      vectfit3(fvals, nodes, xi, weights, opts);
  err_vfit(iter) = norm((fvals - fit).*weights)/norm(exact);
end

figure(1)
semilogy(err_vfit,'bs:')
legend('RKFIT (infinite init poles)','RKFIT (negative init poles)',...
       'RKFIT (bad init poles)','VFIT (negative init poles)')
axis([1, 10, 1e-14, 1e-5])

ffit = @(zz) arrayfun(@(z) ...
                      sum((SER.C).'./(z-diag(SER.A))) + ...
                      SER.D + z*SER.E,zz);
figure(2)
loglog(xx,abs(f(xx) - ffit(xx)), 'b:')
axis([1e-5, 15.7, 1e-15, 5e-4])
set(gca, 'XTick', 10.^(-5:1))
legend('RKFIT error', 'error at evs', 'VFIT error')

%%
% And the same with the "bad" initial poles.

xi = linspace(0, 4, m);
for iter = 1:10
  [SER, xi, rmserr, fit] = ...
      vectfit3(fvals, nodes, xi, weights, opts);
  err_vfit(iter) = norm((fvals - fit).*weights)/norm(exact);
end
figure(1)
semilogy(err_vfit, 'b*:')
legend('RKFIT (infinite init poles)', ...
       'RKFIT (negative init poles)', ...
       'RKFIT (bad init poles)',      ...
       'VFIT (negative init poles)',  ...
       'VFIT (bad init poles)')

%%
% This is the end of this example. The following creates a thumbnail.

figure(2), plot(NaN)

%% References
%
% [1] M. Berljafa and S. Guettel. _A Rational Krylov Toolbox for
%     MATLAB,_ MIMS EPrint 2014.56 (<http://eprints.ma.man.ac.uk/2390/>),
%     Manchester Institute for Mathematical Sciences,
%     The University of Manchester, UK, 2014.
%
% RKT_BIGBREAK
%
% [2] M. Berljafa and S. Guettel. _Generalized rational Krylov
%     decompositions with an application to rational
%     approximation,_ SIAM J. Matrix Anal. Appl., 36(2):894--916,
%     2015.
%
% RKT_BIGBREAK
%
% [3] M. Berljafa and S. Guettel.
%     _The RKFIT algorithm for nonlinear rational approximation,_ 
%     MIMS EPrint 2015.38 (<http://eprints.ma.man.ac.uk/2309/>),
%     Manchester Institute for Mathematical Sciences,
%     The University of Manchester, UK, 2015.
%
% RKT_BIGBREAK
%
% [4] B. Gustavsen. _Improving the pole relocating properties of
%     vector fitting,_ IEEE Trans. Power Del., 21(3):1587--1592,
%     2006.
%
% RKT_BIGBREAK
%
% [5] B. Gustavsen and A. Semlyen. _Rational approximation of
%     frequency domain responses by vector fitting,_ IEEE
%     Trans. Power Del., 14(3):1052--1061, 1999.
