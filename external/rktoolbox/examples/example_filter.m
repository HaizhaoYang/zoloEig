%% Electronic filter design using RKFUN arithmetic
%  Mario Berljafa \and Stefan Guettel
%
%  September 2015
%  
%  Tags: RKFUN, rational functions

%% Introduction
% With Version 2.2 of the RKToolbox the |rkfun| class and its methods 
% have been substantially extended. It is now possible to perform basic arithmetic 
% operations on rational functions that go beyond scalar multiplication
% and addition of constants. For example, one can now add, multiply, 
% divide, or exponentiate |rkfun| objects. Composition of an |rkfun| with a 
% Moebius transform is also possible. 
% Furthermore, the |rkfun| constructor now comes with the option to construct
% an object from a symbolic expression, or a string specifying a rational
% function from the newly introduced |gallery|. We now demonstrate some
% of these features. An application to electronic filter design is given 
% below.

%% The RKFUN gallery
% The |rkfun| class provides a gallery function which allows for the quick 
% construction of some useful rational functions. A list of the
% functions currently implemented can be obtained by typing:

help rkfun.gallery

%%
% For example, we can construct an |rkfun| corresponding to a Moebius
% transform $r(z) = (4z + 3)/(2z - 1)$ as follows:

r = rkfun.gallery('moebius', [4, 3, 2, -1])

%%
% As always with |rkfun| objects, we can perform several computations on |r|, 
% such as computing its roots and poles:

format long e
disp([roots(r) poles(r)])

%% The symbolic constructor
% Symbolic strings can also be used to construct |rkfun| objects, provided 
% that the required MATLAB symbolic toolboxes are installed. Here we
% construct an |rkfun| corresponding to the rational function 
% $s(z) = 3(z-3)(z^2-1)/(z^{5}+1)$:

s = rkfun('3*(z-3)*(z^2-1)/(z^5+1)')

%%
% This works by first finding the roots and poles of the input function
% symbolically and then constructing the |rkfun| with these roots and poles
% numerically. Note that |s| is only of type $(2,4)$ as one of the five 
% denominator roots has been cancelled
% out by symbolic simplifications prior to the numerical construction.
% The constructor should issue an error if the provided string fails to be
% parsed by the symbolic engine, or if it does not represent a rational
% function.

%% Arithmetic operations with rational functions
% Since Version 2.2 of the RKToolbox one can add, multiply, divide, 
% and exponentiate |rkfun| objects. All these operations are implemented via transformations on
% generalized rational Krylov decompositions [1, 2]. 
% For example, here is the product of the two rational functions $r$ and
% $s$ from above:

disp(r.*s)

%%
% Note that it is necessary to use point-wise multiplication |.*|, not 
% the matrix-multiplication operator |*|. This is to be consistent with 
% the MATLAB notation, and also with the notation used in the Chebfun
% system [5] (where the |*| operator returns the inner product of two
% polynomials; something we have not yet implemented for |rkfun| objects). 
%
% RKT_BIGBREAK
%
% It is also possible to compose |rkfun| objects as long as the inner function
% is a Moebius transform, i.e., a rational function of type at most
% $(1,1)$. The rational function $r$ from above _is_ a Moebius transform,
% hence we can form the function $f(z) = s(r(z))^{-1}$:

f = 1./s(r)

%%
% The roots of $f$ should correspond to the poles of $s(r)$, which are the 
% poles of $s$ mapped under the inverse function $r^{-1}$. The latter 
% inverse function is indeed well defined in the whole complex plane 
% as $r$ is an invertible Moebius transform. We can compute it via the 
% |inv| command. Let's verify that the roots/poles agree numerically:

rts = sort(roots(f));
pls = sort(feval(inv(r), poles(s)));
disp([rts, pls])

%% Constructing rational filters
% Rational filter functions are ubiquitous in scientific computing and
% engineering. For example, in signal processing [3, 6] one is interested in 
% deriving rational functions that act as filters on selected frequency 
% bands. Using  the |gallery| of the RKToolbox and some |rkfun| transformations, 
% we can construct meaningful rational filters in just a few lines of 
% MATLAB code. Below is a plot of four popular
% filter types, which are obtained by rational transforms of Chebyshev
% polynomials or the |step| function from |rkfun|'s |gallery|.

x = rkfun;
butterw = 1./(1 + x.^16);
cheby   = rkfun('cheby',8);
cheby1  = 1./(1 + 0.1*cheby.^2);
cheby2  = 1./(1 + 1./(0.1*cheby(1./x).^2));
ellip   = rkfun('step');

figure
subplot(221),ezplot(butterw,[0,2]);    title('Butterworth')
subplot(222),ezplot(cheby1, [0,2],'r');title('Chebyshev type 1')
subplot(223),ezplot(cheby2, [0,2],'g');title('Chebyshev type 2')
subplot(224),ezplot(ellip,  [0,2],'m');title('Elliptic')

%%
% The reader may compare this plot with that at the bottom of the Wikipedia
% page on _electronic filters_ [7].
% For example, the filter |cheby2| involves multiple inverses of a Chebyshev
% polynomial in the transformed variable $x \mapsto x^{-1}$. It has the
% so-called _equiripple property_ in the stopband, which is the region where
% the filter value is close to zero. 
% The elliptic filter |ellip|, also known as Cauer filter [4], 
% has equiripple properties in both the stop- and passbands. Such 
% filters are based on Zolotarev's equioscillating rational functions [8],
% which are also implemented in the |gallery| of the RKToolbox.

%% Limitations
% Although we hope that the new |rkfun| capabilities demonstrated above
% are already useful for many practical purposes, there are still 
% some short-comings one has to be aware of. The main problem is that
% combinations of |rkfun| objects may have degrees higher than theoretically
% necessary, which may lead to an unnecessarily fast growth of parameters.
% For example, when subtracting an |rkfun| from itself,

r = rkfun('(1-x)/(1+x)');
d = r - r

%%
% we currently obtain a type $(2,2)$ rational function, instead of the
% expected type $(0,0)$. This is because the sum (or difference) of two type
% $(1,1)$ rational functions is of type $(2,2)$ _in the worst case_, and we
% currently do not perform any degree reduction on the sum (or difference). 
% (A similar
% problem is encountered with multiplication or division.) The roots and
% poles of the function |d| are meaningless, but the evaluation works
% fine, except when we evaluate at (or nearby) a "legacy" pole:

d([-1, -1+1e-16, 0, 1, 2])

%%
% A numerical degree reduction would probably require the concept of a
% "domain of evaluation." To illustrate, consider the rational function
% $r(z) = \frac{10^{-14}}{z} + \frac{1}{z-1}$. The exact type of this
% function is $(1,2)$, but the residue associated with $z=0$ is tiny, so 
% one may conclude that this pole could be removed.
% However, when evaluating $r$ for $z\approx 0$, the removal of this pole 
% would lead to an inaccurate result (try $z = 10^{-14}$). Hence, reliable
% degree reduction may only be possible when a relevant "domain of 
% evaluation" is specified.

%% References
%
% [1] M. Berljafa and S. Guettel. _Generalized rational Krylov
%     decompositions with an application to rational
%     approximation,_  SIAM J. Matrix Anal. Appl., 36(2):894--916,
%     2015.
%
% RKT_BIGBREAK
%
% [2] M. Berljafa and S. Guettel. 
%     _The RKFIT algorithm for nonlinear rational approximation,_ 
%     MIMS EPrint 2015.38 (<http://eprints.ma.man.ac.uk/2309/>), 
%     Manchester Institute for Mathematical Sciences, 
%     The University of Manchester, UK, 2015.
%
% RKT_BIGBREAK
%
% [3] H. Blinchikoff and A. Zverev. _Filtering in the Time and Frequency 
%     Domains,_ John Wiley & Sons Inc., New York, 1976.
%
% RKT_BIGBREAK
%
% [4] W. Cauer. _Ein Interpolationsproblem mit Funktionen mit positivem 
%     Realteil,_ Math. Z., 38(1):1-44, 1934.
%
% RKT_BIGBREAK
%
% [5] T. A. Driscoll, N. Hale, and L. N. Trefethen. _Chebfun Guide,_  
%     Pafnuty Publications, Oxford, 2014. <http://www.chebfun.org>
%
% RKT_BIGBREAK
% 
% [6] M. Van Valkenburg. _Analog Filter Design,_ Holt, Rinehart and Winston, 
%     1982.
%
% RKT_BIGBREAK
%
% [7] Wikipedia, entry _Electronic filter_ as of 17/09/2015.  
%     <https://goo.gl/vkQwWG>
%
% RKT_BIGBREAK
%
% [8] E. I. Zolotarev. _Application of elliptic functions to questions 
%     of functions deviating least and most from zero,_  
%     Zap. Imp. Akad. Nauk St. Petersburg, 30:1-59, 1877. In Russian.
