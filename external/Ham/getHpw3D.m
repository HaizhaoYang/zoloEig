function [Hfun,That,Vr,Tfun,Vfun,numEigVal] = getHpw2D(npc,K,opt)
%
% Input:
% npc - number of grid points per cell
% K - number of Gaussian wells per dimension
% opt - contains opt.isSC, opt.numDL, opt. sigma, specifying the Gaussian
%       wells in each cell
% opt.isReal - 0: Hfun is the Hamiltonian discretized in the Fourier 
%                 domain, i.e. Hfun is a functional acting on complex 
%                 functions defined in the Fourier domain
%              1: Hfun is a hamiltonian discretized in the real space
%
% Output:
% Hfun -   A function handle representing the Hamiltonian, defined in the Fourier
%          domain or real space, depending on isReal
% That -   The Laplacian operator defined in the Fourier domain
% Vr   -   The potential function defined in the real space
% Tfun -   A function handle representing the Laplacian operator, defined in the Fourier
%          domain or real space, depending on isReal
% Vfun -   A function handle representing the potential term, defined in the Fourier
%          domain or real space, depending on isReal
% numEigVal - the number of eigenvalues, the number of orbitals


%% generate and dicretize the Hamiltonian using pseudospectral method
% If isReal=0, then
% i.e. H = T + F'*V*F,
% where T is the discreitzation of -0.5*Laplacian in the Fourier domain,
% i.e., a diagonal matrix; V is the potential defined in the space domain;
% F is the matrix consists of Fourier basis as columns; F'*v=fftn(v), the
% normalized fft of v.
%
% If isReal=1, then
% i.e. H = F*T*F' + V,
% where T is the discreitzation of -0.5*Laplacian in the Fourier domain,
% i.e., a diagonal matrix; V is the potential defined in the space domain;
% F is the matrix consists of Fourier basis as columns; F'*v=fftn(v), the
% normalized fft of v.

% The Hamiltonian constructed in this code is always positive definite.
if nargin<3
    opt.depth = 1;
    opt.numDL = 0;
    opt.sigma = 0.2;
    opt.isReal = 1;
end
sigma = opt.sigma;
numDL = opt.numDL;
depth = opt.depth;
isReal = opt.isReal;

numEigVal = (K)^2; 
c = 0.5;
h = 1/npc;%step size in space
n = npc*K;% size of grid in 1D space

ks = [0:n/2-1 -n/2:-1]';%grid points in Fourier domain
[k1,k2] = ndgrid(ks);

% construct the Laplacian operator
That = 2*pi^2*(k1.^2 + k2.^2); % diagonal elements of -0.5*Laplacian in Fourier domain, size n^2 in 2D
if (isReal)
    Tfun = @(x) real( reshape( ifft2(reshape(repmat(That,1,numel(x)/n^2),n,n,numel(x)/n^2).*fft2(reshape(x,n,n,numel(x)/n^2))) ,n^2,numel(x)/n^2) );
    % Trfun(x)is -0.5*Laplacian*f in the real space for real functions x
else
    Tfun = @(x) reshape(repmat(That,1,numel(x)/n^2),n^2,numel(x)/n^2).*x;
    % Tfun(x)is -0.5*Laplacian*f in the Fourier domain for complex functions, when x=fft2n(f)
end

% construct the potential energy operator
x=h*(0:npc-1);
v = gaussmf(x,[sigma c]);
v = v(:);
vv = -v*v'*depth;
Vr=repmat(vv,K,K);
if (numDL>0)
    pos = (1:numDL)+sort(round(rand(1,numDL)*(K^2-numDL)));
    [pos1,pos2] = ind2sub([K K],pos);
    for cntDL = 1:numDL
        Vr((1:npc)+(pos1(cntDL)-1)*npc,(1:npc)+(pos2(cntDL)-1)*npc) = 0;
    end
end
if (isReal)
    % Vrfun(x) is Vr acting on x in the real space
    Vfun = @(x) repmat(Vr(:),1,numel(x)/n^2).*x;
else
    % Vfun(x) is V acting on x in the Fourier domain
    Vfun = @(x) reshape(fft2(reshape(repmat(Vr,1,numel(x)/n^2),n,n,numel(x)/n^2).*ifft2(reshape(x,n,n,numel(x)/n^2))),n^2,numel(x)/n^2);
end

% construct the discrete Hamiltonian operator H = T + V
Hfun = @(x) Tfun(x) + Vfun(x);






