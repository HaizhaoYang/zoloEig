function [Hfun,That,Vr,Tfun,Vfun,numEigVal,Eval,n,Evecr,time] = getHpw2D(npc,K,opt)
%
% Input:
% isReal - 0: Hfun is the Hamiltonian discretized in the Fourier domain, i.e.
%             Hfun is a functional acting on complex functions defined in
%             the Fourier domain
%          1: Hfun is a hamiltonian discretized in the real space
%
% Output:
% That -   The Laplacian operator defined in the Fourier domain
% Vr   -   The potential function defined in the real space
% Hfun -   A function handle representing the Hamiltonian, defined in the Fourier
%          domain or real space, depending on isReal
% Tfun -   A function handle representing the Laplacian operator, defined in the Fourier
%          domain or real space, depending on isReal
% Vfun -   A function handle representing the potential term, defined in the Fourier
%          domain or real space, depending on isReal
% numEigVal - the number of eigenvalues, the number of orbitals
% Eval -   the eigenvalues of H
% n    -   the dimension of the problem in 1D
% Evecr -  the eigenvectors of the matrix representation of the Hamiltonian
%          defined in the real space

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
    opt.numDL = 3;
    opt.sigma = 0.2;
end
sigma = opt.sigma;
numAtom = 1;
numDL = opt.numDL;
depth = opt.depth*K^2;
numGaussian1D = K;

numEigVal = (numGaussian1D)^2*numAtom; % number of eigenvalues, orbitos
h = 1/npc;%step size in space
n = npc*K;% size of grid in 1D space
ks = [0:n/2-1 -n/2:-1]';%grid points in Fourier domain
[k1,k2] = ndgrid(ks);

% construct the Laplacian operator
k2diag = 2*pi^2*(k1.^2 + k2.^2); % diagonal elements of -0.5*Laplacian in Fourier domain, size n^2 in 2D
Tfun = @(x) real( reshape( ifft2(reshape(repmat(k2diag,1,numel(x)/n^2),n,n,numel(x)/n^2).*fft2(reshape(x,n,n,numel(x)/n^2))) ,n^2,numel(x)/n^2) );
% Trfun(x)is -0.5*Laplacian*f in the real space for real functions x

% construct the potential energy operator
c = 0.5;
x=h*(0:npc-1);
v = gaussmf(x,[sigma c]);
v = v(:);
vv = -v*v'*depth;
Vr=repmat(vv,numGaussian1D,numGaussian1D);
if (numDL>0)
    pos = (1:numDL)+sort(round(rand(1,numDL)*(K^2-numDL)));
    [pos1,pos2] = ind2sub([K K],pos);
    for cntDL = 1:numDL
        Vr((1:npc)+(pos1(cntDL)-1)*npc,(1:npc)+(pos2(cntDL)-1)*npc) = 0;
    end
end
% Vrfun(x) is Vr acting on x in the real space
Vorgfun = @(x) repmat(Vr(:),1,numel(x)/n^2).*x;


% construct the discrete Hamiltonian operator H = T + V
Horgfun = @(x) Tfun(x) + Vorgfun(x);

%% estimate the spectrum of the original problem by direct method
Hrorgfun = Horgfun;
tic;
[Evecr,Eval] = eigs(Horgfun,n^2,n^2-2,'sr');
time = toc;
[EvecMax,EvalMax] = eigs(Hrorgfun,n^2,2,'lr');%non-symetric case
Evecr = [Evecr EvecMax];
Eval = diag(Eval); EvalMax = diag(EvalMax); Eval = [Eval(:); EvalMax(:)];
Evecr = Evecr*diag(sum(Evecr.^2));

[Eval,ord] = sort(real(Eval));
Evecr = Evecr(:,ord);

Hfun = @(x) Horgfun(x);
Vfun = @(x) Vorgfun(x);
That = k2diag;






