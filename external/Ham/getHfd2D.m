function [H,numEigVal] = getHfd2D(npc,K,opt)
% 4th order finite difference

% Input:
% npc - number of grid points per cell
% K - number of Gaussian wells per dimension
% opt - contains opt.isSC, opt.numDL, opt. sigma, specifying the Gaussian
%       wells in each cell
%
% Output:
% That -   The Laplacian operator defined in the Fourier domain
% Vr   -   The potential function defined in the real space
% H    -   A sparse matrix representing the Hamiltonian by a 9-point finite
%          difference scheme in the real domain
% numEigVal - the number of eigenvalues, the number of orbitals

if nargin < 3
    opt.depth = 1;
    opt.numDL = 0;
    opt.sigma = 0.2;
end

sigma = opt.sigma;
numDL = opt.numDL;
depth = opt.depth*K^2;

numEigVal = K^2;

c = 0.5;
h = 1/npc; % step size in space
n = npc*K; % size of grid in 1D space

% periodic boundary condition
TOneDim = sparse( [1:n 1:n-1 2:n], ...
    [1:n 2:n 1:n-1 ], ...
    [-2*ones(1,n),ones(1,n-1),ones(1,n-1)],n,n);
T = (kron(TOneDim,eye(n)) + kron(eye(n),TOneDim));
% we don't use K^2 in depth, correspondingly, we use npc^2 here,
% i.e., scale the H by 1/K^2;

% construct the potential energy operator
x = h*(0:npc-1)';
v = gaussmf(x,[sigma c]);
vv = -v*v'*depth;
Vr=repmat(vv,K,K);
if (numDL>0)
    pos = (1:numDL)+sort(round(rand(1,numDL)*(K^2-numDL)));
    [pos1,pos2] = ind2sub([K K],pos);
    for cntDL = 1:numDL
        Vr((1:npc)+(pos1(cntDL)-1)*npc,(1:npc)+(pos2(cntDL)-1)*npc) = 0;
    end
end

% Vx = linspace(0,1,size(Vr,1));
% Vy = linspace(0,1,size(Vr,2));
% imagesc(Vx,Vy,Vr);
% axis square;
% colorbar;

V = sparse(1:n^2,1:n^2,Vr(:));
%V = sparse(1:n^2,1:n^2,0.1*rand(1,n^2));

% construct the discrete Hamiltonian operator H = T + V
H = -T/2 + V;

end