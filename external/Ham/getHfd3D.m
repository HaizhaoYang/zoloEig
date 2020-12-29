function [H,numEigVal] = getHfd3D(npc,K,opt)
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
depth = opt.depth*K^3;

numEigVal = K^3;

c = 0.5;
h = 1/npc; % step size in space
n = npc*K; % size of grid in 1D space

% periodic boundary condition
TOneDim = sparse( [1:n 1:n-1 2:n], ...
    [1:n 2:n 1:n-1 ], ...
    [-2*ones(1,n),ones(1,n-1),ones(1,n-1)],n,n);
TTwoDim = (kron(TOneDim,eye(n)) + kron(eye(n),TOneDim));
T = (kron(TTwoDim,eye(n)) + kron(eye(n),TTwoDim));
% we don't use K^2 in depth, correspondingly, we use npc^2 here,
% i.e., scale the H by 1/K^2;

% construct the potential energy operator
x = h*(0:npc-1)';
v = gaussmf(x,[sigma c]);
vvv = -kron(v,kron(v,v))*depth;
Vr=repmat(reshape(vvv,npc,npc,npc),K,K,K);
if (numDL>0)
    pos = (1:numDL)+sort(round(rand(1,numDL)*(K^3-numDL)));
    [pos1,pos2,pos3] = ind2sub([K K K],pos);
    for cntDL = 1:numDL
        Vr((1:npc)+(pos1(cntDL)-1)*npc,(1:npc)+(pos2(cntDL)-1)*npc,(1:npc)+(pos3(cntDL)-1)*npc) = 0;
    end
end
V = sparse(1:n^3,1:n^3,Vr(:));

% construct the discrete Hamiltonian operator H = T + V
H = -T/2 + V;

end