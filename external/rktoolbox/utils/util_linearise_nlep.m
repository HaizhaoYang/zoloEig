function AB = util_linearise_nlep(A, Sigma, Xi, tol, Nmax, cyclic)
% UTIL_LINEARISE_NLEP    Constructs a pencil for the linearisation
%                        of a nonlinear eigenvalue problem.
%
% Input arguments:
% A         - function handle to eigenvalue problem A(lambda) of size n x n
% Sigma     - vector of vertices of the target set
% Xi        - vector vertices of pole set, or cyclically repeated poles 
%             (when 'cyclic' flag is set)
% tol       - tolerance for the construction of the linearisation
% Nmax      - the pencil will be of size at most n x Nmax
% 'cyclic'  - whether or not to repeat the poles in Xi cyclically
% 
% Output arguments:
% AB        - a pencil structure to be used by RAT_KRYLOV
%
% The linearisation constructed by this function is based on rational
% Newton interpolants and described in 
%
% S. Guettel, R. Van Beeumen, K. Meerbergen, and W. Michiels. 
% NLEIGS: A class of fully rational Krylov methods for nonlinear 
% eigenvalue problems, SIAM J. Sci. Comput., 36(6):A2842--A2864, 2014.
%
% CAUSE: It is assumed that none of the poles xi is close to zero. 
% For problems where poles near zero need to be used simply shift the 
% eigenvalue problem to A(lambda+s) and the target set to Sigma-s and
% the poles to Xi-s.
  
  if nargin < 6, cyclic = 0; else cyclic = 1; end
  if nargin < 5, Nmax = 200; end
  if nargin < 4, tol = 1e-8; end
  
  [D,sigma,xi,beta,nrmD,b,control] = sample_nlep(A, Sigma, Xi, tol, Nmax, cyclic);

  N = length(D)-1;
  AB.N = N;
  AB.beta = beta;
  AB.sigma = sigma;
  AB.control = control;
  AB.xi = xi;  
  AB.D = D;
  AB.nrmD = nrmD;
  AB.b = b;

  %%%%%%%%%%%%%%%%%%%%%%%%%
  % Construct the pencil. %
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  AB.isreal = isreal(sigma) && isreal(beta) && isreal(xi);
  if AB.isreal,
    for j = 1:N+1, if ~isreal(D{j}), AB.isreal = false; break, end, end    
  end
  
  AB.eval = @(z) lin_nlep_eval(AB,z);
  AB.multiply = @(rho, eta, x) lin_nlep_multiply(AB, rho, eta, x);
  lin_nlep_solve(); % reset cache
  AB.solve    = @(nu,  mu,  x) lin_nlep_solve(AB, nu, mu, x);
  AB.get_matrices = @() lin_nlep_get_matrices(AB);

end % function


function [AN,BN] = lin_nlep_get_matrices(AB)
% Return linearisation pencil (AN,BN) explicitly. 
% This function should only be called on small n x n problems as the
% linearisation is of size N*n x N*n;
    D = AB.D;
    sigma = AB.sigma;
    xi = AB.xi;
    beta = AB.beta;
    N = AB.N;
    n = size(D{1},1);
    AN = spdiags([sigma(1:N).',[0,beta(1:N-1)].'],-1:0,N,N);
    AN = kron(AN,speye(n));
    for j = 1:N,
      AN(1:n,1+n*(j-1):n*j) = D{j};  
    end
    BN = spdiags([ones(N,1),[0,beta(1:N-1)./xi(1:N-1)].'],-1:0,N,N);
    BN = kron(BN,speye(n));
    AN(1:n, end-n+1:end) = AN(1:n, end-n+1:end) - sigma(N)*D{end}/beta(N);
    BN(1:n, end-n+1:end) = BN(1:n, end-n+1:end) - D{end}/beta(N);

end

function QN = lin_nlep_eval(AB,z)
% Evaluate rational linearization QN(z) at a scalar argument z.

    QN = 0*AB.D{1};
    for j = 0:AB.N,
      QN = QN + AB.b(j,z)*AB.D{j+1};  
    end

end

function y = lin_nlep_multiply(AB, rho, eta, x)
% LIN_NLEP_MULTIPLY
  
  N = length(AB.D)-1;
  n = length(AB.D{1});
  
  % The pencil is of size N*n-by-N*n. The vectors x and y are
  % block-partitioned in N blocks of length n each.  
  assert(size(x, 1) == N*n)
  
  y = zeros(size(x));
  
  % We now compute the first block y(1:n, :) of y.
  if rho ~= 0
    jb = 1;
    je = n;
    for j = 1:N
      jb = (j-1)*n+1;
      je = j*n;
      y(1:n, :) = y(1:n, :) + rho*AB.D{j}*x(jb:je, :);
    end
  end
  
  jb = (N-1)*n+1;
  je = N*n;
  scalar = -rho*AB.sigma(N)/AB.beta(N)+eta/AB.beta(N);
  y(1:n, :) = y(1:n, :) + scalar*AB.D{N+1}*x(jb:je, :);  
  
  % And now the remaining blocks of y.
  for j = 2:N,
    jb = (j-1)*n+1;
    je = j*n;
    
    y(jb:je, :) = ...
        (rho*AB.sigma(j-1)-eta)*x((jb:je)-n, :) + ...
        (AB.beta(j-1)*(rho-eta/AB.xi(j-1)))*x(jb:je, :);
  end
    
end %lin_nlep_multiply


function y = lin_nlep_solve(AB, nu, mu, x)
% LIN_NLEP_SOLVE
  
  persistent CACHE  
  if nargin == 0, % activate/reset cache use   
    CACHE.L = {}; CACHE.U = {}; CACHE.P = {}; CACHE.Q = {}; 
    CACHE.shift = [];
    return
  end
    
  N = length(AB.D)-1;
  n = length(AB.D{1});
  
  % The pencil is of size N*n-by-N*n. The vectors x and y are
  % block-partitioned in N blocks of length n each.
  
  assert(size(x, 1) == N*n)
  
  y = zeros(size(x));    
  
  % Get the first block y(1:n, :) of y.
  % 1) Compute z, using a temporary w.
  % 2) Solve the system y(1:n, :) = A(mu/nu)\z.
  scalar = AB.beta(1)*(nu-mu/AB.xi(1));
  z = x(1:n, :)/nu - AB.D{2}*x((1:n)+n, :)/scalar;
  w = x((1:n)+n, :);
  for j = 2:N-1,
    jb = j*n+1;
    je = (j+1)*n;
    scalar = (nu*AB.sigma(j)-mu)/(AB.beta(j-1)*(nu-mu/AB.xi(j-1)));
    w = x(jb:je, :) - scalar*w;
    z = z - AB.D{j+1}*w / (AB.beta(j)*(nu-mu/AB.xi(j)));
  end
  scalar = (nu*AB.sigma(N)-mu)/(AB.beta(N)*nu);
  scalar = scalar/(AB.beta(N-1)*(nu-mu/AB.xi(N-1)));
  z = z + scalar*AB.D{N+1}*w;

  if isfield(CACHE,'shift'), % store LU factors of Q(shift)
      ind = find(CACHE.shift == mu/nu,1,'first');
      if isempty(ind),        % compute and store factor
         ind = length(CACHE.shift)+1;
         Q = lin_nlep_eval(AB,mu/nu);
         if issparse(Q),
            [CACHE.L{ind},CACHE.U{ind},CACHE.P{ind},CACHE.Q{ind}] = lu(Q);
         else
            [CACHE.L{ind},CACHE.U{ind},CACHE.P{ind}] = lu(Q);
            CACHE.Q{ind} = 1;
         end
         CACHE.shift(ind) = mu/nu;
      end
      y(1:n, :) = CACHE.Q{ind}*(CACHE.U{ind}\(CACHE.L{ind}\(CACHE.P{ind}*z)));
  else
      Q = lin_nlep_eval(AB,mu/nu);    
      y(1:n, :) = Q\z;
  end
  
  % Get the remaining N-1 blocks.
  for j = 2:N
    jb = (j-1)*n+1;
    je = j*n;
    
    a = nu*AB.sigma(j-1)-mu;
    b = AB.beta(j-1)*(nu-mu/AB.xi(j-1));
    
    y(jb:je, :) = (x(jb:je, :) - a*y((jb:je)-n, :))/b;    
  end
  
end % lin_nlep_solve


function [D,sigma,xi,beta,nrmD,b,Control] = sample_nlep(A, Sigma, Xi, tol, Nmax, cyclic, npts)

if nargin < 7, npts = 5000; end
if nargin < 6, cyclic = 0; end
if nargin < 5, Nmax = 200; end
if nargin < 4, tol = 1e-8; end

% For (A, Sigma, Xi, tol, Nmax, cyclic) being provided
% this code produces the corresponding data for the linearization

Control = util_discretise_polygon(Sigma,npts);
Sigma = Control;

if cyclic, % use the poles in Xi as provided, repeated in cyclic fashion
    if any(isnan(Xi)),
        error('All poles in Xi must be finite or infinite numbers.');
    end
else
    Xi = util_discretise_polygon(Xi,npts);
end

beta = []; % we assume b_0 = 1
sigma = Sigma(1); 
D{1} = A(sigma(1));
nrmD(1) = norm(D{1},'fro');
xi = Xi(1);
bControl = 1;
bXi = 1;    
j = 1;

while j <= Nmax,
    
    if abs(xi(j)) < 1e-6, 
        warning('SAMPLE_NLEP: One of the poles in the linearisation is very close to zero. Consider shifting the NLEP for stability.');
    end
    
    % get scaling parameter beta and next sigma
    bControl = bControl.*(Control - sigma(j))./(1 - Control/xi(j));
    [beta(j),ind] = max(abs(bControl));
    bControl = bControl/beta(j);
    bXi = bXi.*(Xi - sigma(j))./(1 - Xi/xi(j))/beta(j);
    sigma(j+1) = Control(ind);
    if cyclic, % cyclic poles
        xi(j+1) = Xi(mod(j,length(Xi))+1);
    else % leja-bagby
        [~,ind] = min(abs(bXi));
        xi(j+1) = Xi(ind);  
    end
    
    % Rational basis functions. See [p. A2847, GVMM14].
    b = @(j, Lambda) arrayfun(@(lambda) ...
    prod((lambda-sigma(1:j))./(beta(1:j).*(1-lambda./xi(1:j)))),Lambda);
    % Next divided difference.
    As = A(sigma(j+1));
    for k = 0:j-1,
        As = As - b(k,sigma(j+1))*D{k+1};
    end
    D{j+1} = As/b(j,sigma(j+1));
    j = j + 1;
    
    nrmD(j) = norm(D{j},'fro');
    if nrmD(j)/nrmD(1) < tol,
        break
    end
    
end
N = length(D) - 1;
sigma = sigma(1:N+1);
xi = xi(1:N); 
% change last pole to infinity
xi(N) = inf;
beta = beta(1:N);
b = @(j, Lambda) arrayfun(@(lambda) ...
    prod((lambda-sigma(1:j))./(beta(1:j).*(1-lambda./xi(1:j)))),Lambda);
D{N+1} = As/b(N,sigma(N+1));
nrmD(N+1) = norm(D{N+1},'fro');

end % sample_nlep
