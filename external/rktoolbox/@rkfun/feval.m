function w = feval(obj, varargin)
%FEVAL    Evaluate RKFUN at scalar or matrix arguments.
%
% Calling syntax: w = feval(obj, A);
%                 w = feval(obj, A, v);
%
% The rational function r represented by the obj variable can
% either be evaluated pointwise (if only A is provided), or as a
% matrix function r(A)*v if both A and v are provided.
%
% This function is used mainly by RKFIT which is described in
%
% [1] M. Berljafa and S. G{\"u}ttel. Generalized rational Krylov
%     decompositions with an application to rational approximation,
%     SIAM J. Matrix Anal. Appl., 36(2):894--916, 2015.

switch nargin
    case 2
        AA = varargin{1};
        if isa(AA, 'rkfun')
            % Composition of two rkfuns.
            if length(AA.coeffs) > 2
                error(['FEVAL: Composition is currently possible only for' ...
                    ' inner rkfuns of type at most [1,1].']);
            end
            % Get Moebius transform parameters (rho*z - eta)/(nu*z - mu).
            if length(AA.coeffs) == 1
                rho = 0; eta = AA.coeffs(1); mu = 1; nu = 0;
            elseif length(AA.coeffs) == 2
                rho = AA.coeffs(1)*AA.K(2, 1) - AA.coeffs(2)*AA.K(1, 1);
                eta = AA.coeffs(1)*AA.H(2, 1) - AA.coeffs(2)*AA.H(1, 1);
                mu  = AA.H(2, 1);
                nu  = AA.K(2, 1);
            end
            w = obj;
            w.K = (rho*obj.K - nu*obj.H);
            w.H = (eta*obj.K - mu*obj.H);
            
            if isreal(w)
                [HH, KK, Q, Z] = qz(w.H(2:end, :), w.K(2:end, :), 'real');
                w.K(1, :) = w.K(1, :)*Z; w.K(2:end, :) = KK;
                w.H(1, :) = w.H(1, :)*Z; w.H(2:end, :) = HH;
                w.coeffs = blkdiag(1, Q)*w.coeffs;
            end
            
            % figure out new type
            [t1,t2] = type(AA);
            if t2 >= 1,
                w.k = 0; % reset k for concatenation with single-pole function
            end
            return
        else
            % Numerical point-wise evaluation.
            [mA, nA] = size(AA);
            AA = AA(:);
            A.isreal = isreal(AA);
            A.multiply = @(rho, eta, x) rho*(AA.*x) - eta*x;
            A.solve    = @(nu,  mu,  x) (nu*AA-mu).\x;
            v = ones(mA*nA, 1);
        end
    case 3
        % Matrix function times a vector.
        A = varargin{1};
        v = varargin{2};
        assert(isstruct(A) || size(A, 1) == size(A, 2), ...
            'FEVAL: A should be a square matrix.');
        assert(isstruct(A) || size(A, 2) == size(v, 1), ...
            'FEVAL: A and v are not conformable.');
        assert(~isempty(A) && ~isempty(v), ...
            'FEVAL: Neither A nor v can be empty.');
end

%{
  for j = 1:size(v, 2)
    if size(obj.K, 1) == 1
      % Constant function of type [0, 0].
      V = v(:, j);
    else
      V = rat_krylov(A, v(:,j), obj.K, obj.H);
    end
    %res = norm(diag(AA)*V*obj.K - V*obj.H)/(norm(AA)*norm(V)*norm(obj.K)+norm(V)*norm(obj.H));
    %disp(['FEVAL: residual in rerunning the decomposition = ' num2str(res) ])
    w(:, j) = V*obj.coeffs;
    % Allocate more space (works for vpa/mp).
    if j == 1 && size(v, 2) > 1
      w(1, size(v, 2)) = 0;
    end
  end
%}
%
N = size(v, 1);
n = size(v, 2);

if isstruct(A)
    AB = A; clear A;
else
    AB.multiply = @(rho, eta, x) rho*(A*x) - eta*x;
    AB.solve    = @(nu, mu, x) (nu*A-mu*speye(N))\x;
    AB.isreal = isreal(A);
end

Pencil.multiply = @(rho, eta, x) multiply_(N,n,AB,rho,eta,x);
Pencil.solve    = @(nu,  mu,  x) solve_(N,n,AB,nu,mu,x);
try
    Pencil.isreal = AB.isreal;
catch
    Pencil.isreal = false;
end

if size(obj.K, 1) == 1
    V = v(:);
else
    V = rat_krylov(Pencil, v(:), obj.K, obj.H);
end

for j = 1:n
    jb = (j-1)*N+1;
    je = j*N;
    
    w(:, j) = V(jb:je, :)*obj.coeffs;
    % Allocate more space (works for vpa/mp).
    if j == 1 && size(v, 2) > 1
        w(1, size(v, 2)) = 0;
    end
end
%}

if nargin == 2, w = reshape(w, mA, nA); end

    function X = multiply_(N, n, AB, a, b, X)
        % Used to represent blkdiag(AB, AB, ..., AB)
        % in order to exploit level 3 operations.
        if ~isempty(obj.prefact)
            prefact = obj.prefact;
            X = a*(prefact.MFB\(prefact.A*reshape(X, N, n))) ...
                - b*reshape(X, N, n);
            X = X(:);
            return;
        end
        X = AB.multiply(a, b, reshape(X, N, n));
        X = X(:);
    end

    function X = solve_(N, n, AB, a, b, X)
        % Used to represent blkdiag(AB, AB, ..., AB)
        % in order to exploit level 3 operations.
        
        if ~isempty(obj.prefact)
            prefact = obj.prefact;
            it = find(prefact.mu(prefact.nu == a) == b);
            if ~isempty(it)
                X = prefact.MF{it}\(prefact.B*reshape(X, N, n));
                X = X(:);
                return;
            end
        end
        X = AB.solve(a, b, reshape(X, N, n));
        X = X(:);
        
    end

end % function

