function [resid,xi,absterm,cnd,pf,Q] = residue(obj,method)
%RESIDUE   Convert a RKFUN into partial fraction form.
%
% [resid,xi,absterm,cnd,pf] = residue(ratfun)
% Converts a rational function with well-separated simple poles into 
% partial fraction form 
%
%   r(z) = absterm + sum_j resid(j)/(z - xi(j))
%
% The rational function should not be of superdiagonal type.
%
% The number cnd is the condition number of the transformation from the 
% rational Krylov basis to the partial fraction basis. 
% This transformation can be extremely badly conditioned and multiple 
% precision arithmetic may be required.
%
% The function handle pf can be used for scalar evaluation of the partial 
% fraction.

if obj.k > 0 && method < 2,
  error('RESIDUE: Conversion is currently not possible for superdiagonal rkfuns.');
end

if nargin < 2, method = 0; end
Q = NaN;

switch method, 
 case 0
  % standard method using transformations on pencil
  [resid,xi,absterm,cnd,Q] = residue_pencil(obj);
 case 1
  % standard method using transformations on pencil
  [resid,xi,absterm,cnd] = residue_pencil_qz(obj);
 case 2
  % use complex contour integration
  cnd = NaN;
  [resid,xi,absterm] = residue_contour(obj);
end

% function handle for evaluation
if isa(absterm,'mp') || isa(absterm,'sym'),
  pf = @(x) absterm;
  for j = 1:length(xi),
    pf = @(x) pf(x) + resid(j)./(x - xi(j));
  end
else 
  pf = @(x) absterm + arrayfun(@(xx) sum(resid./(xx - xi)),x);
end

end % function




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [resid,xi,absterm,cnd,Q] = residue_pencil(obj)
%Partial fraction form via transformation of the pencil (H,K). 
%Uses only elementary transforms so that vpa and mp work.

K = obj.K;
H = obj.H;
coeffs = obj.coeffs;
m = size(H,2);

% bring lower mxm part of (H,K) to diagonal form
Kt = K/K(2:end,:);
Ht = H/K(2:end,:);
[U,D] = eig(Ht(2:end,:));  % done because vpa does not support generalized eig
Q = blkdiag(1,U);
Kt = Q\Kt*U;
Ht = Q\Ht*U;

% zero first row of Kt
X = Q^0; % identity (works for vpa/mp) 
X(1,2:end) = -Kt(1,:);
Kt = X*Kt;
Ht = X*Ht;
Q = Q/X; % keep track of basis transform

% scale columns so that first row of H is all ones
X = diag(1./Ht(1,:));
Kt = Kt*X;
Ht = Ht*X;

% scale rows so that subdiagonal entries of K are all one
X = diag([1;1./diag(Kt(2:end,:))]);
Kt = X*Kt;
Ht = X*Ht;
Q = Q/X;

% get residues and poles
resid = Q\coeffs;
absterm = resid(1);
resid = resid(2:end);
cnd = NaN;
try cnd = cond(double(Q)); catch, end;
%xi = diag(Ht,-1); % bug in MATLAB R2012b VPA, does not return all subdiagonals
xi = diag(Ht(2:end,:));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [resid,xi,absterm,cnd] = residue_pencil_qz(obj)
%Partial fraction form via transformation of the pencil (H,K). 
%Uses ORDQZ.

[R,xi] = get_residues_qz(obj.K, obj.H);
R = [ eye(size(R,1),1) , R ];
resid = R\(obj.coeffs);
absterm = resid(1);
resid = resid(2:end);
cnd = cond(R);
end



function [R,xi] = get_residues_qz(K, H)
% If K and H are upper-Hessenberg (m+1)-by-m matrices related
% to the RAD
%                       A*V*K = V*H,
% then this function returns an (m+1)-by-m matrix R such that
%   V*R = [(A-xi(1)*I)\V(:, 1) | ... | (A-xi(m)*I)\V(:, 1)],
% where xi(j) = H(j+1, j)/K(j+1, j).
%
% It is assumed that all the poles xi(j) are pairwise distinct.
% 
% EXAMPLE
% A = gallery('tridiag', 50);
% b = ones(50, 1);
% xi = [-1, -2, -3];
% [V, K, H] = rat_krylov(A, b, xi);
% R = get_residues_qz(K, H);
% normb = norm(b);
% norm((A-xi(1)*eye(50))\b - V*R(:, 1)*normb)/norm(R(:, 1)*normb)
% norm((A-xi(2)*eye(50))\b - V*R(:, 2)*normb)/norm(R(:, 2)*normb)
% norm((A-xi(3)*eye(50))\b - V*R(:, 3)*normb)/norm(R(:, 3)*normb)
  
  m = size(H, 2);
  xi = diag(H, -1)./diag(K, -1);
  R = zeros(m+1, m);
  
  for j = 1:m
    % At step j, we work on the leading j-by-(j-1) part only.
    % Bring xi(j) to the leading position.     
    index = zeros(1, j); 
    index(end) = 1;
    [~, ~, QL, ZL] = ordqz(H(2:j+1, 1:j), K(2:j+1, 1:j), ...
                           eye(j), eye(j), index);
    % The following 3 lines can be implemented more efficiently.
    QL = blkdiag(1, QL);
    HL = QL*H(1:j+1, 1:j)*ZL;
    KL = QL*K(1:j+1, 1:j)*ZL;
            
    % Locally replace xi(j) with inf.
    [s, c] = clcangl(KL(1, 1), KL(2, 1));
    QL(1:2, :) = [c -s;s' c]*QL(1:2, :);   
    
    % Get the scaling factor gamma.
    HL(1:2, 1) = [c -s;s' c]*HL(1:2, 1);
    KL(1:2, 1) = [c -s;s' c]*KL(1:2, 1);
        
    gamma_inv = KL(1, 1)*QL(1, 1)/(HL(1, 1) - KL(1, 1)*xi(j));
    gamma_inv = gamma_inv + KL(1, 1)*QL(2, 1)/(HL(2, 1));
    gamma_inv = gamma_inv/2;
    
    R(1:j+1, j) = QL(1, 1:j+1)'*gamma_inv;
  end
end

function [s, c] = clcangl(x, y)
%CLC_ANGL computes sin and cos for a plane rotation.
%
% [s, c] = compute_angle(x, y) for scalars x and y computes the sin
% and cos for the plane rotation G given by  
%               |c  -s|
%               |s'  c|,
% such that (G*[x; y]) = [h; 0].
  t = -y/x;
  c = 1/sqrt(1+abs(t)^2);
  s = t'*c;
  % h = hypot(x, y);
  % c = x/h;
  % s = -conj(y)/h;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [resid,xi,absterm] = residue_contour(obj)
%Partial fraction form via contour integration.

absterm = feval(obj,1e290);
xi = poles(obj); xi = xi(:);
%if ~all(isfinite(xi)),
%   error('RESIDUE: The rational functions appears to be of superdiagonal type and ratfun2pfe does not currently handle such cases.');
%end
m = length(xi);
npts  = 100;
pts = exp(2i*pi*(1:npts)/npts);
for j = 1:m, 
    c = xi(j);
    dist = min(abs(xi([1:j-1,j+1:m])-c)); 
    if isempty(dist), dist = 1; end  % just one pole
    r = dist/2;
    if r < 1e-6,
       error('RESIDUE: Poles are not well separated and residue does not currently handle such cases.');
    end
    resid(j,1) = r/npts*sum(sort(feval(obj,c+r*pts).*pts));
end

end