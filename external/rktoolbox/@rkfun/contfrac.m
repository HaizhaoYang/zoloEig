function [h,hh,absterm,cnd,cf,Q] = contfrac(obj)
%CONTFRAC   Convert RKFUN into continued fraction form.
%
% [h,hh,absterm,cnd,cf] = contfrac(obj)
% Convert a diagonal or subdiagonal rkfun of type (m+k,m), with k <= 0,
% into continued fraction form.
% This is done by bringing the upper-Hessenberg pencil (H,K)
% to tridiagonal-and-diagonal form. The vectors h and hh and the number 
% absterm are the coefficients of the continued fraction as follows:
%
%              1
% r_m(z) = -------------------------------------------------------- + absterm.
%          hh(1)*z + 1
%                    ----------------------------------------------
%                    h(1) + 1
%                           ---------------------------------------
%                           hh(2)*z +...+ 1
%                                         -------------------------
%                                         h(m-1) + 1
%                                                  ----------------
%                                                  hh(m)*z + 1/h(m)
%
% If the rational function r_m is superdiagonal of type (m,m-1), then the 
% continued fraction coefficients are such that 
%
% r_m(z) = hh(1)*z + 1
%                    ----------------------------------------------
%                    h(1) + 1
%                           ---------------------------------------
%                           hh(2)*z +...+ 1
%                                         -------------------------
%                                         h(m-1) + 1
%                                                  ----------------
%                                                  hh(m)*z + 1/h(m)
%
% and absterm is set to NaN (as the absolute term is already taken care 
% of by the coefficients in h and hh). 
%
% The number cnd is the condition number of the transformation from the 
% rational Krylov basis to the continued fraction basis. 
% This transformation can be badly conditioned and multiple precision 
% and the use of multiple precision arithmetic is recommended.
%
% The function handle cf allows for the scalar evaluation of the continued
% fraction.
%
% The algorithm is described in 
%
%     V. Druskin, S. Guettel, and L. Knizhnerman. _Compressing 
%     variable-coefficient exterior Helmholtz problems via RKFIT,_
%     MIMS EPrint 2016.53 (<http://eprints.ma.man.ac.uk/2511/>), 
%     Manchester Institute for Mathematical Sciences, 
%     The University of Manchester, UK, 2016.

if obj.k > 1,
    error('CONTFRAC: Conversion is currently only possible for rkfuns with k <= 1');
end

K = obj.K;
H = obj.H;
coeffs = obj.coeffs;
m = size(H,2);
if obj.k <= 0, % transform decomposition so that it starts with [b,r(A)b]
    Q = [ eye(m+1,1) , coeffs ];
else % otherwise start with [r(A)b,b]
    Q = [ coeffs , eye(m+1,1) ];
end
[UU,~,~] = svd(Q);
Q = [ Q , UU(:,3:end) ];
Kt = Q\K;
Ht = Q\H;

% bring lower mxm part of Kt to identity form
Ktt = Kt/Kt(2:end,:);
Htt = Ht/Kt(2:end,:);

if obj.k <= 0, % In this case, Ktt(1,1) is not exactly 0. 
    absterm = -Ktt(1,1); % Compute the negative absterm.
    % now make "exactly" subdiagonal and repeat last two steps
    coeffs(1) = coeffs(1) - absterm;
    Q = [ eye(m+1,1) , coeffs ];
    [UU,~,~] = svd(Q);
    Q = [ Q , UU(:,3:end) ];
    Kt = Q\K;
    Ht = Q\H;
    Ktt = Kt/Kt(2:end,:);
    Htt = Ht/Kt(2:end,:);
else
    absterm = NaN;  % for k == 1 there is no separate absterm (it's in the cf)
end % end of correction of absterm

% zero first row of Ktt 
X = Q^0; % identity (works for vpa/mp) 
X(1,3:end) = -Ktt(1,2:end);
Ktt = X*Ktt;
Htt = X*Htt;
Q = Q/X; % keep track of basis transforms

% zero elements in first row of Htt
X = Q(1:m,1:m)^0; % mxm identity (works for vpa/mp) 
X(1,2:end) = -Htt(1,2:end)/Htt(1,1);
Ktt = Ktt*X;
Htt = Htt*X;

% zero second row of Ktt
X = Q^0; %identity 
X(2,3:end) = -Ktt(2,2:end);
Ktt = X*Ktt;
Htt = X*Htt;
Q = Q/X;

% tridiagonalize lower mxm part of Htt
e1 = Q(1:m,1:m)^0; e1 = e1(:,1);  % first unit vector (works for vpa/mp)
[VV,WW,TT] = util_bilanczos(Htt(2:end,:),e1,e1);    % may break down for subdiag.
Httt = blkdiag(1,WW')*Htt*VV;
Kttt = blkdiag(1,WW')*Ktt*VV;
Q = Q/blkdiag(1,WW');

% diagonal scaling from left and right
r(1,1) = 1/Httt(1,1);
h(1,1) = -1/(Httt(2,1)*r(1));
ell = r; % initialize for vpa/mp
ell(1,1) = 1; ell(2,1) = 1; 
ell(3,1) = 1/(Httt(3,1)*h(1)*r(1));
for j = 2:m,
    r(j,1) = 1/(h(j-1)*ell(j)*Httt(j,j));
    h(j,1) = 1/(-1/h(j-1) - ell(j+1)*Httt(j+1,j)*r(j));
    if j < m,
        ell(j+2,1) = 1/(h(j)*Httt(j+2,j)*r(j));
    end
end 
hh = ell(2:end).*r;
Q = Q/diag(ell);
cnd = cond(double(Q));

% function handle for evaluation
if obj.k <= 0,
    cf = @(x) 0*x;
    for j = m:-1:1,
       cf = @(x) 1 ./ (hh(j)*x + 1./(h(j) + cf(x))); 
    end
    cf = @(x) absterm + cf(x);
else % obj.k == 1
    cf = @(x) Inf*x;
    for j = m:-1:1,
        cf = @(x) hh(j)*x + 1./(h(j) + 1./cf(x)); 
    end
end
    

end

