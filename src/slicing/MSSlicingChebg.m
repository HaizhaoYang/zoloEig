function slicing = MSSlicingChebg( A, B, a, b, a0, b0, slicing, opt, Ys, W )

kcol = opt.kcol;
ratio = opt.ratio;
cheborder = opt.cheborder;

if nargin < 8
    Ys = [];
end

if nargin < 9
    W = randn( size(A,1), kcol );
end

% Hutchinson's method.
[AW, Ys] = chebJacSliceg( A, B, W, cheborder, a, b, a0, b0, 1, Ys );

tr = sumel(conj(W).*AW)/kcol;

if tr < kcol*ratio % with high probability, this will work
    slicing(end+1) = b;
else
    midab = (a+b)/2;
    slicing = MSSlicingChebg( A, B, a, midab, a0, b0, slicing, opt, Ys, W );
    slicing = MSSlicingChebg( A, B, midab, b, a0, b0, slicing, opt, Ys, W );
end

end