function zz = util_discretise_polygon(z,npts)
% ZZ = UTIL_DISCRETISE_POLYGON(Z,NPTS)
% Discretise a polygon with vertices given by the (complex) entries of Z.
% Several polygons can be separated by NaN entries in Z. The function 
% also supports infinite vertices in which case a logarithmic spacing 
% is used on the corresponding edge.

if nargin < 2, npts = 5e3; end
if nargin < 1, z = 0; end
if length(z) == 1,
    if isinf(z),
        zz = inf(1,npts);
    else
        zz = z + exp(2i*pi*(0:npts-1)/npts);
    end
    return
end

% estimate total length to be discretized, all finite edges get a 
% share of npts points; all infinite pieces get npts/10 points by 
% default.

dz = diff(z); dz = abs(dz(isfinite(dz)));
L = sum(dz);
zz = [];
for j = 1:length(z)-1,
    if isfinite(z(j)) && isfinite(z(j+1)),
        np = ceil(1+npts*abs(z(j+1)-z(j))/L);
        %zz = [ zz , linspace(z(j),z(j+1),np) ];
        zz = [ zz , z(j)+(z(j+1)-z(j))*(1-cos(pi*(0:np)/np))/2 ];
    end
    if isinf(z(j)) && isfinite(z(j+1)),
        d = sign(real(z(j))) + 1i*sign(imag(z(j)));
        np = ceil(1+npts/10);
        zz = [ zz , d*logspace(-12,12,np)+z(j+1) ];
    end
    if isfinite(z(j)) && isinf(z(j+1)),
        d = sign(real(z(j+1))) + 1i*sign(imag(z(j+1)));
        np = ceil(1+npts/10);
        zz = [ zz , d*logspace(-12,12,np)+z(j) ];
    end
end
zz = unique(zz);

