function [Uout,eigvals] = zologeig(A, B, varargin)

[nc,r,reltol,itmax,verbose,itsolacc] = zolooptions(varargin{:});

n = size(A,1);
nc = max(nc,ceil(sqrt(N)));

%========================================================================
% Zolo Eig
if verbose
    fprintf('-----------------------------------------\n');
    tottime = cputime;
    partime = cputime;
    fprintf('Slicing Calculation ...');
end

% find the largest spectrum range
[~,evsa] = eigs(A,B,1,'sa');
[~,evla] = eigs(A,B,1,'la');

a0 = evsa - 0.1;
b0 = evla + 0.1;

% roughly divide the spectrum range into pieces
opt = [];
opt.kcol = nc;
opt.cheborder = 50;
opt.ratio = 0.9;
slicing = MSSlicingChebg(A,B,a0,b0,a0,b0,a0,opt);

nSlices = numel(slicing) - 1;

maxgap = 0;
mingap = Inf;
aint = cell(1,nSlices);
bint = cell(1,nSlices);
for cnt = 1 : nSlices
    if cnt == 1
        aint{cnt} = [evsa-1, evsa];
    else
        aint{cnt} = bint{cnt-1};
    end
    if cnt == nSlices
        bint{cnt} = [evla,evla+1];
    else
        shift = slicing(cnt+1);
        mlt = 1;
        while mlt < n
            evshift = eigs(A-shift*B,B,4*mlt,'sm');
            if max(evshift)-min(evshift) > 1e-8
                break;
            end
            fprintf('Repeated %d eigenvalues = %e', mlt*4, shift);
            mlt = mlt*2;
        end
        evshift = sort(abs(evshift)) + shift;
        gap = evshift(2:end) - evshift(1:end-1);
        [gapval,posg] = max(gap);
        if gapval > maxgap
            maxgap = gapval;
        end
        if gapval < mingap
            mingap = gapval;
        end
        bint{cnt} = [evshift(posg),evshift(posg+1)];
    end
end

if verbose
    partime = cputime-partime;
    fprintf('\b\b\b      %12.2f secs\n',partime);
    fprintf('    Number of slicing            %d\n', nSlices );
        fprintf('    Max gap           %.2e\n',maxgap);
        fprintf('    Min gap           %.2e\n',mingap);
end

if nargout == 2
    Uout = zeros(n,n);
end
eigvals = zeros(n,1);

% apply zoloeigs in each spectrum range
offset = 0;
for cnt = 1:nSlices
    if nargout == 2
        [ Usub, evsub ] = zologeigs( A, B, aint{cnt}, bint{cnt}, nc, ...
            r, reltol, itmax, verbose, itsolacc );
        evsub = diag(evsub);
        Uout(:,offset+(1:size(Usub,2))) = Usub;
    else
        evsub = zologeigs( A, B, aint{cnt}, bint{cnt}, nc, ...
            r, reltol, itmax, verbose, itsolacc );
    end
    eigvals(offset+(1:numel(evsub))) = evsub;
    offset = offset + numel(evsub);
end

if verbose
    tottime = cputime-tottime;
    fprintf('Total time              %12.2f secs\n',tottime);
    fprintf('-----------------------------------------\n');
end

if nargout == 1
    Uout = eigvals;
    return;
end

end