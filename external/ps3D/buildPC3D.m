function[P,C]=buildPC3D(v,h,boxsize,allv,stencil,reach,spread)
n=size(v,1);
N=n*n*n;
numelallv=numel(allv);

[boxreach1,boxreach2,boxreach3]=ndgrid(min([boxsize,1:boxsize-1],boxsize:-1:1));
boxreach=min(min(min(boxreach1,boxreach2),boxreach3),reach);
% boxreach=ones(boxsize,boxsize);

I=zeros(0,1);
J=zeros(0,1);
SP=zeros(0,1);
SC=zeros(0,1);
idx=reshape(1:N,[n,n,n]);

allvrealpad=[-inf;real(allv);inf];
vreal=real(v(:));
[~,pos]=histc(vreal,allvrealpad);
flagnear=abs(allvrealpad(pos)-vreal)<abs(allvrealpad(pos+1)-vreal);
pos(flagnear)=pos(flagnear)-1;
pos=reshape(pos,[n,n,n]);

vlarge=zeros(n+2*reach,n+2*reach,n+2*reach);
vlarge(reach+1:reach+n,reach+1:reach+n,reach+1:reach+n)=v;
vlarge([1:reach,reach+n+1:2*reach+n],:,:)=vlarge([n+1:n+reach,1+reach:2*reach],:,:);
vlarge(:,[1:reach,reach+n+1:2*reach+n],:)=vlarge(:,[n+1:n+reach,1+reach:2*reach],:);
vlarge(:,:,[1:reach,reach+n+1:2*reach+n])=vlarge(:,:,[n+1:n+reach,1+reach:2*reach]);

idxlarge=zeros(n+2*reach,n+2*reach,n+2*reach);
idxlarge(reach+1:reach+n,reach+1:reach+n,reach+1:reach+n)=idx;
idxlarge([1:reach,reach+n+1:2*reach+n],:,:)=idxlarge([n+1:n+reach,1+reach:2*reach],:,:);
idxlarge(:,[1:reach,reach+n+1:2*reach+n],:)=idxlarge(:,[n+1:n+reach,1+reach:2*reach],:);
idxlarge(:,:,[1:reach,reach+n+1:2*reach+n])=idxlarge(:,:,[n+1:n+reach,1+reach:2*reach]);

for bigi3=1:boxsize:n
    for bigi2=1:boxsize:n
        for bigi1=1:boxsize:n
            for smalli3=1:boxsize
                for smalli2=1:boxsize
                    for smalli1=1:boxsize
                        i3=bigi3+smalli3-1;
                        i2=bigi2+smalli2-1;
                        i1=bigi1+smalli1-1;
                        localreach=boxreach(smalli1,smalli2,smalli3);
                        localvin=vlarge(i1+reach-localreach:i1+reach+localreach,i2+reach-localreach:i2+reach+localreach,i3+reach-localreach:i3+reach+localreach);
                        localpos=pos(i1,i2,i3);
                        localq=stencil(localpos).q{localreach};
                        localc=stencil(localpos).c{localreach};
                        localidxj=idxlarge(i1+reach-localreach:i1+reach+localreach,i2+reach-localreach:i2+reach+localreach,i3+reach-localreach:i3+reach+localreach);
                        localidxi=ones(numel(localidxj),1)*idxlarge(i1+reach,i2+reach,i3+reach);
                        localp=localq+localc.*(localvin(:)-allv(localpos));
                        I=[I;localidxi(:)];
                        J=[J;localidxj(:)];
                        SP=[SP;localp(:)];
                        SC=[SC;localc(:)];
                    end
                end
            end
        end
    end
end
P=sparse(I,J,SP,N,N);
C=sparse(I,J,SC,N,N);
end
