function[allv,stencil]=allvstencil3D(v,h,nlocal,numelallv,reach,spread)
vreal=real(v(:));
vmin=min(vreal);
vmax=max(vreal);
%vmax=min(max(vreal),0);

allv=linspace(vmin,vmax,numelallv);
allv=pi^2*(4*unique(ceil(1/(4*pi^2)*allv(:)))-2);

numelallv=numel(allv);
%allv=allv+1i;

stencil(1:numelallv)=struct('q',[],'c',[]);
for countallv=1:numelallv
    stencil(countallv)=onevstencilconv3D(allv(countallv),h,nlocal,reach,spread);
end
