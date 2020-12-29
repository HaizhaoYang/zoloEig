function[stencil]=onevstencilconv2D(onev,h,n,reach,spread,ecut)
%n is local n.

if nargin < 6
    ecut = Inf;
end

kk=[(0:(n+mod(n,2))/2-1),(-(n-mod(n,2))/2:-1)];
kk=(2*pi*kk).^2;
[kk1,kk2]=ndgrid(kk);
ksq=kk1+kk2;
ksq(ksq>ecut)=ecut;
uhat=1./((1/(n*h)^2)*ksq+onev);

u=ifft2(uhat);
ureach=u(mymod(1-2*reach:1+2*reach,n),mymod(1-2*reach:1+2*reach,n));

uhatsq=abs(uhat).^2;
udot=ifft2(uhatsq);
udotreach=udot(mymod(1-2*reach:1+2*reach,n),mymod(1-2*reach:1+2*reach,n));

Bdot=zeros(1+2*reach,1+2*reach,1+2*reach,1+2*reach);
Breach=zeros(1+2*reach,1+2*reach,1+2*reach,1+2*reach);
for i2=0:2*reach
    for i1=0:2*reach
        Bdot(:,:,i1+1,i2+1)=udotreach(1+2*reach-i1:1+4*reach-i1,1+2*reach-i2:1+4*reach-i2);
        Breach(:,:,i1+1,i2+1)=ureach(1+2*reach-i1:1+4*reach-i1,1+2*reach-i2:1+4*reach-i2);
    end
end
Bdot=reshape(Bdot,(1+2*reach)^2,(1+2*reach)^2);
Breach=reshape(Breach,(1+2*reach)^2,(1+2*reach)^2);

q=cell(reach,1);
c=cell(reach,1);
for countreach=1:reach
    tempidxin=false(1+2*reach,1+2*reach);
    tempidxin(1+reach+(-countreach:countreach),1+reach+(-countreach:countreach))=true;
    Bin=Breach(tempidxin,tempidxin);
    tempspread=min(countreach,spread);
    tempidxspread=false(1+2*reach,1+2*reach);
    tempidxspread(1+reach+(-tempspread:tempspread),1+reach+(-tempspread:tempspread))=true;
    Bspread=Breach(tempidxspread,tempidxin);
    Bnotspreaddot=Bdot(tempidxin,tempidxin)-Bspread'*Bspread;
    Bnotspreaddot=0.5*(Bnotspreaddot+Bnotspreaddot');
    
%     idxinin=false(1+2*countreach,1+2*countreach);
%     idxinin(1+countreach,1+countreach)=true;
    
%     matc=Bin\Bnotspreaddot/Bin;
%     temp=matc(~idxinin,~idxinin)\matc(~idxinin,idxinin);
%     tempc=zeros((1+2*countreach)^2,1);
%     tempc(idxinin)=1;
%     tempc(~idxinin)=-temp;
%     tempq=Bin\tempc;
    
    [V,D]=eig(Bnotspreaddot);
    [~,I]=min(diag(D));
    tempq=V(:,I);
    tempc=Bin*tempq;

%     temp=Bnotspreaddot(~idxinin,~idxinin)\Bnotspreaddot(~idxinin,idxinin);
%     tempq=zeros((1+2*countreach)^2,1);
%     tempq(idxinin)=1;
%     tempq(~idxinin)=-temp;
%     tempc=Bin*tempq;

    q{countreach}=tempq;
    c{countreach}=tempc;
end
stencil=struct('q',{q},'c',{c});
end
