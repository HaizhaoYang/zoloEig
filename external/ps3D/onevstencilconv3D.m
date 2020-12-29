function[stencil]=onevstencilconv3D(onev,h,n,reach,spread)
%n is local n.

kk=[(0:(n+mod(n,2))/2-1),(-(n-mod(n,2))/2:-1)];
kk=(2*pi*kk).^2;
[kk1,kk2,kk3]=ndgrid(kk);
uhat=1./((1/(n*h)^2)*(kk1+kk2+kk3)+onev);

u=ifftn(uhat);
ureach=u(mymod(1-2*reach:1+2*reach,n),mymod(1-2*reach:1+2*reach,n),mymod(1-2*reach:1+2*reach,n));

uhatsq=abs(uhat).^2;
udot=ifftn(uhatsq);
udotreach=udot(mymod(1-2*reach:1+2*reach,n),mymod(1-2*reach:1+2*reach,n),mymod(1-2*reach:1+2*reach,n));

Bdot=zeros(1+2*reach,1+2*reach,1+2*reach,1+2*reach,1+2*reach,1+2*reach);
Breach=zeros(1+2*reach,1+2*reach,1+2*reach,1+2*reach,1+2*reach,1+2*reach);
for i3=0:2*reach
    for i2=0:2*reach
        for i1=0:2*reach
            Bdot(:,:,:,i1+1,i2+1,i3+1)=udotreach(1+2*reach-i1:1+4*reach-i1,1+2*reach-i2:1+4*reach-i2,1+2*reach-i3:1+4*reach-i3);
            Breach(:,:,:,i1+1,i2+1,i3+1)=ureach(1+2*reach-i1:1+4*reach-i1,1+2*reach-i2:1+4*reach-i2,1+2*reach-i3:1+4*reach-i3);
        end
    end
end
Bdot=reshape(Bdot,(1+2*reach)^3,(1+2*reach)^3);
Breach=reshape(Breach,(1+2*reach)^3,(1+2*reach)^3);

q=cell(reach,1);
c=cell(reach,1);
for countreach=1:reach
    tempidxin=false(1+2*reach,1+2*reach,1+2*reach);
    tempidxin(1+reach+(-countreach:countreach),1+reach+(-countreach:countreach),1+reach+(-countreach:countreach))=true;
    Bin=Breach(tempidxin,tempidxin);
    tempspread=min(countreach,spread);
    tempidxspread=false(1+2*reach,1+2*reach,1+2*reach);
    tempidxspread(1+reach+(-tempspread:tempspread),1+reach+(-tempspread:tempspread),1+reach+(-tempspread:tempspread))=true;
    Bspread=Breach(tempidxspread,tempidxin);
    Bnotspreaddot=Bdot(tempidxin,tempidxin)-Bspread'*Bspread;
    Bnotspreaddot=0.5*(Bnotspreaddot+Bnotspreaddot');
    
%     idxinin=false(1+2*countreach,1+2*countreach,1+2*countreach);
%     idxinin(1+countreach,1+countreach,1+countreach)=true;

%     matc=Bin\Bnotspreaddot/Bin;
%     temp=matc(~idxinin,~idxinin)\matc(~idxinin,idxinin);
%     tempc=zeros((1+2*countreach)^3,1);
%     tempc(idxinin)=1;
%     tempc(~idxinin)=-temp;
%     tempq=Bin\tempc;
    
    [V,D]=eig(Bnotspreaddot);
    [~,I]=min(diag(D));
    tempq=V(:,I);
    tempc=Bin*tempq;
    
%     temp=Bnotspreaddot(~idxinin,~idxinin)\Bnotspreaddot(~idxinin,idxinin);
%     tempq=zeros((1+2*countreach)^3,1);
%     tempq(idxinin)=1;
%     tempq(~idxinin)=-temp;
%     tempc=Bin*tempq;
    
    q{countreach}=tempq;
    c{countreach}=tempc;
end
stencil=struct('q',{q},'c',{c});
end
