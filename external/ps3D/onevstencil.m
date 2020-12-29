function[stencil]=onevstencil(onev,h,n,reach,spread)
%n is local n.
f=zeros([n,n,n]);
f(1,1,1)=1;
u=applyG3D(f,onev,n,h);
u=reshape(u,[n,n,n]);
u=fftshift(u);
N=n^3;
mid=n/2+1;
B=zeros([n,n,n,1+2*reach,1+2*reach,1+2*reach]);
for i3=-reach:reach
for i2=-reach:reach
    for i1=-reach:reach
        B(:,:,:,i1+1+reach,i2+1+reach,i3+1+reach)=u(mymod(-i1+(1:n),n),mymod(-i2+(1:n),n),mymod(-i3+(1:n),n));
    end
end
end
B=reshape(B,N,(1+2*reach)^3);
idxreach=false([n,n,n]);
idxreach(mid+(-reach:reach),mid+(-reach:reach),mid+(-reach:reach))=true;
Bnotreach=B(~idxreach,:);
Breach=B(idxreach,:);
q=cell(reach,1);
c=cell(reach,1);
for countreach=1:reach
    tempidxin=false([1+2*reach,1+2*reach,1+2*reach]);
    tempidxin(1+reach+(-countreach:countreach),1+reach+(-countreach:countreach),1+reach+(-countreach:countreach))=true;
    Bout=[Breach(~tempidxin,tempidxin);Bnotreach(:,tempidxin)];
    Bin=Breach(tempidxin,tempidxin);
%     [~,~,vsvd]=svd(Bnotindot);
    idxinin=false([1+2*countreach,1+2*countreach,1+2*countreach]);
    idxinin(1+countreach,1+countreach,1+countreach)=true;
    temp=Bout(:,~idxinin)\Bout(:,idxinin);
    tempq=zeros((1+2*countreach)^3,1);
    tempq(idxinin)=1;
    tempq(~idxinin)=-temp;
    q{countreach}=tempq;
    c{countreach}=Bin*q{countreach};
end
stencil=struct('q',{q},'c',{c});
end