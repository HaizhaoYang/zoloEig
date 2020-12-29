function[invP]=setupinvP3D(P,n,boxsize)
%n/boxsize needs to be a power of 2
invP(1:2*(n/boxsize)^3)=struct('A',[],'B',[],'C',[],'D',[],...
    'idxin',[],'idxbd',[]);
count=0;
L=log2(n/boxsize)+1;
idx=zeros([n+1,n+1,n+1]);
idx(1:n,1:n,1:n)=reshape(1:n^3,[n,n,n]);
idx(n+1,:,:)=idx(1,:,:);
idx(:,n+1,:)=idx(:,1,:);
idx(:,:,n+1)=idx(:,:,1);

%building correction matrix
%bdr1,bdr2,bdr3 will be taken care by others.
[bdr1,bdr2,bdr3]=ndgrid([false(boxsize,1);true]);
idxbdr1=find(bdr1);
idxbdr2=find(bdr2);
idxbdr3=find(bdr3);

for l=1
    cursize=2^(l-1)*boxsize;
    curnum=n/cursize;
    buildgrid=true(cursize+1,1);
    [grid1,grid2,grid3]=ndgrid(buildgrid);
    curgrid=grid1|grid2|grid3;
    curin=false([cursize+1,cursize+1,cursize+1]);
    curin(2:cursize,2:cursize,2:cursize)=true;
    curbd=~curin;
    curgridin=curgrid&curin;
    curgridbd=curgrid&curbd;
    lclidx=zeros([cursize+1,cursize+1,cursize+1]);
    curmatsize=sum(curgrid(:));
    lclidx(curgrid)=1:curmatsize;
    lclidxin=lclidx(curgridin);
    lclidxbd=lclidx(curgridbd);
    boxidx=zeros([curnum,curnum,curnum]);
    for bigi3=1:curnum
        for bigi2=1:curnum
            for bigi1=1:curnum
                i3=(bigi3-1)*cursize+1:bigi3*cursize+1;
                i2=(bigi2-1)*cursize+1:bigi2*cursize+1;
                i1=(bigi1-1)*cursize+1:bigi1*cursize+1;
                glbidx=idx(i1,i2,i3);
                glbidxin=glbidx(curgridin);
                glbidxbd=glbidx(curgridbd);
                curmat=full(P(glbidx,glbidx));
                %correcting matrix
                curmat(idxbdr1,idxbdr1)=0;
                curmat(idxbdr2,idxbdr2)=0;
                curmat(idxbdr3,idxbdr3)=0;
                A=inv(curmat(lclidxin,lclidxin));
                B=curmat(lclidxbd,lclidxin);
                C=A*curmat(lclidxin,lclidxbd);
                D=curmat(lclidxbd,lclidxbd)-B*C;
                count=count+1;
                invP(count).A=A;
                invP(count).B=B;
                invP(count).C=C;
                invP(count).D=D;
                invP(count).idxin=glbidxin;
                invP(count).idxbd=glbidxbd;
                boxidx(bigi1,bigi2,bigi3)=count;
                clear A B C D;
            end
        end
    end
end

for l=2:L-1
    cursize=2^(l-1)*boxsize;
    curnum=n/cursize;
    buildgrid=[true;false(cursize/2-1,1);true;false(cursize/2-1,1);true];
    [grid1,grid2,grid3]=ndgrid(buildgrid);
    curgrid=grid1|grid2|grid3;
    curin=false([cursize+1,cursize+1,cursize+1]);
    curin(2:cursize,2:cursize,2:cursize)=true;
    curbd=~curin;
    curgridin=curgrid&curin;
    curgridbd=curgrid&curbd;
    lclidx=zeros([cursize+1,cursize+1,cursize+1]);
    curmatsize=sum(curgrid(:));
    lclidx(curgrid)=1:curmatsize;
    lclidxin=lclidx(curgridin);
    lclidxbd=lclidx(curgridbd);
    childrenidx=cell([2,2,2]);
    for c3=1:2
        for c2=1:2
            for c1=1:2
                childgrid=false([cursize+1,cursize+1,cursize+1]);
                childgrid(1+(c1-1)*(cursize/2):1+c1*(cursize/2),...
                    1+(c2-1)*(cursize/2):1+c2*(cursize/2),...
                    1+(c3-1)*(cursize/2):1+c3*(cursize/2))=true;
                childrenidx{c1,c2,c3}=lclidx(curgrid&childgrid);
            end
        end
    end
    prevboxidx=boxidx;
    boxidx=zeros([curnum,curnum,curnum]);
    for bigi3=1:curnum
        for bigi2=1:curnum
            for bigi1=1:curnum
                i3=(bigi3-1)*cursize+1:bigi3*cursize+1;
                i2=(bigi2-1)*cursize+1:bigi2*cursize+1;
                i1=(bigi1-1)*cursize+1:bigi1*cursize+1;
                glbidx=idx(i1,i2,i3);
                glbidxin=glbidx(curgridin);
                glbidxbd=glbidx(curgridbd);
                curmat=zeros(curmatsize,curmatsize);
                for c3=1:2
                    for c2=1:2
                        for c1=1:2
                            curmat(childrenidx{c1,c2,c3},childrenidx{c1,c2,c3})=...
                                curmat(childrenidx{c1,c2,c3},childrenidx{c1,c2,c3})+...
                                invP(prevboxidx(2*bigi1+c1-2,2*bigi2+c2-2,2*bigi3+c3-2)).D;
                            invP(prevboxidx(2*bigi1+c1-2,2*bigi2+c2-2,2*bigi3+c3-2)).D=[];
                        end
                    end
                end
                A=inv(curmat(lclidxin,lclidxin));
                B=curmat(lclidxbd,lclidxin);
                C=A*curmat(lclidxin,lclidxbd);
                D=curmat(lclidxbd,lclidxbd)-B*C;
                count=count+1;
                invP(count).A=A;
                invP(count).B=B;
                invP(count).C=C;
                invP(count).D=D;
                invP(count).idxin=glbidxin;
                invP(count).idxbd=glbidxbd;
                boxidx(bigi1,bigi2,bigi3)=count;
                clear A B C D;
            end
        end
    end
end

for l=L
    cursize=2^(l-1)*boxsize;
    curnum=n/cursize;
    buildgrid=[true;false(cursize/2-1,1);true;false(cursize/2-1,1);true];
    [grid1,grid2,grid3]=ndgrid(buildgrid);
    curgrid=grid1|grid2|grid3;
    curin=false([cursize+1,cursize+1,cursize+1]);
    curin(1:cursize,1:cursize,1:cursize)=true;
    curgridin=curgrid&curin;
    lclidx=zeros([cursize+1,cursize+1,cursize+1]);
    curmatsize=sum(curgridin(:));
    lclidx(curgridin)=1:curmatsize;
    lclidx(cursize+1,:,:)=lclidx(1,:,:);
    lclidx(:,cursize+1,:)=lclidx(:,1,:);
    lclidx(:,:,cursize+1)=lclidx(:,:,1);
    childrenidx=cell([2,2,2]);
    for c3=1:2
        for c2=1:2
            for c1=1:2
                childgrid=false([cursize+1,cursize+1,cursize+1]);
                childgrid(1+(c1-1)*(cursize/2):1+c1*(cursize/2),...
                    1+(c2-1)*(cursize/2):1+c2*(cursize/2),...
                    1+(c3-1)*(cursize/2):1+c3*(cursize/2))=true;
                childrenidx{c1,c2,c3}=lclidx(curgrid&childgrid);
            end
        end
    end
    prevboxidx=boxidx;
    for bigi3=1:curnum
        for bigi2=1:curnum
            for bigi1=1:curnum
                i3=(bigi3-1)*cursize+1:bigi3*cursize+1;
                i2=(bigi2-1)*cursize+1:bigi2*cursize+1;
                i1=(bigi1-1)*cursize+1:bigi1*cursize+1;
                glbidx=idx(i1,i2,i3);
                glbidxin=glbidx(curgridin);
                curmat=zeros(curmatsize,curmatsize);
                for c3=1:2
                    for c2=1:2
                        for c1=1:2
                            curmat(childrenidx{c1,c2,c3},childrenidx{c1,c2,c3})=...
                                curmat(childrenidx{c1,c2,c3},childrenidx{c1,c2,c3})+...
                                invP(prevboxidx(2*bigi1+c1-2,2*bigi2+c2-2,2*bigi3+c3-2)).D;
                            invP(prevboxidx(2*bigi1+c1-2,2*bigi2+c2-2,2*bigi3+c3-2)).D=[];
                        end
                    end
                end
                A=inv(curmat);
                count=count+1;
                invP(count).A=A;
                invP(count).idxin=glbidxin;
                clear A;
            end
        end
    end
end
invP=invP(1:count);
end