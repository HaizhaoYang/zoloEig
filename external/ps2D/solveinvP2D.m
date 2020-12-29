function[x]=solveinvP2D(invP,x)
numel_invP=numel(invP);
for invP_count=1:numel_invP-1
    idxin=invP(invP_count).idxin;
    idxbd=invP(invP_count).idxbd;
    x(idxin)=invP(invP_count).A*x(idxin);
    x(idxbd)=x(idxbd)-invP(invP_count).B*x(idxin);
end
for invP_count=numel_invP
    idxin=invP(invP_count).idxin;
    x(idxin)=invP(invP_count).A*x(idxin);
end
for invP_count=numel_invP-1:-1:1
    idxin=invP(invP_count).idxin;
    idxbd=invP(invP_count).idxbd;
    x(idxin)=x(idxin)-invP(invP_count).C*x(idxbd);
end
end