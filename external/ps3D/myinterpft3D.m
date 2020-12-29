function[v]=myinterpft3D(v,n)
v=interpft(v,n);
v=permute(v,[2,3,1]);
v=interpft(v,n);
v=permute(v,[2,3,1]);
v=interpft(v,n);
v=permute(v,[2,3,1]);
end