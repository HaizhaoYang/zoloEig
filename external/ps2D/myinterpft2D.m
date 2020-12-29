function[v]=myinterpft2D(v,n)
v=interpft(v,n);
v=permute(v,[2,1]);
v=interpft(v,n);
v=permute(v,[2,1]);
end