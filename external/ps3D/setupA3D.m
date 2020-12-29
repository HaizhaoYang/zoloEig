function[A]=setupA3D(v,h)
n=size(v,1);
% N=n*n*n;
A1D=setupA1D(zeros(n,1),h);
A=kron(eye(n),kron(eye(n),A1D))+kron(eye(n),kron(A1D,eye(n)))+kron(A1D,kron(eye(n),eye(n)))+diag(v(:));
end