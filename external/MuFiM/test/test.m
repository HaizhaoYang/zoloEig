A = getHfd2D(4,4);

tic;
AMF = Multifrontal(A);
toc;

n = size(A,1);

tic;
AinvX = AMF\randn(n,1);
toc;

Aapprox = AMF*eye(size(A));
disp('MF*I')
norm(A-Aapprox)/norm(full(A))

Aapprox = eye(size(A))*AMF;
disp('I*MF')
norm(A-Aapprox)/norm(full(A))

Ainvapprox = eye(size(A))/AMF;
disp('I/MF')
norm(eye(size(A))-Ainvapprox*A)

Ainvapprox = AMF\eye(size(A));
disp('MF\I')
norm(eye(size(A))-Ainvapprox*A)


A = A + speye(size(A))*1i;

AMF = Multifrontal(A);

Aapprox = AMF*eye(size(A));
disp('MF*I')
norm(A-Aapprox)/norm(full(A))

Aapprox = eye(size(A))*AMF;
disp('I*MF')
norm(A-Aapprox)/norm(full(A))

Ainvapprox = eye(size(A))/AMF;
disp('I/MF')
norm(eye(size(A))-Ainvapprox*A)

Ainvapprox = AMF\eye(size(A));
disp('MF\I')
norm(eye(size(A))-Ainvapprox*A)

A = randn(5^4);
A = A+A';

AMF = Multifrontal(A);

Aapprox = AMF*eye(size(A));
disp('MF*I')
norm(A-Aapprox)/norm(full(A))

Aapprox = eye(size(A))*AMF;
disp('I*MF')
norm(A-Aapprox)/norm(full(A))

Ainvapprox = eye(size(A))/AMF;
disp('I/MF')
norm(eye(size(A))-Ainvapprox*A)

Ainvapprox = AMF\eye(size(A));
disp('MF\I')
norm(eye(size(A))-Ainvapprox*A)