n = 1000;

fprintf('matrix size:   n = %4.0f \n', n);
fprintf(['Residuals and orthogonality measure should be of' ...
    ' the order of 10^(-15).\n']);

fprintf('=====Symmetric eigendecomposition=====\n');
A = randn(n); A=A'+A;
tic;
[V,D] = zoloeig(A);
toc;
fprintf('ZOLOEIG QR :       residual = %9.2e, orthogonality = %9.2e\n', ...
    norm(A-V*diag(D)*V','fro')/norm(A,'fro'), ...
    norm(V'*V-eye(n),'fro')/sqrt(n));

tic;
[V,D] = zoloeiginv(A);
toc;
fprintf('ZOLOEIG INV:       residual = %9.2e, orthogonality = %9.2e\n', ...
    norm(A-V*diag(D)*V','fro')/norm(A,'fro'), ...
    norm(V'*V-eye(n),'fro')/sqrt(n));

tic;
[V,D] = eig(A);
toc;
fprintf(['MATLAB eig:       residual = %9.2e,' ...
    ' orthogonality = %9.2e\n\n'], ...
    norm(A-V*D*V','fro')/norm(A,'fro'), ...
    norm(V'*V-eye(n),'fro')/sqrt(n));
