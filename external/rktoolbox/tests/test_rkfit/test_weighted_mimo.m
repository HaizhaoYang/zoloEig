function check = test_weighted_mimo()
  N = 40;
  A = gallery('lehmer', N);
  v = sqrt(1:N).';
  tol = 1e-13;

  F{1} = eye(N)/(A+2*speye(N));
  F{2} = F{1}/(A+3*speye(N));
  F{3} = F{2}/(A+4*speye(N));

  W{1} = diag(v.^0.2);
  W{2} = diag(v.^0.25);
  W{3} = diag(v.^1.3);

  xi = inf(3, 1);
  param.D = W;
  param.real = 1;
  param.reduction = 0;
  [xi, ratfun, misfit] = rkfit(F, A, v, xi, param);

  n1 = norm(sort(xi, 'descend')-[-2, -3, -4])*1e-5;
  n2 = misfit(end);

  nrm_denom = norm([norm(W{1}*F{1}*v) norm(W{3}*F{3}*v) norm(W{2}*F{2}*v)]);
  n3 = norm([norm(W{1}*(ratfun{1}(A, v)-F{1}*v)) ...
	     norm(W{2}*(ratfun{2}(A, v)-F{2}*v)) ...
	     norm(W{3}*(ratfun{3}(A, v)-F{3}*v))])/nrm_denom;

  check = [n1 n2 n3] < tol;
end
