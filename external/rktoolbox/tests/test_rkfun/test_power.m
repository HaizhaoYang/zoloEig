function check = test_power
  
  x = [-.5, 1, 5, pi, 1i, 2+3i, 2-3i];
  
  I  = rkfun();
  C  = 2*I.^0;
  R1 = I.^3;
  R2 = I.^-2;
  R3 = R1.*R2 - I;
  
  c1 = norm(I(x)-x, inf);
  c2 = norm(C(x)-2*ones(size(x)), inf);
  c3 = norm(R1(x).*R2(x) - x, inf);
  c4 = norm(R3(x), inf);  
  
  R = rkfun.nodes2rkfun([10i, 20i],-12);
  C  = (2.*R).^0;
  R1 = R.^-5;
  R2 = R.^4;
  R3 = R.*(R2.*R1) - (2*C)./2;
  
  c5 = norm(C(x)-ones(size(x)), inf);
  c6 = norm(R1(x).*R2(x) - 1./R(x), inf);
  c7 = norm(R3(x), inf);
  
  check = [c1 c2 c3 c4 c5 c6 c7] < 5e-15;
  
end
