function check = test_arithm_cplx()
  tol = 1e-13;
  
  x = [7+(1:10)/2 (11:20)+1i*(1:10)/3];
  
  a = rkfun.nodes2rkfun([1i, 2, -1i], [1-5i, 1+5i, 8i, -8i]);  
  b = rkfun.nodes2rkfun([6, 4],       [-1, -2, -3]);
  
  c1 = a + b; c2 = b + a;
  n1 = norm((a(x)+b(x)-c1(x))./(a(x)+b(x)), inf);
  n2 = norm((c2(x)-c1(x))./(a(x)+b(x)),     inf);

  c1  = a - b; c2  = b - a;
  n3 = norm((a(x)-b(x)-c1(x))./(a(x)-b(x)), inf);
  n4 = norm((c1(x)+c2(x))./(a(x)-b(x)),     inf);

  c1 = a .* b; c2 = b .* a;  
  n5 = norm((a(x).*b(x)-c1(x))./(a(x).*b(x)), inf);
  n6 = norm((c2(x)-c1(x))./(a(x).*b(x)),      inf);

  c1 = a ./ b; c2 = b ./ a; 
  n7 = norm((a(x)./b(x)-c1(x))./(a(x)./b(x)),  inf);
  n8 = norm((c2(x).*c1(x)-ones(1, length(x))), inf);

  c1 = a.^-3; c2 = (a.^3).^-1;
  n9 = norm((a(x).^-3-c1(x))./(a(x).^-3),  inf);
  n0 = norm((c1(x)-c2(x))./(c1(x).^3),     inf);

  check = [n1 n2 n3 n4 n5 n6 n7 n8 n9 n0] < tol;
end