function check = test_arithm()
  ev = [-1, 0, pi, -1+2i, -1-2i];

  tol = 1e-12;

  % Check arithmetic operations.
  for run = 1:2
    if run == 1, x = ev;
    else         x = rkfun; end

    r = -x+2-3i;
    r = (r./2).*(-4);
    r = r.*(r.*x-7+2i);
    r = r./(1+x.^3-2*x);

    if run == 1, v1 = r;
    else         v2 = r(ev); end
  end

  check(1) = norm(v1-v2, inf);

  % Check composition.
  m = rkfun.gallery('moebius', [2 2-2i 8 1-3i]);
  s = r(m);
  check(2) = norm(s(ev)-r(m(ev)), inf);
  check = check < tol;
end
