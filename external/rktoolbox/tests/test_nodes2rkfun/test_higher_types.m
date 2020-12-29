function check = test_higher_types()

tol = 1e-14;

ev = pi - exp(1)*1i;
rts = [];
r = rkfun.nodes2rkfun(rts);
c1 = norm(type(r))==0 && r(ev)==1;

pls = [];
r = rkfun.nodes2rkfun(rts,pls);
c2 = norm(type(r))==0 && r(ev)==1;

scl = -pi+2i; 
r = rkfun.nodes2rkfun(rts,pls,scl);
c3 = norm(type(r))==0 && r(ev)==scl;

rts = [ 1 , -2, 7i+3, pi-1i , 8, 2 ];
r = rkfun.nodes2rkfun(rts,[],scl);
rev = scl*prod(ev-rts);
c4 = norm(type(r)-[6,0])==0 && abs(r(ev)-rev)/abs(rev)<tol;

pls = [ -3i-1, pi, 0, 0 , -3i+1 ];
r = rkfun.nodes2rkfun(rts,pls,scl);
rev = scl*prod(ev-rts)/prod(ev-pls);
c5 = norm(type(r)-[6,5])==0 && abs(r(ev)-rev)/abs(rev)<tol;

rts = [ 0, -1+1i, 7, -3, inf, -1-1i , 0, 2i, -2i ];
pls = [ inf, inf, 0, -1i, 1i ];
r = rkfun.nodes2rkfun(rts,pls);
rev = 1*prod(ev-rts(isfinite(rts)))/prod(ev-pls(isfinite(pls)));
c6 = norm(type(r)-[8,3])==0 && abs(r(ev)-rev)/abs(rev)<tol && isreal(r);

check = [c1 c2 c3 c4 c5 c6];

end

