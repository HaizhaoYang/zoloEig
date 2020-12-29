function [func] = Ainvfunc(A,B,sigma)
func = @(x) (A-sigma*B)\x;
end