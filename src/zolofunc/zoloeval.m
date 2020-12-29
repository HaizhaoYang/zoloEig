function yval = zoloeval(xval,l,r,scale)

if nargin < 4
    scale = 0;
end

func = zolofunc(l,r,scale);
yval = func(xval);

end