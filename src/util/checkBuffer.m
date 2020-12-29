function [a1,b1] = checkBuffer(a1,b1,SIGMA1)
if SIGMA1 < a1
    if b1-a1<a1-SIGMA1
        b1 = a1;
        a1 = SIGMA1;
    end
else if SIGMA1>b1
        if b1-a1<SIGMA1-b1
            a1 = b1;
            b1 = SIGMA1;
        end
    end
end
end