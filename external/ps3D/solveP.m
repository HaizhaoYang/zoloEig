function u = solveP(N,C, f)
    
    u = f;
    cntc = numel(C);
    for g=1:cntc
        rdidx = C(g).rdidx;
        skidx = C(g).skidx;
        u(skidx) = u(skidx) - C(g).X * u(rdidx);
    end
    for g=1:cntc
        rdidx = C(g).rdidx;
        skidx = C(g).skidx;
        u(rdidx) = C(g).Arrinv * u(rdidx);
    end
    for g=cntc:-1:1
        rdidx = C(g).rdidx;
        skidx = C(g).skidx;
        u(rdidx) = u(rdidx) - C(g).Y * u(skidx);
    end
end

