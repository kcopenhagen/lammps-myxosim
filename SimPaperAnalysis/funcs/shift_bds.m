function bds = shift_bds(bds,x,y,boxSize,sz)
    xs = [bds.xs];
    ys = [bds.ys];
    xs = xs - x;
    ys = ys - y;
    xs = xs + sz/2;
    ys = ys + sz/2;

    xs(xs<0) = xs(xs<0) + boxSize;
    xs(xs>boxSize) = xs(xs>boxSize) - boxSize;
    
    ys(ys<0) = ys(ys<0) + boxSize;
    ys(ys>boxSize) = ys(ys>boxSize) - boxSize;

    del = zeros(size(xs));
    del(xs<0) = 1;
    del(xs>sz) = 1;
    del(ys<0) = 1;
    del(ys>sz) = 1;
    
    xst = num2cell(xs);
    yst = num2cell(ys);
    [bds.xs] = xst{:};
    [bds.ys] = yst{:};
    bds(del==1) = [];

end