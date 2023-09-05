function msk = spline2mask(cspl,boxSize,res,d)
    XYcal = boxSize/res;
    img = zeros(res);
    
    xs = cspl(:,1)/XYcal;
    ys = cspl(:,2)/XYcal;
    [xs,ys] = sortspline(xs,ys,res);

    xs2 = xs(1:end-1)+diff(xs)/2;
    ys2 = ys(1:end-1)+diff(ys)/2;
    xx = [xs; xs2];
    yy = [ys; ys2];
    %[xs, ys] = sortspline(xx,yy,res);

    inds = sub2ind(size(img),ceil(ys),ceil(xs));
    img(inds) = 1;
    msk = imdilate(img,strel('disk',d));
    msk = round(imgaussfilt(msk,d/2));

end