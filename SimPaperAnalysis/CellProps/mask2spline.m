function [cspl, L] = mask2spline(imt,dstep,n)
%%  Create a spline along the centerline of a cell mask (imt). 
%       dstep = point spacing (default = 0.396).
%       n = smoothing size (default = 21).
    if (nargin<3)
        n = 21;
    end
    
    if (nargin<2)
        dstep = 0.396;
    end

    indxs = find(imt);
    [y, x] = ind2sub(size(imt),indxs);
    dx = x - x';
    dy = y - y';
    dr = sqrt(dx.^2 + dy.^2);
    [~,mind] = max(dr,[],'all');
    [cid,~] = ind2sub(size(dr),mind);
    [~,sortedinds] = sort(dr(cid,:));
    y = y(sortedinds);
    x = x(sortedinds);

    yy = smoothdata(y,'gaussian',n);
    xx = smoothdata(x,'gaussian',n);

    dxx = xx(2:end) - xx(1:end-1);
    dyy = yy(2:end) - yy(1:end-1);
    drr = sqrt(dxx.^2 + dyy.^2);

    cdist = 0;
    L = 0;
    slx = [xx(1)];
    sly = [yy(1)];

    for k = 1:numel(xx)-1

        dslx = slx(end) - xx(k);
        dsly = sly(end) - yy(k);

        dslx2 = slx(end) - xx(k+1);
        dsly2 = sly(end) - yy(k+1);

        dslr = sqrt(dslx^2 + dsly^2);
        dslr2 = sqrt(dslx2^2 + dsly2^2);

        if (dslr2 > dstep) && (dslr < dstep)
            frac = (dstep - dslr)/(dslr2 - dslr);
            ndx = frac*(xx(k+1) - xx(k));
            ndy = frac*(yy(k+1) - yy(k));
            L = L + sqrt(((xx(k) + ndx) - slx(end))^2 ...
                + ((yy(k) + ndy) - sly(end))^2);

            slx = [slx; xx(k) + ndx];
            sly = [sly; yy(k) + ndy];
        end
    end

    dxx = diff(xx);
    dyy = diff(yy);
    drr = sqrt(dxx.^2 + dyy.^2);
    L = sum(drr);
    cspl = [slx, sly];
end
