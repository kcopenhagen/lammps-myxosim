function [xs2, ys2] = sortspline(xs,ys,res)
    dxs = xs - xs';
    dys = ys - ys';
    dxs(dxs>res/2) = dxs(dxs>res/2) - res;
    dys(dys>res/2) = dys(dys>res/2) - res;
    dxs(dxs<-res/2) = dxs(dxs<-res/2) + res;
    dys(dys<-res/2) = dys(dys<-res/2) + res;

    drs = sqrt(dxs.^2 + dys.^2);
    [~,i] = max(drs,[],'all');
    [~,ii] = ind2sub(size(drs),i);

    [~,I] = sort(drs(ii,:));
    xs2 = xs(I);
    ys2 = ys(I);

end

