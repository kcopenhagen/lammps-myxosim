
function [inds, rimsz] = makeringelement(r)
    x = -r-1:r+1;
    y = -r-1:r+1;
    [xx,yy] = meshgrid(x,y);
    dr = sqrt(xx.^2 + yy.^2);
    seneigh = (abs(dr-r) < 0.5);

    % seb = strel('disk',r);
    % sea = strel('disk',r-1);
    % sea = padarray(sea.Neighborhood,[1 1]);
    % seneigh = seb.Neighborhood - sea;
    [y, x] = find(seneigh);
    x = x-numel(seneigh(:,1))/2;
    y = y-numel(seneigh(:,1))/2;
    a = atan2(y,x);
    [~,ind] = sort(a);
    x = x(ind);
    y = y(ind);

    x = x+numel(seneigh(:,1))/2;
    y = y+numel(seneigh(:,1))/2;
    inds = sub2ind(size(seneigh),y,x);
    inds = [inds; inds(1)];
    rimsz = size(seneigh,1);
    
end
    