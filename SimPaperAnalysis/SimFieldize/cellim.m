function img = cellim(bds,boxSize,XYcal,amt,cd,cg,cs)
%% bds is the data read in with loadsimdata2.
% boxSize is boxSize
% res is how many pixels on a side for the final image.
% amt is the shift (use min([bds.xs]) min([bds.ys]) on recentered data.).

    res = round(boxSize/XYcal);
    if (nargin < 4)
        amt = [0 0];
    end

    if (nargin < 5)
        cd = round(0.22*res/boxSize);
    end
    if (nargin < 6)
        cg = round(0.1*res/boxSize);
    end
    if (nargin < 7)
        cs = 2;
    end

    cellimg = zeros(res);
    indx = round(([bds.xs] - amt(1))*res/boxSize)+1;
    indy = round(([bds.ys] - amt(2))*res/boxSize)+1;
    
    gd = ones(size(indx));
    gd(indx>res) = 0;
    gd(indy>res) = 0;

    indx = indx(gd==1);
    indy = indy(gd==1);

    inds = sub2ind([res res], indy,indx);
    cellimg(inds) = 1;
    cellimg = imdilate(cellimg, strel('disk',cd));
    cellimg = round(imgaussfilt(cellimg, cs));
    cellimg = imgaussfilt(cellimg, cg);
    img = cellimg;
end