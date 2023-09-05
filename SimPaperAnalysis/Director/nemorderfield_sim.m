function S = nemorderfield_sim(df,rr)
% Calculates local order field of image in fpath, at time t.
%    Saves into analysis/order folder as an array of floats.
%% Calc order    
    %rr = round(r/XYcal);
    se = strel('disk',rr);
    sea = strel('disk',rr-1);
    sea = padarray(sea.Neighborhood,[1 1]);
    mf = se.Neighborhood;%-sea;
    
    norm = conv2(ones(size(df)),mf,'same');

    dfield2 = 2*df;
    dx2 = cos(dfield2);
    dy2 = sin(dfield2);

    mdx2 = conv2(dx2,mf,'same')./norm;
    mdy2 = conv2(dy2,mf,'same')./norm;

    mdf = atan2(mdy2,mdx2)/2;

    %Prefered direction (n).
    nx = cos(mdf);
    ny = sin(mdf);

    %Director field mean over areas.
    dx = cos(df);
    dy = sin(df);

    mdxdx = conv2(dx.*dx,mf,'same')./norm;
    mdxdy = conv2(dx.*dy,mf,'same')./norm;
    mdydy = conv2(dy.*dy,mf,'same')./norm;

    S = (2*(nx.*nx.*mdxdx+2*nx.*ny.*mdxdy+ny.*ny.*mdydy)-1);

end