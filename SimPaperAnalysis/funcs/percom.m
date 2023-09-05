function [xcom, ycom] = percom(bds,boxSize)
    xst = [bds.xs]*2*pi/boxSize;
    xscos = cos(xst);
    xssin = sin(xst);

    xcom = mod(atan2(sum(xssin),sum(xscos)),2*pi)*boxSize/(2*pi);

    yst = [bds.ys]*2*pi/boxSize;
    yscos = cos(yst);
    yssin = sin(yst);

    ycom = mod(atan2(sum(yssin),sum(yscos)),2*pi)*boxSize/(2*pi);

end
    