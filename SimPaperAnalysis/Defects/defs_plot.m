function defs_plot(bds,adefs)
    figure
    L = 3;
    XYcal = 0.2;
    boxSize = 100;
    res = 752;
    sXYcal = boxSize/res;
    cim = cellim(bds,100,res,[0 0], 1,1,0.5);
    %plot([bds.xs],[bds.ys],'k.','MarkerSize',2);
    p = pcolor(0:sXYcal:boxSize-sXYcal,0:sXYcal:boxSize-sXYcal,cim);
    p.EdgeColor ='none';
    colormap gray
    axis off
    axis equal
    hold on
    pdefs = adefs([adefs.q]>0);
    ndefs = adefs([adefs.q]<0);
    plot([pdefs.x],[pdefs.y],'r.','MarkerSize',20);
    plot([ndefs.x],[ndefs.y],'b.','MarkerSize',20);
    for p = 1:numel(pdefs)
        plot([pdefs(p).x pdefs(p).x+L*pdefs(p).dx],...
            [pdefs(p).y pdefs(p).y+L*pdefs(p).dy],'r','LineWidth',2)
    end
    for n = 1:numel(ndefs)
        ang = atan2(ndefs(n).dy,ndefs(n).dx);
        ang2 = ang+2*pi/3;
        ang3 = ang+4*pi/3;
        xs = [ndefs(n).x; ndefs(n).x+L*cos(ang); ndefs(n).x; ...
            ndefs(n).x+L*cos(ang2); ndefs(n).x; ndefs(n).x+L*cos(ang3)];
        ys = [ndefs(n).y; ndefs(n).y+L*sin(ang); ndefs(n).y; ...
            ndefs(n).y+L*sin(ang2); ndefs(n).y; ndefs(n).y+L*sin(ang3)];
        plot(xs,ys,'b','LineWidth',2);
    end

end
        