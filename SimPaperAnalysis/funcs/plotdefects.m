function plotdefects(ax,adefs,XYcal)

    hold(ax,'on');
    L = 30;
    pdefs = adefs(abs([adefs.q]-0.5)< 0.1);
    ndefs = adefs(abs([adefs.q]+0.5)<0.1);
    plot(ax,[pdefs.x]/XYcal,[pdefs.y]/XYcal,'r.','MarkerSize',20);
    plot(ax,[ndefs.x]/XYcal,[ndefs.y]/XYcal,'b.','MarkerSize',20);

    for i = 1:numel(pdefs)
        plot(ax,[pdefs(i).x/XYcal pdefs(i).x/XYcal+L*pdefs(i).dx],...
            [pdefs(i).y/XYcal pdefs(i).y/XYcal+L*pdefs(i).dy],'r','LineWidth',2);
    end
    for i = 1:numel(ndefs)
        defang = atan2(ndefs(i).dy,ndefs(i).dx);
        for rot = 0:(2*pi/3):(4*pi/3)
            defangt = defang + rot;
            plot(ax,[ndefs(i).x/XYcal ndefs(i).x/XYcal + L*cos(defangt)],...
                [ndefs(i).y/XYcal ndefs(i).y/XYcal + L*sin(defangt)],...
                'b','LineWidth',2);
        end
    end
end