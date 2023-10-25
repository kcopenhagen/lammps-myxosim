function defvid(fpath,savepath,dataset)
    addpath(genpath("~/SimPaperAnalysis"));
        
    switch dataset
        case "mbpila"
            L = 40;
        case "mbwt"
            L = 40;
        case "kcwt"
            L = 20;
        case "sim"
            L = 20;
    end

    v = VideoWriter(savepath,'Motion JPEG AVI');
    v.Quality = 95;
    v.FrameRate = 12;
    open(v);
    switch dataset
        case "sim"
            fig = figure('Units','pixels','Position',[0 0 500 500]);
        otherwise
            fig = figure('Units','pixels','Position',[0 0 1024/2 768/2]);
    end
    ax = axes(fig,'Position',[0 0 1 1]);

    switch dataset
        case "kcwt"
            ts = getts(fpath);
            minr = 7;
            rsm = 9;
            XYcal = 0.133;
            dsm = 1;
        case "mbwt"
            files = dir(fullfile(fpath,"img"));
            files = files(arrayfun(@(x) x.name(1)~='.',files));
            T = readtable(fullfile(fpath,'frameTimeZ.csv'));
            ts = T.time-T.time(1);
            minr = 7;
            rsm = 17;
            XYcal = 0.072;
            dsm = 1;
        case "mbpila"
            files = dir(fullfile(fpath,"img"));
            files = files(arrayfun(@(x) x.name(1)~='.',files));
            T = readtable(fullfile(fpath,'frameTimeZ.csv'));
            ts = T.time-T.time(1);
            minr = 7;
            rsm = 17;
            XYcal = 0.072;
            dsm = 1;
        case "sim"
            [files,boxSize] = getframes(fpath);
            ts = 30*(1:numel(files));
            XYcal = 0.072;
            rsm = 17;
            dsm = 1;
            minr = 7;
    end

    for t = 1:numel(ts)
        switch dataset
            case "mbwt"
                lraw = imread(fullfile(files(t).folder,files(t).name));
            case "mbpila"
                lraw = imread(fullfile(files(t).folder,files(t).name));
            case "kcwt"
                lraw = laserdata(fpath,t);
            case "sim"
                bds = loadsimdata(fullfile(files(t).folder,files(t).name));
                lraw = cellim(bds,boxSize,XYcal);
        end
        l = processl(lraw,dataset);

        % Get all defects.
        dirf = dfield_expt(l,rsm);
        holes = findholes(l,9);
        if dataset == "kcwt"
            holes = loaddata(fpath,t,'covid_layers','int8')==0;
        end
        holes = imdilate(holes,strel('disk',5));
        adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
        
        imagesc(l);
        hold on
        colormap gray

        pdefs = adefs([adefs.q] > 0.1);
        ndefs = adefs([adefs.q] < -0.1);

        plot([pdefs.x]/XYcal,[pdefs.y]/XYcal,'r.','MarkerSize',20);
        plot([ndefs.x]/XYcal,[ndefs.y]/XYcal,'b.','MarkerSize',20);

        for i = 1:numel(pdefs)
            plot([pdefs(i).x/XYcal pdefs(i).x/XYcal+L*pdefs(i).dx],...
                [pdefs(i).y/XYcal pdefs(i).y/XYcal+L*pdefs(i).dy],'r','LineWidth',2);
        end
        for i = 1:numel(ndefs)
            defang = atan2(ndefs(i).dy,ndefs(i).dx);
            for rot = 0:(2*pi/3):(4*pi/3)
                defangt = defang + rot;
                plot([ndefs(i).x/XYcal ndefs(i).x/XYcal + L*cos(defangt)],...
                    [ndefs(i).y/XYcal ndefs(i).y/XYcal + L*sin(defangt)],...
                    'b','LineWidth',2);
            end
        end
        F = getframe(fig);
        writeVideo(v,F);
        hold off
    end
    close(v)
end
    