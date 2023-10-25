%% 

disp(dataset)
vsperdef = 6000;

dx = 0.5;

xbins = -20:dx:20;
ybins = -20:dx:20;

pvx = zeros(numel(ybins),numel(xbins));
pvy = pvx;
pM2x = pvx;
pM2y = pvx;
pctx = pvx;
pcty = pvx;
nvx = pvx;
nvy = pvx;
nM2x = pvx;
nM2y = pvx;
nctx = pvx;
ncty = pvx;

switch dataset
    case {"mbwt", "mbpila", "mbflavo"}
        [fpath1,fpath2,fpath3] = exptpaths(dataset);
        N = numel(fpath3);
        minr = 7;
        rsm = 17;
        XYcal = 0.072;
        dsm = 1;
    case "kcwt"
        [fpath1,fpath2,~] = exptpaths(dataset);
        N = numel(fpath2);
        minr = 7;
        rsm = 9;
        XYcal = 0.133;
        dsm = 1;
    case "sim"
        eps = 5; d = 0.7; ll = 7; rho = 0.3; v = 5; Kagar = 500;
        Kstiff = 100; rev = 8; runs = 0:2; fp = "PoissonRevs";
        XYcal = 0.133;
        dsm = 1;
        minr = 7;
        N = numel(runs);
        rsm = 9;
        
end

for f = 1:N
    switch dataset
        case {"mbwt", "mbpila", "mbflavo"}
            files = dir(fullfile(fpath1,fpath2,fpath3(f),"img"));
            files = files(arrayfun(@(x) x.name(1)~='.',files));
            T = readtable(fullfile(fpath1,fpath2,fpath3(f),'frameTimeZ.csv'));
            ts = T.time-T.time(1);
        case "kcwt"
            files = dir(fullfile(fpath1,fpath2(f)));
            ts = getts(fullfile(fpath1,fpath2(f)));
        case "sim"
            fpath = simname(eps,d,ll,rho,v,Kagar,Kstiff,rev,runs(f),fp);
            [files,boxSize] = getframes(fpath);
            ts = (0:numel(files) - 1)*30;
    end

    OFlow = opticalFlowFarneback('NumPyramidLevels',3,'NumIterations',1,...
            'NeighborhoodSize',5,'FilterSize',21);

    switch dataset
        case {"mbwt", "mbpila", "mbflavo"}
            lraw = imread(fullfile(files(1).folder,files(1).name));
        case "kcwt"
            lraw = laserdata(files(1).folder,1);
        case "sim"
            bds = loadsimdata(fullfile(files(1).folder,files(1).name));
            lraw = cellim(bds,boxSize,XYcal,[0 0],3,2,2);
    end
    l = processl(lraw,dataset);
    sz = size(l);
    xs = (1:sz(2))*XYcal;
    ys = (1:sz(1))*XYcal;
    [xx,yy] = meshgrid(xs,ys);

    for t= 1:numel(ts)

        lp = l;
        switch dataset
            case {"mbwt", "mbpila", "mbflavo"}
                lraw = imread(fullfile(files(t).folder,files(t).name));
            case "kcwt"
                lraw = laserdata(files(1).folder,t);
            case "sim"
                bds = loadsimdata(fullfile(files(t).folder,files(t).name));
                lraw = cellim(bds,boxSize,XYcal,[0 0],2,1,0.5);
        end
        l = processl(lraw,dataset);

        % Find drift
        h = xcorr_fft(l,lp);
        p = xcorrpeak(h);
        sz = size(l);
        drx = sz(2)/2 - p(1);
        dry = sz(1)/2 - p(2);
        
        % Drift correct.
        if drx>0
            lcr = l(:,(ceil(drx)+1):end);
            lcr = padarray(lcr,[0 ceil(drx)],'replicate','post');
        elseif drx<0
            lcr = l(:,1:end-ceil(-drx));
            lcr = padarray(lcr,[0 ceil(-drx)],'replicate','pre');
        end
        
        if dry>0
            lcr = lcr((ceil(dry)+1):end,:);
            lcr = padarray(lcr,[ceil(dry) 0],'replicate','post');
        elseif dry<0
            lcr = lcr(1:end-ceil(-dry),:);
            lcr = padarray(lcr,[ceil(-dry) 0],'replicate','pre');
        end

        % Calculate flow.
        flow = estimateFlow(OFlow,lcr);

        if t == 1
            dt = (ts(2) - ts(1))/60;
        else
            dt = (ts(t) - ts(t-1))/60;
        end

        Vx = flow.Vx * XYcal / dt;
        Vy = flow.Vy * XYcal / dt;
        Vang = atan2(Vy,Vx);
        V = sqrt(Vx.^2 + Vy.^2);

        % Get all defects.    
        dirf = dfield_expt(lcr,rsm);
        holes = findholes(lcr,9);
        adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);

        % Loop through all defects.
        for i = 1:numel(adefs)
            dx = xx - adefs(i).x;
            dy = yy - adefs(i).y;
            dr = sqrt(dx.^2 + dy.^2);
            thet = atan2(dy,dx);
            if (adefs(i).q > 0.1)
                rotthet = thet - atan2(adefs(i).dy,adefs(i).dx);
                pts = (abs(dr.*cos(rotthet)) < 20).*(abs(dr.*sin(rotthet))<20);
                pts(holes) = 0;
                inds = find(pts);
                if (numel(inds) > vsperdef)
                    testpts = randperm(numel(inds),vsperdef);
                else
                    testpts = randperm(numel(inds));
                end
                rotvs = Vang - atan2(adefs(i).dy,adefs(i).dx);

                for ii = 1:numel(testpts)
                    xind = discretize(dr(inds(testpts(ii)))*...
                        cos(rotthet(inds(testpts(ii)))),xbins);
                    yind = discretize(dr(inds(testpts(ii)))*...
                        sin(rotthet(inds(testpts(ii)))),ybins);
                    
                    [pvx(yind,xind), pM2x(yind,xind),...
                        pctx(yind,xind)] = ...
                        running_mean(pvx(yind,xind),...
                        pM2x(yind,xind),pctx(yind,xind),...
                        V(inds(testpts(ii))).*cos(rotvs(inds(testpts(ii)))));                
                    [pvy(yind,xind), pM2y(yind,xind),...
                        pcty(yind,xind)] = ...
                        running_mean(pvy(yind,xind),...
                        pM2y(yind,xind),pcty(yind,xind),...
                        V(inds(testpts(ii))).*sin(rotvs(inds(testpts(ii)))));
                    
                end

            elseif (adefs(i).q < -0.1)
                for rot = 0:2*pi/3:4*pi/3
                    rotthet = thet - atan2(adefs(i).dy,adefs(i).dx)+rot;
                    pts = (abs(dr.*cos(rotthet)) < 20).*(abs(dr.*sin(rotthet))<20);
                    pts(holes) = 0;
                    inds = find(pts);
                    if (numel(inds) > vsperdef)
                        testpts = randperm(numel(inds),vsperdef);
                    else
                        testpts = randperm(numel(inds));
                    end                    
                    rotvs = Vang - atan2(adefs(i).dy,adefs(i).dx)+rot;
    
                    for ii = 1:numel(testpts)
                        xind = discretize(dr(inds(testpts(ii)))*...
                            cos(rotthet(inds(testpts(ii)))),xbins);
                        yind = discretize(dr(inds(testpts(ii)))*...
                            sin(rotthet(inds(testpts(ii)))),ybins);
                        [nvx(yind,xind), nM2x(yind,xind),...
                            nctx(yind,xind)] = ...
                            running_mean(nvx(yind,xind),...
                            nM2x(yind,xind),nctx(yind,xind),...
                            V(inds(testpts(ii))).*cos(rotvs(inds(testpts(ii)))));                
                        [nvy(yind,xind), nM2y(yind,xind),...
                            ncty(yind,xind)] = ...
                            running_mean(nvy(yind,xind),...
                            nM2y(yind,xind),ncty(yind,xind),...
                            V(inds(testpts(ii))).*sin(rotvs(inds(testpts(ii)))));  
                    end
                end
            end
        end

        save(fullfile('~/lammps-myxosim/SimPaperAnalysis/Figures/Defects/Flows/',[char(dataset), 't.mat']),...
            'pvx','pvy','pM2x','pM2y','pctx','pcty',...
            'nvx','nvy','nM2x','nM2y','nctx','ncty');

    end
end


save(fullfile('~/lammps-myxosim/SimPaperAnalysis/Figures/Defects/Flows/',dataset),...
    'pvx','pvy','pM2x','pM2y','pctx','pcty',...
    'nvx','nvy','nM2x','nM2y','nctx','ncty');
