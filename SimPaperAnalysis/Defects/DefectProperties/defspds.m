%datasets = ["kcwt"; "mbwt"; "mbpila"; "sim"];
datasets = ["sim"];
for dataseti = 1:numel(datasets)
    dataset = datasets(dataseti);
    disp(dataset)
    minstep = 5;
    switch dataset
        case "kcwt"
            [fpath1,fpath2,~] = exptpaths(dataset);
            N = numel(fpath2);
            minr = 2.4;
            rsm = 13;
            XYcal = 0.133;
            dsm = 1;
        case "mbwt"
            [fpath1,fpath2,fpath3] = exptpaths(dataset);
            N = numel(fpath3);
            minr = 2.4;
            rsm = 17;
            XYcal = 0.072;
            dsm = 1;
        case "mbpila"
            [fpath1,fpath2,fpath3] = exptpaths(dataset);
            N = numel(fpath3);
            minr = 2.4;
            rsm = 17;
            XYcal = 0.072;
            dsm = 1;
        case "sim"
            eps = 5; d = 0.7; l = 7; rho = 0.3; v = 5; Kagar = 500;
            Kstiff = 100; rev = 8; runs = 0:2; fp = "PoissonRevs";
            N = numel(runs);
            minr = 2.4;   % Minumum defect spacing. (um)
            rsm = 1;    % Director field smoothing for simulated director fields (um).
            XYcal = 0.2; % XY calibration for simulation fields (um / pix).
            dsm = 1;    % Gradient smoothing for calculating defect directions (um).
    end
    
    pvs = [];
    nvs = [];
    eds = 0:0.5:10;
    edx = 0.25:0.5:10;
    for f = 1:N
        switch dataset
            case "kcwt"
                fpath = fullfile(fpath1,fpath2(f));
                ts = getts(fpath);
            case "mbwt"
                fpath = fullfile(fpath1,fpath2,fpath3(f),"img");
                files = dir(fpath);
                files = files(arrayfun(@(x)x.name(1)~='.',files));
                T = readtable(fullfile(fpath1,fpath2,fpath3(f),'frameTimeZ.csv'));
                ts = T.time-T.time(1);
            case "mbpila"
                fpath = fullfile(fpath1,fpath2,fpath3(f),"img");
                files = dir(fpath);
                files = files(arrayfun(@(x)x.name(1)~='.',files));
                T = readtable(fullfile(fpath1,fpath2,fpath3(f),'frameTimeZ.csv'));
                ts = T.time-T.time(1);
            case "sim"
                fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,runs(f),fp);
                [files,boxSize] = getframes(fpath);
                ts = 1:numel(files);
                ts = 30*(ts-1);
        end
    
        deftrks = [];
        
        drift = [0 0];
    
        for t = 1:numel(ts)
            if (t > 1)
                lp = l;
            end
            switch dataset
                case "sim"
                    fname = fullfile(files(t).folder,files(t).name);
                    bds = loadsimdata(fname);
                    adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm);
                case "mbwt"
                    fname = fullfile(files(t).folder,files(t).name);
                    l = imread(fname);
                    holes = findholes(l,9);
                    l = rgb2gray(l);
                    l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                    dirf = dfield_expt(l,rsm);
                    adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
                case "kcwt"
                    l = laserdata(fpath,t);
                    dirf = dfield_expt(l,rsm);
                    holes = loaddata(fpath,t,'covid_layers','int8')==0;
                    holes = imdilate(holes,strel('disk',5));
                    adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
                case "mbpila"
                    fname = fullfile(files(t).folder,files(t).name);
                    l = imread(fname);
                    holes = findholes(l,9);
                    l = rgb2gray(l);
                    l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                    dirf = dfield_expt(l,rsm);
                    adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
            end
            if (t > 1)
                h = xcorr_fft(l,lp);
                p = xcorrpeak(h);
                sz = size(l);
                drx = sz(2)/2 - p(1);
                dry = sz(1)/2 - p(2);
                drift = [drift; drx dry];
            end
            ids = num2cell(-1*ones(numel(adefs)));
            defts = num2cell(ts(t)*ones(numel(adefs)));
            deftts = num2cell(t*ones(numel(adefs)));
            [adefs.id] = ids{:};
            [adefs.t] = defts{:};
            [adefs.tt] = deftts{:};
            
            deftrks = [deftrks; adefs];
        end
    
        cdefs = deftrks([deftrks.tt] == 1);
        ids = num2cell(1:numel(cdefs));
        [deftrks(1:numel(cdefs)).id] = ids{:};
    
        for t = 2:numel(ts)
            cdefs = deftrks([deftrks.tt] == t);
            ldefs = deftrks([deftrks.tt] == t-1);
            if (numel(cdefs) > 0) && (numel(ldefs) > 0)
                dx = [cdefs.x] - [ldefs.x]';
                dy = [cdefs.y] - [ldefs.y]';
                dr = sqrt(dx.^2 + dy.^2);
                dq = [cdefs.q] - [ldefs.q]';
                
        
                [~,I] = min(dr,[],1);
                if numel(ldefs) == 1
                    I = ones(1,numel(cdefs));
                end
                J = 1:numel(cdefs);
        
                [~,II] = min(dr',[],1);
                if numel(cdefs) == 1
                    II = ones(1,numel(ldefs));
                end
                JJ = 1:numel(ldefs);
        
                trks = (II(I) == J);
        
                pts = sub2ind(size(dr),I(trks),J(trks));
                gd = (abs(dq(pts)) < 0.1).*(dr(pts)<minstep)==1;
                
                trks = find(trks);
                trks = trks(gd);
        
                ntrks = num2cell([ldefs(I(trks)).id]);
                [cdefs(J(II(I(trks)))).id] = ntrks{:};
            end
            ntrks = num2cell(max([deftrks.id])+1:max([deftrks.id])+sum([cdefs.id]==-1));
            [cdefs([cdefs.id]==-1).id] = ntrks{:};
            deftrks([deftrks.tt] == t) = cdefs;
        end
    
        defids = unique([deftrks.id]);
    
        for i = 1:numel(defids)
            cdef = deftrks([deftrks.id] == defids(i));
            if (numel(cdef) > 3)
                driftx = drift([cdef(2:end).tt],1);
                drifty = drift([cdef(2:end).tt],2);
    
                dt = diff([cdef.t])/60;
                vx = diff([cdef.x])-driftx'*XYcal;
                vy = diff([cdef.y])-drifty'*XYcal;
                v = sqrt(vx.^2 + vy.^2)./dt;
                if (abs(cdef(1).q - 0.5) < 0.1)
                    pvs = [pvs; v'];
                elseif (abs(cdef(1).q + 0.5) < 0.1)
                    nvs = [nvs; v'];
                end
            end
        end
        plot(edx,histcounts(pvs,eds,'Normalization','pdf'),'r')
        hold on
        plot(edx,histcounts(nvs,eds,'Normalization','pdf'),'b')
        drawnow
    end
    save(dataset,'pvs','nvs');
end
