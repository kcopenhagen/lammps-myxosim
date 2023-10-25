datasets = ["kcwt","mbwt","mbpila"];
%datasets = ["sim"];
ndefdata = [];
for datai = 1:numel(datasets)
    dataset = datasets(datai);
    disp(dataset)
    
    switch dataset
        case "sim"
            boxSize = 100;
            minr = 2.4;   % Minumum defect spacing. (um)
            rsm = 1;    % Director field smoothing for simulated director fields (um).
            XYcal = 0.2; % XY calibration for simulation fields (um / pix).
            dsm = 1;    % Gradient smoothing for calculating defect directions (um).
        case "mbwt"
            minr = 2.4;
            rsm = 17;
            XYcal = 0.072;
            dsm = 1;
        case "kcwt"
            minr = 2.4;
            rsm = 13;
            XYcal = 0.133;
            dsm = 1;
        case "mbpila"
            minr = 2.4;
            rsm = 17;
            XYcal = 0.072;
            dsm = 1;
    end
    
    fpaths = [];
    switch dataset
        case "kcwt"
            [fpath1, fpath2, ~] = exptpaths(dataset);
            N = numel(fpath2);
            for f = 1:N
                fpaths = [fpaths; string(fullfile(fpath1,fpath2(f)))];
            end
        case "mbwt"
            [fpath1, fpath2, fpath3] = exptpaths(dataset);
            N = numel(fpath3);
            for f = 1:N
                fpaths = [fpaths; string(fullfile(fpath1,fpath2,fpath3(f)))];
            end
        case "mbpila"
            [fpath1, fpath2, fpath3] = exptpaths(dataset);
            N = numel(fpath3);
            for f = 1:N
                fpaths = [fpaths; string(fullfile(fpath1,fpath2,fpath3(f)))];
            end
        case "sim"
            eps = 5;
            d = 0.7;
            l = 7;
            rho = 0.3;
            v = 5;
            Kagar = 500;
            Kstiff = 100;
            rev = 8;
            runs = 0:2;
            fp = "PoissonRevs";
            N = numel(runs);
            for f = 1:N
                fpaths = [fpaths; string(simname(eps,d,l,rho,v,Kagar,Kstiff,rev,runs(f),fp))];
            end
        otherwise
            disp("Nope")
    end
    
    for f = 1:N
        switch dataset
            case "sim"
                files = dir(fullfile(fpaths(f),"cells"));
                ts = 30*(0:120);%(numel(files)-1));
            case "kcwt"
                files = dir(fullfile(fpaths(f),"Laser"));
                ts = getts(fpaths(f));
                ts = ts - ts(1);
            case "mbwt"
                files = dir(fullfile(fpaths(f),"img"));
                T = readtable(fullfile(fpaths(f),'frameTimeZ.csv'));
                ts = T.time-T.time(1);
            case "mbpila"
                files = dir(fullfile(fpaths(f),"img"));
                T = readtable(fullfile(fpaths(f),'frameTimeZ.csv'));
                ts = T.time-T.time(1);
            otherwise
                "Nope"
        end
        files = files(arrayfun(@(x) x.name(1)~='.',files));
        ndefs = [];
    
        for t = 1:121%numel(files)
            fname = fullfile(files(t).folder,files(t).name);
            switch dataset
                case "sim"
                    bds = loadsimdata(fname);
                    adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm);
                    A = boxSize*boxSize;
                case "mbwt"
                    l = imread(fname);
                    holes = findholes(l,9);
                    l = rgb2gray(l);
                    l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                    dirf = dfield_expt(l,rsm);
                    adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
                    sz = size(l);
                    A = (sz(1)*sz(2)-sum(holes,'all'))*XYcal*XYcal;
                case "kcwt"
                    l = laserdata(fname);
                    %holes = findholes(l,XYcal);
                    dirf = dfield_expt(l,rsm);
                    holes = loaddata(fname,'covid_layers','int8')==0;
                    adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
                    sz = size(l);
                    A = (sz(1)*sz(2)-sum(holes,'all'))*XYcal*XYcal;
                case "mbpila"
                    l = imread(fname);
                    holes = findholes(l,9);
                    l = rgb2gray(l);
                    l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                    dirf = dfield_expt(l,rsm);
                    adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
                    sz = size(l);
                    A = (sz(1)*sz(2)-sum(holes,'all'))*XYcal*XYcal;
            end
            ndefs = [ndefs; numel(adefs)/A];
        end
        switch dataset
            case "kcwt"
                col = 'm';
            case "mbwt"
                col = 'r';
            case "mbpila"
                col = 'g';
            case "sim"
                col = 'b';
        end
        plot(ts,ndefs,col)
        hold on
        drawnow
        ndefdatat = struct('ndefs',ndefs,'dataset',dataset,'ts',ts);
        ndefdata = [ndefdata; ndefdatat];
    end
end