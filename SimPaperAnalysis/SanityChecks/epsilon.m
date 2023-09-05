
epss = 1:2:21;
d = 0.7;
l = 7;
rho = 0.21;
v = 5;
Kagar = 300;
Kstiff = 15;
rev = 8;
runs = 0:2;
props = [];
sXYcal = 0.1;

for e = 1:numel(epss)
    eps = epss(e);
    runmeanv = 0;
    runstdv = 0;
    ctv = 0;

    runmeancw = 0;
    runstdcw = 0;
    ctcw = 0;

    for r = 1:numel(runs)
        run = runs(r);
        fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run);
        [files,boxSize] = getframes(fpath);
        
        for t = 1:10:numel(files)
            bds = loadsimdata(fullfile(files(t).folder,files(t).name));
            if (rand() < 0.1)
                spds = sqrt([bds.vx].^2 + [bds.vy].^2 + [bds.vz].^2);
                for i = 1:numel(spds)
                    ctv = ctv + 1;
                    prevmeanv = runmeanv;
                    runmeanv = runmeanv + (spds(i)-runmeanv)/ctv;
                    runstdv = runstdv + (spds(i)-runmeanv)*(spds(i)-prevmeanv);
                end
            end

            cim = cellim(bds,100,1000);
            sz = size(cim);
            n = min(sz);
            F = fft2(cim,n,n);

            ns = -n/2:n/2-1;
            Ft = abs(fftshift(F));
            
            [xx,yy] = meshgrid(ns,ns);
            rr = sqrt(xx.^2 + yy.^2);
            inds = round(rr/4)+1;
            
            rs = accumarray(inds(:),rr(:),[],@mean);
            ffs = accumarray(inds(:),Ft(:),[],@mean);
            
            minprom = 1000;
            [~,loc] = findpeaks(ffs,rs,'MinPeakProminence',minprom,'MinPeakWidth',5);
            while isempty(loc) && minprom>50
                minprom = minprom - 50;
                [~,loc] = findpeaks(ffs,rs,'MinPeakProminence',minprom,'MinPeakWidth',5);
            end

            cw = n/loc*sXYcal;
            ctcw = ctcw + 1;
            prevmeancw = runmeancw;
            runmeancw = runmeancw + (cw-runmeancw)/ctcw;
            runstdcw = runstdcw + (cw-runmeancw)*(cw-prevmeancw);
        end
    end
    propt = struct('eps',eps,'meanv',runmeanv,'stdv',sqrt(runstdv/ctv),...
        'meancw',runmeancw,'stdcw',sqrt(runstdcw/ctcw));
    props = [props; propt];
end

%%

errorbar([props.eps],[props.meanv],[props.stdv],'LineWidth',2);
set(gca,'FontSize',24)
figure
errorbar([props.eps],[props.meancw],[props.stdcw],'LineWidth',2);
set(gca,'FontSize',24)

