%% Simulation defect interactions.
minr = 5;   % Minumum defect spacing. (um)
rsm = 1;    % Director field smoothing for simulated director fields (um).
XYcal = 0.2; % XY calibration for simulation fields (um / pix).
dsm = 1;    % Gradient smoothing for calculating defect directions (um).

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

drpp = [];
drpn = [];
drnn = [];

for r = 1:numel(runs)
    fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,runs(r),fp);
    [files,boxSize] = getframes(fpath);
    for t = 1:numel(files)
        fname = fullfile(files(t).folder,files(t).name);
        bds = loadsimdata(fname);
        adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm);

        pdefs = adefs([adefs.q]>0.1);
        dxs = perdr([pdefs.x],[pdefs.x],boxSize);
        dys = perdr([pdefs.y],[pdefs.y],boxSize);
        drs = sqrt(dxs.^2 + dys.^2);
        temp = tril(drs);
        drpp = [drpp; temp(temp>0.01)];
        
        ndefs = adefs([adefs.q]<-0.1);
        dxs = perdr([ndefs.x],[ndefs.x],boxSize);
        dys = perdr([ndefs.y],[ndefs.y],boxSize);
        drs = sqrt(dxs.^2 + dys.^2);
        temp = tril(drs);
        drnn = [drnn; temp(temp>0.01)];

        dxs = perdr([pdefs.x],[ndefs.x],boxSize);
        dys = perdr([pdefs.y],[pdefs.x],boxSize);
        drs = sqrt(dxs.^2 + dys.^2);
        temp = tril(drs);
        drpn = [drpn; temp(temp>0.01)];

    end
end

%%

binsz = 0.5;

drbins = 0:binsz:50;
xs = (drbins(1)+binsz/2):binsz:drbins(end);

ppgr = histcounts(drpp,drbins)./xs;
pngr = histcounts(drpn,drbins)./xs;
nngr = histcounts(drnn,drbins)./xs;

plot(-log(ppgr),'r','LineWidth',2);
hold on
plot(-log(pngr),'Color',[1 0.5 1],'LineWidth',2)
plot(-log(nngr),'b','LineWidth',2);

