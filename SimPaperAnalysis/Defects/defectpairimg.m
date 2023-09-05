%% Experiment image.

XYcal = 0.133;
fpath = '/tigress/kc32/Monolayerpaperdata/Data/190117KC2';

load(fullfile(fpath,'adefs.mat'));
ts = getts(fpath);

for t = 1:numel(ts)
    cdefs = adefs([adefs.ts] == t);

    dx = [cdefs.x] - [cdefs.x]';
    dy = [cdefs.y] - [cdefs.y]';

    dr = sqrt(dx.^2+dy.^2);
    dr(dr==0) = 999;


    if min(dr,[],'all')<50
        img = laserdata(fpath,t);
        imagesc(img./imgaussfilt(img,64));
        colormap gray
        hold on
        plot([cdefs.x],[cdefs.y],'.');
        keyboard;
    end
end

%%
t = 15;

fig = figure();
ax = axes(fig,'Position',[0 0 1 1],'Color','k');
img = laserdata(fpath,t);
dirf = loaddata(fpath,t,'dfield','float');
im = imagesc(dirf);
colormap orientcmap
im.AlphaData = rescale(img./(imgaussfilt(img,64)));
set(gcf,'Color','k')
set(gcf, 'InvertHardCopy', 'off'); 

axis equal
set(gca,'xlim',xlims,'ylim',ylims)
axis off
hold on
plot([xlims(2)-5/0.133 xlims(2)]-5,[ylims(2)-10 ylims(2)-10],'w','LineWidth',2)

%% Simulation image.
eps = 5;
d = 0.7;
l = 7;
rho = 0.3;
v = 5;
Kagar = 500;
Kstiff = 100;
rev = 8;
run = 1;
fp = 'PoissonRevs';

sfpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run,fp);

[files,boxSize] = getframes(sfpath);
res = 752;

sXYcal = boxSize/res;

for t = 1:numel(files)
    
    bds = loadsimdata(fullfile(files(t).folder,files(t).name));
    dirf = dfield_sim(bds,boxSize,res);
    S = nemorderfield_sim(dirf,5);
    defs = finddefects_sim(dirf,S,35,boxSize,res);
    dx = [defs.x] - [defs.x]';
    dy = [defs.y] - [defs.y]';

    dr = sqrt(dx.^2+dy.^2);
    dr(dr==0) = 999;


    if min(dr,[],'all')<10
        img = cellim(bds,boxSize,res,[0 0],1,1,0.5);
        imagesc(img);
        hold on
        colormap gray
        axis equal
        plot([defs.x]/sXYcal,[defs.y]/sXYcal,'.');
        keyboard;
    end
end

%%
simim = cellim(bds,100,res,[0 0], 1,1,0.5);
sXYcal = 100/res;
fig = figure();
ax = axes(fig,'Position',[0 0 1 1],'Color','k');
im = imagesc(dirf);
colormap orientcmap
im.AlphaData = rescale(simim,0.3,0.7);
set(gcf,'Color','k')
set(gcf, 'InvertHardCopy', 'off');

axis equal
sxlims = [0 (xlims(2)-xlims(1))*XYcal/sXYcal]+xoffset;
sylims = [0 (ylims(2)-ylims(1))*XYcal/sXYcal]+yoffset;
set(gca,'xlim',sxlims,'ylim',sylims)
axis off
hold on
plot([sxlims(2)-5/sXYcal sxlims(2)]-5,[sylims(2)-10 sylims(2)-10],'w','LineWidth',2)
