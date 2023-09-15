tic
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
frac = 1;
imsz = 40;
dr = 0.2;

xbins = -imsz/2:dr:imsz/2;
ybins = xbins;

pdefvx = zeros(numel(ybins),numel(xbins));
pdefvxM2 = pdefvx;
pdefvy = pdefvx;
pdefvyM2 = pdefvx;
pdefct = pdefvx;
ndefvx = pdefvx;
ndefvxM2 = ndefvx;
ndefvy = pdefvx;
ndefvyM2 = ndefvy;
ndefct = pdefvx;

for r = 1:numel(runs)
    run = runs(r);
    fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run,fp);
    
    [files,boxSize] = getframes(fpath);

    for t = 1:numel(files)
        if (rand() < frac)
            fname = fullfile(files(t).folder,files(t).name);
            bds = loadsimdata(fname);
            adefs = finddefs_sim(bds,boxSize,5,1,0.2,1);
            for i = 1:numel(adefs)
                if (adefs(i).q > 0.1)
                    [cxs,cys,cmux,cmuy,cvx,cvy] = recenter_bds(bds,adefs(i),boxSize,imsz);
                    cxinds = discretize(cxs,xbins);
                    cyinds = discretize(cys,ybins);
                    inds = sub2ind(size(pdefvx),cyinds,cxinds);
                    pdefct(inds) = pdefct(inds) + 1;
                    vxdelta = cvx - pdefvx(inds);
                    vydelta = cvy - pdefvy(inds);
                    pdefvx(inds) = pdefvx(inds) + vxdelta./pdefct(inds);
                    pdefvy(inds) = pdefvy(inds) + vydelta./pdefct(inds);
                    vxdelta2 = cvx - pdefvx(inds);
                    vydelta2 = cvy - pdefvy(inds);
                    pdefvxM2(inds) = pdefvxM2(inds) + vxdelta.*vxdelta2;
                    pdefvyM2(inds) = pdefvyM2(inds) + vydelta.*vydelta2;

                elseif (adefs(i).q < -0.1)
                    cdef = adefs(i);
                    for rot = 0:2*pi/3:(2*pi-0.1)
                        cdefang = atan2(cdef.dy,cdef.dx);
                        cdefang = cdefang + rot;
                        cdef.dx = cos(cdefang);
                        cdef.dy = sin(cdefang);
                        
                        [cxs,cys,cmux,cmuy,cvx,cvy] = recenter_bds(bds,cdef,boxSize,imsz);
                        cxinds = discretize(cxs,xbins);
                        cyinds = discretize(cys,ybins);
                        inds = sub2ind(size(ndefvx),cyinds,cxinds);
                        ndefct(inds) = ndefct(inds) + 1;
                        vxdelta = cvx - ndefvx(inds);
                        vydelta = cvy - ndefvy(inds);
                        ndefvx(inds) = ndefvx(inds) + vxdelta./ndefct(inds);
                        ndefvy(inds) = ndefvy(inds) + vydelta./ndefct(inds);
                        vxdelta2 = cvx - ndefvx(inds);
                        vydelta2 = cvy - ndefvy(inds);
                        ndefvxM2(inds) = ndefvxM2(inds) + vxdelta.*vxdelta2;
                        ndefvyM2(inds) = ndefvyM2(inds) + vydelta.*vydelta2;
                    end
                end
            end
        end
    end
end
toc

%%
qudr = 10;
imsm = 2;
mpdefvx = imgaussfilt(pdefvx,imsm);
mpdefvy = imgaussfilt(pdefvy,imsm);

varpdefvx = imgaussfilt(pdefvxM2./pdefct,imsm);
varpdefvy = imgaussfilt(pdefvyM2./pdefct,imsm);

figure
pvs = sqrt(mpdefvx.^2 + mpdefvy.^2);
pvvar = sqrt(varpdefvx.^2 + varpdefvy.^2);
p = pcolor(xbins,ybins,pvs);
p.EdgeColor = 'none';
colormap(inferno)
axis off
axis equal
hold on

[xx,yy] = meshgrid(xbins,ybins);

quiver(xx(1:qudr:end,1:qudr:end),yy(1:qudr:end,1:qudr:end),...
    mpdefvx(1:qudr:end,1:qudr:end),mpdefvy(1:qudr:end,1:qudr:end),'w','LineWidth',2)
plot([-imsz imsz]/2,[0 0],'g--');
plot([0 0],[-imsz imsz]/2,'g--');

figure
p = pcolor(xbins,ybins,pvvar);
p.EdgeColor = 'none';
colormap(inferno)
axis equal
axis off
set(gca,'clim',[0 7]);

mndefvx = imgaussfilt(ndefvx,imsm);
mndefvy = imgaussfilt(ndefvy,imsm);
varndefvx = imgaussfilt(ndefvxM2./ndefct,imsm);
varndefvy = imgaussfilt(ndefvyM2./ndefct,imsm);

figure
nvs = sqrt(mndefvx.^2 + mndefvy.^2);
nvvar = sqrt(varndefvx.^2 + varndefvy.^2);

n = pcolor(xbins,ybins,nvs);
n.EdgeColor = 'none';

colormap(inferno)
axis off
axis equal
hold on

quiver(xx(1:qudr:end,1:qudr:end),yy(1:qudr:end,1:qudr:end),...
    mndefvx(1:qudr:end,1:qudr:end),mndefvy(1:qudr:end,1:qudr:end),'w','LineWidth',2)
plot([-imsz imsz]/2,[0 0],'g--');
plot([0 0],[-imsz imsz]/2,'g--');

figure
p = pcolor(xbins,ybins,nvvar);
p.EdgeColor = 'none';
colormap(inferno)
axis equal
axis off
set(gca,'clim',[0 7]);

%% Experiments

[fpath1, fpath2, ~] = exptpaths("kcwt");

for r = 1:numel(fpath2)

