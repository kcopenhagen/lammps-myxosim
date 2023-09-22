%% Simulation flow fields.
minr = 5;   % Minumum defect spacing. (um)
rsm = 1;    % Director field smoothing for simulated director fields (um).
XYcal = .2; % XY calibration for simulation fields (um / pix).
dsm = 1;    % Gradient smoothing for calculating defect directions (um).

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
            adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm);
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
tic
[fpath1, fpath2, fpath3] = exptpaths("mbpila");
rsm = 5;
XYcal = 0.072;
dsm = 1;
minr = 5;
imsz = 40;
dr = 0.2;

xbins = -imsz/2:dr:imsz/2;
ybins = xbins;

epdefvx = zeros(numel(ybins),numel(xbins));
epdefvxM2 = epdefvx;
epdefvy = epdefvx;
epdefvyM2 = epdefvx;
epdefct = epdefvx;
endefvx = epdefvx;
endefvxM2 = endefvx;
endefvy = epdefvx;
endefvyM2 = endefvy;
endefct = epdefvx;

for r = 1:numel(fpath3)
    
    fpath = fullfile(fpath1,fpath2,fpath3(r));
    fpath = fullfile(fpath,'img');
    files = dir(fpath);
    files = files(~[files.isdir]);

    for t = 2:numel(files)
        l = laserdata(fpath,t);
        dirf = dfield_expt(l,rsm);
        holes = (loaddata(fpath,t,'covid_layers','int8') == 0);
        holes(1:10,:) = true;
        holes(:,1:10) = true;
        holes(end-9:end,:) = true;
        holes(:,end-9:end) = true;
        holes = imdilate(holes,strel('disk',6));
        adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);

        Vx = loaddata(fpath,t,'flows/Vx','float')*60;
        Vy = loaddata(fpath,t,'flows/Vy','float')*60;

        for i = 1:numel(adefs)
            % if (adefs(i).q > 0.1)
            %     [cxs,cys,cvx,cvy] = recenter_eVs(Vx,Vy,adefs(i),XYcal,imsz);
            %     cxinds = discretize(cxs,xbins);
            %     cyinds = discretize(cys,ybins);
            %     inds = sub2ind(size(epdefvx),cyinds,cxinds);
            %     epdefct(inds) = epdefct(inds) + 1;
            %     vxdelta = cvx - epdefvx(inds);
            %     vydelta = cvy - epdefvy(inds);
            %     epdefvx(inds) = epdefvx(inds) + vxdelta./epdefct(inds);
            %     epdefvy(inds) = epdefvy(inds) + vydelta./epdefct(inds);
            %     vxdelta2 = cvx - epdefvx(inds);
            %     vydelta2 = cvy - epdefvy(inds);
            %     epdefvxM2(inds) = epdefvxM2(inds) + vxdelta.*vxdelta2;
            %     epdefvyM2(inds) = epdefvyM2(inds) + vydelta.*vydelta2;
            if (adefs(i).q < -0.1)
                for rot = 0:2*pi/3:(2*pi-0.1)
                    cdef = adefs(i);
                    cdefang = atan2(cdef.dy,cdef.dx);
                    cdefang = cdefang + rot;
                    cdef.dx = cos(cdefang);
                    cdef.dy = sin(cdefang);
                    
                    [cxs,cys,cvx,cvy] = recenter_eVs(Vx,Vy,cdef,XYcal,imsz);
                    cxinds = discretize(cxs,xbins);
                    cyinds = discretize(cys,ybins);
                    inds = sub2ind(size(endefvx),cyinds,cxinds);
                    endefct(inds) = endefct(inds) + 1;
                    vxdelta = cvx - endefvx(inds);
                    vydelta = cvy - endefvy(inds);
                    endefvx(inds) = endefvx(inds) + vxdelta./endefct(inds);
                    endefvy(inds) = endefvy(inds) + vydelta./endefct(inds);
                    vxdelta2 = cvx - endefvx(inds);
                    vydelta2 = cvy - endefvy(inds);
                    endefvxM2(inds) = endefvxM2(inds) + vxdelta.*vxdelta2;
                    endefvyM2(inds) = endefvyM2(inds) + vydelta.*vydelta2;
                end
            end
        end
    end
end
toc

%% Experiment images

qudr = 10;
imsm = 2;
empdefvx = imgaussfilt(epdefvx,imsm);
empdefvy = imgaussfilt(epdefvy,imsm);

evarpdefvx = imgaussfilt(epdefvxM2./epdefct,imsm);
evarpdefvy = imgaussfilt(epdefvyM2./epdefct,imsm);

figure
epvs = sqrt(empdefvx.^2 + empdefvy.^2);
epvvar = sqrt(evarpdefvx.^2 + evarpdefvy.^2);
p = pcolor(xbins,ybins,epvs);
p.EdgeColor = 'none';
colormap(summer)
axis off
axis equal
hold on

[xx,yy] = meshgrid(xbins,ybins);

quiver(xx(1:qudr:end,1:qudr:end),yy(1:qudr:end,1:qudr:end),...
    empdefvx(1:qudr:end,1:qudr:end),empdefvy(1:qudr:end,1:qudr:end),'k','LineWidth',2)
plot([-imsz imsz]/2,[0 0],'g--');
plot([0 0],[-imsz imsz]/2,'g--');

figure
p = pcolor(xbins,ybins,epvvar);
p.EdgeColor = 'none';
colormap(inferno)
axis equal
axis off
set(gca,'clim',[0 7]);

emndefvx = imgaussfilt(endefvx,imsm);
emndefvy = imgaussfilt(endefvy,imsm);
evarndefvx = imgaussfilt(endefvxM2./endefct,imsm);
evarndefvy = imgaussfilt(endefvyM2./endefct,imsm);

figure
envs = sqrt(emndefvx.^2 + emndefvy.^2);
envvar = sqrt(evarndefvx.^2 + evarndefvy.^2);

n = pcolor(xbins,ybins,envs);
n.EdgeColor = 'none';

colormap(inferno)
axis off
axis equal
hold on

quiver(xx(1:qudr:end,1:qudr:end),yy(1:qudr:end,1:qudr:end),...
    emndefvx(1:qudr:end,1:qudr:end),emndefvy(1:qudr:end,1:qudr:end),'w','LineWidth',2)
plot([-imsz imsz]/2,[0 0],'g--');
plot([0 0],[-imsz imsz]/2,'g--');

figure
p = pcolor(xbins,ybins,envvar);
p.EdgeColor = 'none';
colormap(inferno)
axis equal
axis off
set(gca,'clim',[0 7]);


%% Saving
fname = fullfile('~/lammps-myxosim/SimPaperAnalysis/Figures/Data','flows.mat');

save(fname,'pdefvx_kcwt',epdefvx,'pdefvy_kcwt',epdefvy,...
    'pdefvxM2_kcwt',epdefvxM2,'pdefvyM2_kcwt',epdefvyM2,'pdefct_kcwt',epdefct,...
    'ndefvx_kcwt',endefvx,'ndefvy_kcwt',endefvy,...
    'ndefvxM2_kcwt',endefvxM2,'ndefvyM2_kcwt',endefvyM2,'ndefct_kcwt',endefct,...
    'pdefvx_sim',pdefvx,'pdefvy_sim',pdefvy,...
    'pdefvxM2_sim',pdefvxM2,'pdefvyM2_sim',pdefvyM2,'pdefct_sim',pdefct,...
    'ndefvx_sim',ndefvx,'ndefvy_sim',ndefvy,...
    'ndefvxM2_sim',ndefvxM2,'ndefvyM2_sim',ndefvyM2,'ndefct_sim',ndefct);

