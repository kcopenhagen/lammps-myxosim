tic
disp(dataset)         % kcwt, mbpila, mbwt, sim.
frac = 1;
dx = 0.5;

xbins = -20:dx:20;
ybins = -20:dx:20;

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

files = getallframes(dataset);

%%

ppmeanphisx = zeros(numel(ybins),numel(xbins));
ppM2phisx = ppmeanphisx;
ppctx = ppmeanphisx;
pnmeanphisx = ppmeanphisx;
pnM2phisx = ppmeanphisx;
pnctx = ppmeanphisx;
npmeanphisx = ppmeanphisx;
npM2phisx = ppmeanphisx;
npctx = ppmeanphisx;
nnmeanphisx = ppmeanphisx;
nnM2phisx = ppmeanphisx;
nnctx = ppmeanphisx;

ppmeanphisy = zeros(numel(ybins),numel(xbins));
ppM2phisy = ppmeanphisy;
ppcty = ppmeanphisy;
pnmeanphisy = ppmeanphisy;
pnM2phisy = ppmeanphisy;
pncty = ppmeanphisy;
npmeanphisy = ppmeanphisy;
npM2phisy = ppmeanphisy;
npcty = ppmeanphisy;
nnmeanphisy = ppmeanphisy;
nnM2phisy = ppmeanphisy;
nncty = ppmeanphisy;

for f = 1:numel(files)
    if (rand() < frac)
        fname = files(f);
        switch dataset
            case "sim"
                bds = loadsimdata(fname);
                adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm);
            case "mbwt"
                l = imread(fname);
                holes = findholes(l,9);
                l = rgb2gray(l);
                l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                dirf = dfield_expt(l,rsm);
                adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
            case "kcwt"
                l = laserdata(fname);
                %holes = findholes(l,XYcal);
                dirf = dfield_expt(l,rsm);
                holes = loaddata(fname,'covid_layers','int8')==0;
                adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
            case "mbpila"
                l = imread(fname);
                holes = findholes(l,9);
                l = rgb2gray(l);
                l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                dirf = dfield_expt(l,rsm);
                adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
        end

        for i = 1:numel(adefs)
            for j = 1:numel(adefs)
                if (i ~= j)
                    phii = atan2(adefs(i).dy,adefs(i).dx);
                    phij = atan2(adefs(j).dy,adefs(j).dx);
                    
                    switch dataset
                        case "sim"
                            rx = perdr(adefs(j).x,adefs(i).x,boxSize);
                            ry = perdr(adefs(j).y,adefs(i).y,boxSize);
                        otherwise
                            rx = adefs(j).x - adefs(i).x;
                            ry = adefs(j).y - adefs(i).y;
                    end

                    r = sqrt(rx.^2 + ry.^2);
                    labtheta = atan2(ry,rx);


                    if (adefs(i).q > 0.1 && adefs(j).q > 0.1)
                        theta = labtheta - phii;
                        while (theta>pi)
                            theta = theta - 2*pi;
                        end
                        while (theta<-pi)
                            theta = theta + 2*pi;
                        end
                        phi = phij - phii;
                        while (phi>pi)
                            phi = phi - 2*pi;
                        end
                        while (phi<-pi)
                            phi = phi + 2*pi;
                        end
        
                        xind = discretize(r*cos(theta),xbins);
                        yind = discretize(r*sin(theta),ybins);

                        if ~isnan(xind) && ~isnan(yind)
                            [ppmeanphisx(yind,xind), ppM2phisx(yind,xind),...
                                ppctx(yind,xind)] = ...
                                running_mean(ppmeanphisx(yind,xind),...
                                ppM2phisx(yind,xind),ppctx(yind,xind),cos(phi));
                            [ppmeanphisy(yind,xind), ppM2phisy(yind,xind),...
                                ppcty(yind,xind)] = ...
                                running_mean(ppmeanphisy(yind,xind),...
                                ppM2phisy(yind,xind),ppcty(yind,xind),sin(phi));
                        end
                        if sum(isnan(ppmeanphisx),'all')>0
                            keyboard
                        end

                    elseif (adefs(i).q > 0.1 && adefs(j).q < -0.1)
                        theta = labtheta - phii;
                        while (theta>pi)
                            theta = theta - 2*pi;
                        end
                        while (theta<-pi)
                            theta = theta + 2*pi;
                        end

                        phi = phij - phii;
                        phi = phi*3;

                        while (phi>pi)
                            phi = phi - 2*pi;
                        end
                        while (phi<-pi)
                            phi = phi + 2*pi;
                        end

                        xind = discretize(r*cos(theta),xbins);
                        yind = discretize(r*sin(theta),ybins);

                        if ~isnan(xind) && ~isnan(yind)
                            [pnmeanphisx(yind,xind), pnM2phisx(yind,xind),...
                                pnctx(yind,xind)] = ...
                                running_mean(pnmeanphisx(yind,xind),...
                                pnM2phisx(yind,xind),pnctx(yind,xind),cos(phi));
                            [pnmeanphisy(yind,xind), pnM2phisy(yind,xind),...
                                pncty(yind,xind)] = ...
                                running_mean(pnmeanphisy(yind,xind),...
                                pnM2phisy(yind,xind),pncty(yind,xind),sin(phi));
                        end
                    
                    elseif (adefs(i).q < -0.1 && adefs(j).q > 0.1)
                        for rot = 1:3
                            phii = phii + 2*pi/3;
                            theta = labtheta - phii;
                            while (theta>pi)
                                theta = theta - 2*pi;
                            end
                            while (theta<-pi)
                                theta = theta + 2*pi;
                            end
    
                            phi = phij - phii;
                            while (phi>pi)
                                phi = phi - 2*pi;
                            end
                            while (phi<-pi)
                                phi = phi + 2*pi;
                            end
    
                            xind = discretize(r*cos(theta),xbins);
                            yind = discretize(r*sin(theta),ybins);
    
                            if ~isnan(xind) && ~isnan(yind)
                                [npmeanphisx(yind,xind), npM2phisx(yind,xind),...
                                    npctx(yind,xind)] = ...
                                    running_mean(npmeanphisx(yind,xind),...
                                    npM2phisx(yind,xind),npctx(yind,xind),cos(phi));
                                [npmeanphisy(yind,xind), npM2phisy(yind,xind),...
                                    npcty(yind,xind)] = ...
                                    running_mean(npmeanphisy(yind,xind),...
                                    npM2phisy(yind,xind),npcty(yind,xind),sin(phi));
                            end
                        end
                        
                    elseif (adefs(i).q < -0.1 && adefs(j).q < -0.1)
                        for rot = 1:3
                            phii = phii+2*pi/3;
                            
                            theta = labtheta - phii;
                            while (theta>pi)
                                theta = theta - 2*pi;
                            end
                            while (theta<-pi)
                                theta = theta + 2*pi;
                            end
    
                            phi = phij - phii;
                            phi = phi*3;
                            while (phi>pi)
                                phi = phi - 2*pi;
                            end
                            while (phi<-pi)
                                phi = phi + 2*pi;
                            end
    
                            xind = discretize(r*cos(theta),xbins);
                            yind = discretize(r*sin(theta),ybins);
    
                            if ~isnan(xind) && ~isnan(yind)
                                [nnmeanphisx(yind,xind), nnM2phisx(yind,xind),...
                                    nnctx(yind,xind)] = ...
                                    running_mean(nnmeanphisx(yind,xind),...
                                    nnM2phisx(yind,xind),nnctx(yind,xind),cos(phi));
                                [nnmeanphisy(yind,xind), nnM2phisy(yind,xind),...
                                    nncty(yind,xind)] = ...
                                    running_mean(nnmeanphisy(yind,xind),...
                                    nnM2phisy(yind,xind),nncty(yind,xind),sin(phi));
                            end
                        end
                        
                    end
                end
            end
        end
    end
end

toc
%% Plots

pos = [0 0 1000 500];
fig = figure('Position',pos);
ax1 = axes(fig,'Position',[0 0 0.25 0.5]);
ax2 = axes(fig,'Position',[0.25 0 0.25 0.5]);
ax3 = axes(fig,'Position',[0.5 0 0.25 0.5]);
ax4 = axes(fig,'Position',[0.75 0 0.25 0.5]);
ax5 = axes(fig,'Position',[0 0.5 0.25 0.5]);
ax6 = axes(fig,'Position',[0.25 0.5 0.25 0.5]);
ax7 = axes(fig,'Position',[0.5 0.5 0.25 0.5]);
ax8 = axes(fig,'Position',[0.75 0.5 0.25 0.5]);

imagesc(ax1,ppctx)
imagesc(ax2,pnctx)
imagesc(ax3,npctx)
imagesc(ax4,nnctx)

imagesc(ax5,atan2(ppmeanphisy,ppmeanphisx))
imagesc(ax6,atan2(pnmeanphisy,pnmeanphisx)/3)
imagesc(ax7,atan2(npmeanphisy,npmeanphisx))
imagesc(ax8,atan2(nnmeanphisy,nnmeanphisx)/3)

axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8];

for i = 5:8
    colormap(axs(i),colorcet('C8'));
end
for i = 1:8
    set(axs(i),'Visible','off');
end

set(ax5,'clim',[-pi pi]);
set(ax6,'clim',[-pi/3 pi/3]);
set(ax7,'clim',[-pi pi]);
set(ax8,'clim',[-pi/3 pi/3]);

%%
save(fullfile('~/lammps-myxosim/SimPaperAnalysis/Figures/Defects/allpairs/',dataset),...
    'ppctx','ppmeanphisx','ppmeanphisy','ppM2phisx','ppM2phisy',...
    'pnctx','pnmeanphisx','pnmeanphisy','pnM2phisx','pnM2phisy',...
    'npctx','npmeanphisx','npmeanphisy','npM2phisx','npM2phisy',...
    'nnctx','nnmeanphisx','nnmeanphisy','nnM2phisx','nnM2phisy');