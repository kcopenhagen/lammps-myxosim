tic
dataset = "mbwt"; % kcwt, mbpila, mbwt.
frac = 0.01;
dx = 0.5;
xbins = -20:dx:20;
ybins = -20:dx:20;

switch dataset
    case "sim"
        boxSize = 100;
        minr = 3;   % Minumum defect spacing. (um)
        rsm = 1;    % Director field smoothing for simulated director fields (um).
        XYcal = 0.2; % XY calibration for simulation fields (um / pix).
        dsm = 1;    % Gradient smoothing for calculating defect directions (um).
    case "mbwt"
        minr = 3;
        rsm = 21;
        XYcal = 0.072;
        dsm = 1;
    case "kcwt"
        minr = 3;
        rsm = 13;
        XYcal = 0.133;
        dsm = 1;
    case "mbpila"
        minr = 3;
        rsm = 21;
        XYcal = 0.072;
        dsm = 1;
end

files = getallframes(dataset);

%%
pp = [];
pn = [];
np = [];
nn = [];

for f = 1:numel(files)
    if (rand() < frac)
        fname = files(f);
        switch dataset
            case "sim"
                bds = loadsimdata(fname);
                adefs = finddefs_sim(bds,boxSize,minr,rsm,XYcal,dsm);
            case "mbwt"
                l = imread(fname);
                l = rgb2gray(l);
                l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                holes = false(size(l));
                dirf = dfield_expt(l,rsm);
                adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
            case "kcwt"
                l = laserdata(fname);
                dirf = dfield_expt(l,rsm);
                holes = loaddata(fname,'covid_layers','int8')==0;
                adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
            case "mbpila"
                l = imread(fname);
                l = rgb2gray(l);
                l = imgaussfilt(imsharpen(l,'Radius',3,'Amount',3),2);
                holes = false(size(l));
                dirf = dfield_expt(l,rsm);
                adefs = finddefs_expt(dirf,holes,minr,XYcal,dsm);
        end

        theta0s = atan2([adefs.dy],[adefs.dx]);
        [theta3s, theta0s] = meshgrid(theta0s,theta0s);
        rxs = perdr([adefs.x],[adefs.x],boxSize);
        rys = perdr([adefs.y],[adefs.y],boxSize);
        rrs = sqrt(rxs.^2 + rys.^2);
        theta1s = atan2(rys,rxs);
        theta2s = theta1s - theta0s;
        theta2s(theta2s>pi) = theta2s(theta2s>pi) - 2*pi;
        theta2s(theta2s<-pi) = theta2s(theta2s<-pi) + 2*pi;
        theta4s = theta3s - theta0s;
        theta4s(theta4s>pi) = theta4s(theta4s>pi) - 2*pi;
        theta4s(theta4s<-pi) = theta4s(theta4s<-pi) + 2*pi;
        
        [~,tempqs] = meshgrid([adefs.q],[adefs.q]);
        totq = [adefs.q] + [adefs.q]';
        qdiffs = totq+tempqs;
        


        pp = [pp; struct('r',num2cell(rrs(qdiffs>1)),...
            'theta',num2cell(theta2s(qdiffs>1)),...
            'phi',num2cell(theta4s(qdiffs>1)))];
        pn = [pn; struct('r',num2cell(rrs(abs(qdiffs-0.5)<0.1)),...
            'theta',num2cell(theta2s(abs(qdiffs-0.5)<0.1)),...
            'phi',num2cell(theta4s(abs(qdiffs-0.5)<0.1)))];
        np = [np; struct('r',num2cell(rrs(abs(qdiffs+0.5)<0.1)),...
            'theta',num2cell(theta2s(abs(qdiffs+0.5)<0.1)),...
            'phi',num2cell(theta4s(abs(qdiffs+0.5)<0.1)))];
        nn = [nn; struct('r',num2cell(rrs(qdiffs<-1)),...
            'theta',num2cell(theta2s(qdiffs<-1)),...
            'phi',num2cell(theta4s(qdiffs<-1)))];
    end
end

toc
%% Rotate and duplicate negative defects.



%%
% For density plots.
dr = 1;
rbins = 0:dr:20;
dth = (2*pi/40);
thetabins = -pi:dth:pi;

[rbs, tbs] = meshgrid(rbins',thetabins');
xs = rbs.*cos(tbs);
ys = rbs.*sin(tbs);

%Plot pair interactions.
figure

scatter([pp.r].*cos([pp.theta]),[pp.r].*sin([pp.theta]),3,[pp.phi],'filled')
set(gca,'xlim',[-20 20],'ylim',[-20 20]);
colorcet('R2')
axis equal

figure
rs = discretize([pp.r],rbins);
thetas = discretize([pp.theta],thetabins);
thetas = thetas(~isnan(rs));
rs = rs(~isnan(rs));
ct = ones(size(rs));
cts = accumarray([rs' thetas'],ct,[],@sum);
p = pcolor(xs',ys',padarray(cts,[1 1],1,'post'));

figure
scatter([pn.r].*cos([pn.theta]),[pn.r].*sin([pn.theta]),8,[pn.phi],'filled')
set(gca,'xlim',[-20 20],'ylim',[-20 20]);
colorcet('R2')
axis equal

figure
rs = discretize([pn.r],rbins);
thetas = discretize([pn.theta],thetabins);
thetas = thetas(~isnan(rs));
rs = rs(~isnan(rs));
ct = ones(size(rs));
cts = accumarray([rs' thetas'],ct,[],@sum);
p = pcolor(xs',ys',padarray(cts,[1 1],1,'post'));

figure
scatter([np.r].*cos([np.theta]),[np.r].*sin([np.theta]),8,[np.phi],'filled')
set(gca,'xlim',[-20 20],'ylim',[-20 20]);
colorcet('R2')
axis equal

figure
rs = discretize([np.r],rbins);
thetas = discretize([np.theta],thetabins);
thetas = thetas(~isnan(rs));
rs = rs(~isnan(rs));
ct = ones(size(rs));
cts = accumarray([rs' thetas'],ct,[],@sum);
p = pcolor(xs',ys',padarray(cts,[1 1],1,'post'));
%
figure
scatter([nn.r].*cos([nn.theta]),[nn.r].*sin([nn.theta]),8,[nn.phi],'filled')
set(gca,'xlim',[-20 20],'ylim',[-20 20]);
colorcet('R2')
axis equal

figure
rs = discretize([nn.r],rbins);
thetas = discretize([nn.theta],thetabins);
thetas = thetas(~isnan(rs));
rs = rs(~isnan(rs));
ct = ones(size(rs));
cts = accumarray([rs' thetas'],ct,[],@sum);
p = pcolor(xs',ys',padarray(cts,[1 1],1,'post'));
%%
drx = 0.2;
xbins = -20:drx:20;
ybins = -20:drx:20;

xs = discretize([pp.r].*cos([pp.theta]),xbins);
ys = discretize([pp.r].*sin([pp.theta]),ybins);

ys(isnan(xs)) = [];
xs(isnan(xs)) = [];
xs(isnan(ys)) = [];
ys(isnan(ys)) = [];
ct = ones(size(xs));

cts = accumarray([xs' ys'], ct);
cts(cts>numel(ct)/100) = 0;

figure
imagesc(imgaussfilt(cts,3));

xs = discretize([pn.r].*cos([pn.theta]),xbins);
ys = discretize([pn.r].*sin([pn.theta]),ybins);

ys(isnan(xs)) = [];
xs(isnan(xs)) = [];
xs(isnan(ys)) = [];
ys(isnan(ys)) = [];
ct = ones(size(xs));

cts = accumarray([xs' ys'], ct);
cts(cts>numel(ct)/100) = 0;

figure
imagesc(imgaussfilt(cts,3));

xs = discretize([np.r].*cos([np.theta]),xbins);
ys = discretize([np.r].*sin([np.theta]),ybins);

ys(isnan(xs)) = [];
xs(isnan(xs)) = [];
xs(isnan(ys)) = [];
ys(isnan(ys)) = [];
ct = ones(size(xs));

cts = accumarray([xs' ys'], ct);
%cts(cts>numel(ct)/100) = 0;

figure
imagesc(cts);