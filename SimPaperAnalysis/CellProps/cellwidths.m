%% My WT data.

test_n = 100;

[fpath1, fpath2, ~] = exptpaths("kcwt");
XYcal = 0.133;
acellw = [];

for f = 1:numel(fpath2)
    fpath = fullfile(fpath1,fpath2(f));
    ts = getts(fpath);
    cellw = [];
    for test_n = 1:numel(ts)
        % t = randi(numel(ts));
        t = test_n;
        img = laserdata(fpath,t);
        img = img./imgaussfilt(img,64);
        proc_img = img;
        proc_img(img>1.05) = 1.05;
        proc_img(img<0.95) = 0.95;
        proc_img = rescale(proc_img,0, 1);
        sz = size(proc_img);
        n = min(sz);
        F = fft2(proc_img,n,n);
        
        ns = -n/2:n/2-1;
        Ft = abs(fftshift(F));
        % imagesc(Ft,'XData',ns,'YData',ns);
        % set(gca,'clim',[0 10000])
        % axis equal
        
        [xx,yy] = meshgrid(ns,ns);
        rr = sqrt(xx.^2 + yy.^2);
        inds = round(rr)+1;
        
        rs = accumarray(inds(:),rr(:),[],@mean);
        ffs = accumarray(inds(:),Ft(:),[],@mean);
        
        minprom = 1000;
        [~,loc] = findpeaks(ffs,rs,'MinPeakProminence',minprom,'MinPeakWidth',15);
        while isempty(loc) && minprom>50
            minprom = minprom - 50;
            [~,loc] = findpeaks(ffs,rs,'MinPeakProminence',minprom,'MinPeakWidth',15);
        end
        cellw = [cellw; n/loc*XYcal];
    end
    acellw = [acellw; cellw];
    % histogram(cellw,eds)
    % hold on
    % drawnow
end

kcwtcw = acellw;
%% Matt wt.

[fpath1, fpath2, fpath3] = exptpaths("mbwt");
tests_n = 100;
XYcal = 0.072;
eds = 0.6:0.01:0.8;
acellw = [];
for f = 1:numel(fpath3)
    fpath = fullfile(fpath1,fpath2,fpath3s(f));
    frames = dir(fullfile(fpath,'img'));
    frames = frames(~[frames.isdir]);
    cellw = [];
    for test_n = 1:numel(frames)
        % t = randi(numel(frames));
        t = test_n;
        fname = fullfile(frames(t).folder,frames(t).name);
        img = loadrgbimg(fname,'true');
        proc_img = img;
        proc_img(img>1.05) = 1.05;
        proc_img(img<0.95) = 0.95;
        proc_img = rescale(proc_img,0, 1);
        sz = size(proc_img);
        n = min(sz);
        F = fft2(proc_img,n,n);
        
        ns = -n/2:n/2-1;
        Ft = abs(fftshift(F));
        % imagesc(Ft,'XData',ns,'YData',ns);
        % set(gca,'clim',[0 10000])
        % axis equal
        
        [xx,yy] = meshgrid(ns,ns);
        rr = sqrt(xx.^2 + yy.^2);
        inds = round(rr)+1;
        
        rs = accumarray(inds(:),rr(:),[],@mean);
        ffs = accumarray(inds(:),Ft(:),[],@mean);
        
        [~,loc] = findpeaks(ffs,rs,'MinPeakProminence',500,'MinPeakWidth',20);
        
        cellw = [cellw; n/loc*XYcal];
    end
    acellw = [acellw; cellw];
    % histogram(cellw,eds)
    % hold on
    % drawnow
end

mbwtcw = acellw;

%% Matt pilA.
[fpath1, fpath2, fpath3] = exptpaths("mbpila");
tests_n = 100;
XYcal = 0.076;
eds = 0.6:0.01:0.8;
acellw = [];
for f = 1:numel(fpath3)
    fpath = fullfile(fpath1,fpath2,fpath3(f));
    frames = dir(fullfile(fpath,'img'));
    frames = frames(~[frames.isdir]);
    cellw = [];
    for test_n = 1:tests_n
        % t = randi(numel(frames));
        t = test_n;
        fname = fullfile(frames(t).folder,frames(t).name);
        img = loadrgbimg(fname,'true');
        proc_img = img;
        proc_img(img>1.05) = 1.05;
        proc_img(img<0.95) = 0.95;
        proc_img = rescale(proc_img,0, 1);
        sz = size(proc_img);
        n = min(sz);
        F = fft2(proc_img,n,n);
        
        ns = -n/2:n/2-1;
        Ft = abs(fftshift(F));
        % imagesc(Ft,'XData',ns,'YData',ns);
        % set(gca,'clim',[0 10000])
        % axis equal
        
        [xx,yy] = meshgrid(ns,ns);
        rr = sqrt(xx.^2 + yy.^2);
        inds = round(rr)+1;
        
        rs = accumarray(inds(:),rr(:),[],@mean);
        ffs = accumarray(inds(:),Ft(:),[],@mean);
        
        [~,loc] = findpeaks(ffs,rs,'MinPeakProminence',500,'MinPeakWidth',20);
        
        cellw = [cellw; n/loc*XYcal];
    end
    acellw = [acellw; cellw];
    % histogram(cellw,eds)
    % hold on
    % drawnow
end

mbpilcw = acellw;
%%  Simulations
sXYcal = 0.1;
tic
eps = 5;
d = 0.7;
l = 7;
rho = 0.21;
v = 5;
Kagar = 250;
Kstiff = 33;
rev = 8;
run = 0;
props = [];
nts = 100;

for e = 1:numel(eps)
    fpath = simname(eps(e), d, l, rho, v, Kagar, Kstiff, rev, run);

    [files, boxSize] = getframes(fpath);
    cws = [];
    for ni = 1:nts
        ctime = randi(numel(files));
        bds = loadsimdata(fullfile(files(ctime).folder,files(ctime).name));
        cim = cellim(bds,boxSize,1000);
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
        cws = [cws; n/loc*sXYcal];
    end
    propt = struct('eps',eps(e),'cws',cws);
    props = [props; propt];
end
toc
%%

save(fullfile('~/SimPaperAnalysis','Figures','Data','cellwidths.mat'),...
    'kcwtcw','mbwtcw','mbpilcw','props');

%%

eds = 0.6:0.02:0.9;
xs = (eds(1)+(eds(2)-eds(1))/2):(eds(2)-eds(1)):eds(end);

ys = histcounts(kcwtcw,eds,'Normalization','pdf');
plot(xs,ys,'m','LineWidth',2);
hold on
ys = histcounts(mbpilcw*0.072/0.076,eds,'Normalization','pdf');
plot(xs,ys,'g','LineWidth',2);
ys = histcounts(mbwtcw*0.072/0.076,eds,'Normalization','pdf');
plot(xs,ys,'r','LineWidth',2);
simcw = props(1).cws;
ys = histcounts(simcw,eds,'Normalization','pdf');
plot(xs,ys,'b','LineWidth',2);

set(gca,'FontSize',24);
xlabel('Cell widths ({\mu}m)');
ylabel('PDF');