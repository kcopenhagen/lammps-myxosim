%% My WT.

XYcal = 0.072;
fpath = '/tigress/kc32/train/msk';
files = dir(fpath);
files = files(~[files.isdir]);
gd = arrayfun(@(x) x.name(1) ~= '.',files);
files = files(gd);
kcdens = [];

for t = 1:numel(files)

    im = load(fullfile(files(t).folder,files(t).name));
    cmsks = im.msk2cls == 1;
    CC = bwconncomp(cmsks);
    R = regionprops(CC,'PixelIdxList','Area');
    R = R([R.Area]>60);
    kcdens = [kcdens; numel(R)/(size(cmsks,1)*size(cmsks,2)*XYcal*XYcal)];
end

%% Matt's pilA
fpaths = '/scratch/gpfs/kc32/matt-monolayers/mb-pila/2022-06-23/';
runs = ["trial1"; "trial2"; "trial3"; "trial4"];

mbpiladens = [];
XYcal = 0.072;

for r = 1:numel(runs)
%%
    fpath = fullfile(fpaths,runs(r),'seg');
    files = dir(fpath);
    files = files(~[files.isdir]);
    testfiles = files(rand(size(files))<0.25);
%%
    for fi = 1:numel(testfiles)
        %%
        fname = fullfile(testfiles(fi).folder,testfiles(fi).name);
        [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(fname);
        
        f = memmapfile(fname, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
        
        
        tmp = f.Data.d;
        tmpt = tmp(:,:,1).*(1-tmp(:,:,2));
        cmsks = tmpt>0.1;
        CC = bwconncomp(cmsks);
        Rt = regionprops(CC,'PixelIdxList','Area');
        R = Rt([Rt.Area]>42);
        %%
        mbpiladens = [mbpiladens; numel(R)/(size(cmsks,1)*size(cmsks,2)*XYcal*XYcal)];
    end
end

%% Matt's WT.

fpaths = '/scratch/gpfs/kc32/matt-monolayers/mb-wt/2022-11-20/';
runs = ["trial1"; "trial2"; "trial3"];

mbwtdens = [];
XYcal = 0.076;

for r = 1:numel(runs)
    fpath = fullfile(fpaths,runs(r),'seg');
    files = dir(fpath);
    files = files(~[files.isdir]);
    testfiles = files(rand(size(files))<0.25);
    for fi = 1:numel(testfiles)
        fname = fullfile(testfiles(fi).folder,testfiles(fi).name);
        [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(fname);
        
        f = memmapfile(fname, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);
        
        %%
        tmp = f.Data.d;
        cmsks = tmp>0.5;
        CC = bwconncomp(cmsks);
        R = regionprops(CC,'PixelIdxList','Area');
        R = R([R.Area]>60);
        
        mbwtdens = [mbwtdens; numel(R)/(size(cmsks,1)*size(cmsks,2)*XYcal*XYcal)];
    end
end

%% Simulations.

eps = 5;
d = 0.7;
l = 7;
rho = 0.21;
v = 5;
Kagar = 200;
Kstiff = 15;
rev = 8;
runs = 0:3;
props = [];

simdens = [];
for r = 1:numel(runs)
    run = runs(r);
    fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run);
    [files, boxSize] = getframes(fpath);
    testfiles = files(rand(size(files))<0.1);
    t = 1;
    %for t = 1:numel(testfiles)
        bds = loadsimdata(fullfile(testfiles(t).folder,testfiles(t).name));
        aids = unique([bds.mol]);
        simdens = [simdens; numel(aids)/(boxSize*boxSize)];
    %end
end
%%
save(fullfile('~/SimPaperAnalysis/Figures/Data/'))

%%
eds = 0:0.01:0.3;
xs = (eds(2) - eds(1))/2 + eds(1):(eds(2)-eds(1)):eds(end);

mbwt = histcounts(mbwtdens, eds, 'Normalization','pdf');
mbpila = histcounts(mbpiladens, eds, 'Normalization','pdf');
kc = histcounts(kcdens, eds, 'Normalization', 'pdf');
sim = histcounts(simdens, eds, 'Normalization', 'pdf');

plot(xs,mbwt,'r','LineWidth',2);
hold on
plot(xs,mbpila,'g','LineWidth',2);
plot(xs,kc,'m','LineWidth',2);
plot(xs,sim,'b','LineWidth',2);
set(gca,'FontSize',24);
