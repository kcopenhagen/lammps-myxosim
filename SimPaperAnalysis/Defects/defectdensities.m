%% Defect count

%% Experiments
exptnames = ["190109KC3", "190109KC3b", "190111KC2", "190116KC2", "190117KC1", "190117KC2", "190117KC2a", "190117KC2b"]; 
eds = 0.0001:0.0002:0.004;
andefs = [];
for f = 1:numel(exptnames)
    fpath = fullfile('/tigress/kc32/Monolayerpaperdata/Data',exptnames(f));
    load(fullfile(fpath,'adefs.mat'));
    ts = getts(fpath);
    ndefs = [];
    for t = 1:numel(ts)
        ndefs = [ndefs; sum([adefs.ts]==t)];
    end
    A = 1024*768*0.133*0.133;
    histogram(ndefs/A,eds)
    hold on
    drawnow
    andefs = [andefs; ndefs];
end

%% Simulations

eps = 5;
d = 0.7;
l = 7;
rho = 0.3;
v = 5;
Kagar = 500;
Kstiff = 100;
revs = [1 2 4 6 8 10 12 14 16 18 20 100 999];
runs = 0:2;
res = 1389;
props = [];

for r = 1:numel(revs)
    ndefs = [];
    for r2 = 1:numel(runs)
        run = runs(r2)
        fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,revs(r),run,'PoissonRevs');
        [files,boxSize] = getframes(fpath);
        files = files(rand(size(files))<0.05);
        for t = 1:numel(files)
            bds = loadsimdata(fullfile(files(t).folder,files(t).name));
            dirf = dfield_sim(bds,boxSize,res);
            S = nemorderfield_sim(dirf,11);
            defs = finddefects_sim(dirf,S,35,boxSize,res);
            ndefs = [ndefs; numel(defs)];
        end
    end
    A = 100*100;
    propt = struct('rev',revs(r),'defdens',ndefs/A);
    props = [props; propt];
end

%%
meany = [];
stdy = [];
xs = [];
for i = 1:numel(props)
    xs = [xs; props(i).rev];
    meany = [meany; mean(props(i).defdens)];
    stdy = [stdy; std(props(i).defdens)];
end

errorbar(xs,meany,stdy,'b','LineWidth',2);
hold on
plot([xs(1) xs(end)],[mean(andefs/A) mean(andefs/A)],'m','LineWidth',2);
plot([xs(1) xs(end)],(mean(andefs/A)-std(andefs/A))*[1 1],'m--','LineWidth',2);
plot([xs(1) xs(end)],(mean(andefs/A)+std(andefs/A))*[1 1],'m--','LineWidth',2);
set(gca,'FontSize',24);