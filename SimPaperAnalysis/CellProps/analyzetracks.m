%% Velocities

fnames = ["first.mat","third.mat","fourth.mat","fifth.mat","sixth.mat",...
    "seventh.mat","eighth.mat"];
  
avs = [];
for f = 2:numel(fnames)
load(fnames(f),'tracks','drifts','fpath');

%%

eds = 0:20;
xs = 0.5:20;
ts = getts(fpath);
XYcal = 0.133;

drx = diff(drifts(:,1))*XYcal;
dry = diff(drifts(:,2))*XYcal;
drx = [0; drx];
dry = [0; dry];
cvs = [];
for tr = 1:numel(tracks)

    [~,sinds] = sort(tracks(tr).t);
    % xcorr = tracks(tr).x(sinds) - drifts(tracks(tr).t(sinds),1);
    % ycorr = tracks(tr).y(sinds) - drifts(tracks(tr).t(sinds),2);
    % 
    % xcorr = smooth(xcorr,3);
    % ycorr = smooth(ycorr,3);
    % dx = diff(xcorr)*XYcal;
    % dy = diff(ycorr)*XYcal;
    dx = diff(tracks(tr).x(sinds))*XYcal;
    dy = diff(tracks(tr).y(sinds))*XYcal;
    
    dr = sqrt(dx.^2 + dy.^2);
    tts = ts(tracks(tr).t(sinds));
    dt = diff(tts);
    v = dr./dt*60;
    cvs = [cvs; v];
    avs = [avs; v];
end
ys = histcounts(cvs,eds);

plot(xs,ys,'LineWidth',2);
hold on

end

%% Simulations

eps = 5;
d = 0.7;
l = 7;
rho = 0.21;
v = 5;
Kagar = 250;
Kstiffs = 33;
rev = 8;
runs = 0;

fpath = simname(eps,d,l,rho,v,Kagar,Kstiff,rev,run);
[files, boxSize] = getframes(fpath);
asimvs = [];
for f = 1:20:numel(files)
    bds = loadsimdata(fullfile(files(f).folder,files(f).name));
    v = sqrt(([bds.vx].^2+[bds.vy].^2+[bds.vz].^2))';
    asimvs = [asimvs; v(rand(size(v))<0.1)];
end

%%

simys = histcounts(asimvs,eds,'Normalization','pdf');
exptys = histcounts(avs,eds,'Normalization','pdf');

plot(xs,exptys,'m','LineWidth',2);
hold on
plot(xs,simys,'b','LineWidth',2);
set(gca,'FontSize',24)

%% Reversal rates
eds = 1:3:60;
xs = 2.5:3:59;
arevs = [];
amidrevs = [];
aleadrevs = [];
atailrevs = [];
cols = ["#0072BD", "#D95319", "#EDB120", "#7E2F8E", "#77AC30", "#4DBEEE", "#A2142F"];

for f = 1:numel(fnames)
    load(fnames(f),'tracks','fpath');
    ts = getts(fpath);

    midrevs = [];
    leadrevs = [];
    tailrevs = [];

    for tr = 1:numel(tracks)

        trts = tracks(tr).t;
        trts = sort(trts);
        revts = ts([trts(1); tracks(tr).revs; trts(end)]);
        revts = sort(revts);
        times = diff(revts)/60;
        leadrevs = [leadrevs; times(1)];
        if (numel(times) == 1)
            arevs = [arevs; times(1) 1 1];
        else
            arevs = [arevs; times(1) 1 0];
            arevs = [arevs; times(end) 0 1];
        end
        if (numel(times)>2)
            for t = 2:numel(times)-1
                arevs = [arevs; times(t) 0 0];
            end
        end


        if numel(times)>1
            tailrevs = [tailrevs; times(end)];
        end
        times([1 end]) = [];
        midrevs = [midrevs; times];
    end
    % 
    % ys1 = histcounts(leadrevs,eds,'Normalization','pdf');
    % ys3 = histcounts(tailrecs,eds,'Normalization','pdf');
    % ys2 = histcounts(midrevs,eds,'Normalization','pdf');
    % plot(xs,ys1,'--','Color',cols(f),'LineWidth',2)
    % hold on
    % plot(xs,ys2,'Color',cols(f),'LineWidth',2)
    amidrevs = [amidrevs; midrevs];
    aleadrevs = [aleadrevs; leadrevs];
    atailrevs = [atailrevs; tailrevs];

end
%%
figure
histogram(aleadrevs,eds,'FaceColor','r');
hold on
histogram(atailrevs,eds,'FaceColor','b');
histogram(amidrevs,eds,'FaceColor','g');

%%

