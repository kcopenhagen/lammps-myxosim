%% Height histograms.

%% Experiments

exptnames = ["190109KC3", "190109KC3b", "190111KC2", "190116KC2", "190117KC1", "190117KC2", "190117KC2a", "190117KC2b"]; 
eds = 0.5:0.01:1;
fracpiled = [];

for f = 1:numel(exptnames)
    fpath = (fullfile('/tigress/kc32/Monolayerpaperdata/Data',exptnames(f)));
    ts = getts(fpath);
    ct = 0;
    tests = 1:numel(ts);
    temp = [];
    for t = 1:numel(tests)
        try
            heights = loaddata(char(fpath),tests(t),'covid_layers','int8');
            temp = [temp; sum(heights==1, 'all')/numel(heights)];
            fracpiled = [fracpiled; sum(heights==1, 'all')/numel(heights)];
        end
    end
end

%% Sim heights

eps = 5;
d = 0.7;
l = 7;
rho = 0.21;
v = 5;
Kagar = 100:50:500;
Kstiff = 15;
rev = 8;
runs = 0:3;

res = 100;
type = "beads";
props = [];
for k = 1:numel(Kagar)
    ashs = [];
    for r = 1:numel(runs)
        fpath = simname(eps,d,l,rho,v,Kagar(k),Kstiff,rev,runs(r));
        [files, boxSize] = getframes(fpath);
        for t = 5:50:numel(files)-5
            if type == "him"
                bdst = [];
                for tt = t
                    bdst = [bdst loadsimdata(fullfile(files(tt).folder,files(tt).name))];
                end
                hs = heightim(bdst,boxSize,res,1);
                hs = medfilt2(hs,[3 3]);
                ashs = [ashs; (sum(abs(hs(:) - 0.3) < 0.15))/numel(hs(:))];
            elseif type == "beads"
                bds = loadsimdata(fullfile(files(t).folder,files(t).name));
                mols = unique([bds.mol]);
                nmono = 0;
                ntot = 0;
                for i = 1:numel(mols)
                    ccell = bds([bds.mol]==mols(i));
                    nmonot = sum(abs([ccell.zs] - 0.3) < 0.15);
                    ntott = numel(ccell);
                    if (nmonot == 1)
                        nmonot = nmonot - 1;
                        ntott = ntott - 1;
                    end
                    nmono = nmono + nmonot;
                    ntot = ntot + ntott;
                end
                % ashs = [ashs; (sum(abs([bds.zs] - 0.3)<0.15))/numel(bds)];
                ashs = [ashs; nmono / ntot];
            end
        end
    end
    propt = struct('Kagar',Kagar(k),'hs',ashs);
    props = [props; propt];
end

%%

eds = 0.5:0.02:1;
xs = (eds(1)+(eds(2)-eds(1))/2):(eds(2)-eds(1)):eds(end);

cols = sky(numel(props));
for k = 1:numel(props)
    ys = histcounts(props(k).hs,eds,'Normalization','pdf');
    plot(xs,ys,'LineWidth',2,'Color',cols(k,:));
    hold on
end
ys = histcounts(fracpiled,eds,'Normalization','pdf');
plot(xs,ys,'LineWidth',2,'Color','m')
set(gca,'FontSize',24)
colormap sky
c = colorbar;
c.FontSize = 24;
%c.Title.String = "K_{agar}";
c.TickLabels = [100 300 500];