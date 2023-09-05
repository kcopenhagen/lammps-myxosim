Ksprs = 33:33:231;
runs = 0:4;
fpaths = '/scratch/gpfs/kc32/lammps-myxo-data/20230701';
dr = 1;
eds = 0:dr:100;
rs = dr/2:dr:eds(end);
pncols = copper(numel(Ksprs));
for k = 1:numel(Ksprs)
    Ksprs(k)
    ppdrs = [];
    nndrs = [];
    pndrs = [];
    for runn = 1:numel(runs)
        fpath = fullfile(fpaths,sprintf('Kspr%d',Ksprs(k)),sprintf('run_00%d',runs(runn)));
        load(fullfile(fpath,'adefs.mat'));
        ts = unique([adefs.t]);
        for t = 1:5:numel(ts)
            cdefs = adefs([adefs.t]==ts(t));
            pdefs = cdefs([cdefs.q]>0);
            ndefs = cdefs([cdefs.q]<0);
            
            ppdx = [pdefs.x] - [pdefs.x]';
            ppdy = [pdefs.y] - [pdefs.y]';
            ppdr = sqrt(ppdx.^2+ppdy.^2);
            ppdrs = [ppdrs; ppdr(ppdr>0)];

            nndx = [ndefs.x] - [ndefs.x]';
            nndy = [ndefs.y] - [ndefs.y]';
            nndr = sqrt(nndx.^2+nndy.^2);
            nndrs = [nndrs; nndr(nndr>0)];

            pndx = [pdefs.x] - [ndefs.x]';
            pndy = [pdefs.y] - [ndefs.y]';
            pndr = sqrt(pndx.^2+pndy.^2);
            pndrs = [pndrs; pndr(pndr>0)];
        end
        
    end
    pndrct = histcounts(pndrs,eds);
    plot(rs,-log(pndrct./rs),'LineWidth',2,'Color',pncols(k,:));

    hold on
    drawnow
end