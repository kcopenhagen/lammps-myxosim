

fpath = '/scratch/gpfs/kc32/lammps-myxo-data/20230601/isof/rev8/run_000';

[files,boxSize] = getframes(fpath);

t = randi(numel(files));

bds = loadsimdata(fullfile(files(t).folder,files(t).name));
polf = pfield(bds,100,1000);
cim = cellim(bds,100,1000);
%%

fig = figure();
ax = axes(fig,'Position',[0 0 1 1],'Color','k');
im = imagesc(polf);
colormap orientcmap
im.AlphaData = (cim/2 + 0.25);
set(gcf,'Color','k')
set(gcf, 'InvertHardCopy', 'off'); 

axis equal
axis off
xc = 850;
yc = 950;
hold on
plot([0 100]+xc,[0 0]+yc,'w','LineWidth',3)