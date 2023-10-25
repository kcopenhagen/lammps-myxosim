kcwt = load("kcwt");
mbwt = load("mbwt");
mbpila = load("mbpila");
sim = load("sim");
%%

ppct = sim.ppctx;
pnct = sim.pnctx;
npct = sim.npctx;
nnct = sim.nnctx;

ppmeanphisx = sim.ppmeanphisx;
ppmeanphisy = sim.ppmeanphisy;
pnmeanphisx = sim.pnmeanphisx;
pnmeanphisy = sim.pnmeanphisy;
npmeanphisx = sim.npmeanphisx;
npmeanphisy = sim.npmeanphisy;
nnmeanphisx = sim.nnmeanphisx;
nnmeanphisy = sim.nnmeanphisy;

ppvarx = sim.ppM2phisx./sim.ppctx;
ppvary = sim.ppM2phisy./sim.ppctx;
pnvarx = sim.pnM2phisx./sim.pnctx;
pnvary = sim.pnM2phisy./sim.pnctx;
npvarx = sim.npM2phisx./sim.npctx;
npvary = sim.npM2phisy./sim.npctx;
nnvarx = sim.nnM2phisx./sim.nnctx;
nnvary = sim.nnM2phisy./sim.nnctx;

%%

ppct = kcwt.ppctx + mbwt.ppctx + mbpila.ppctx;
pnct = kcwt.pnctx + mbwt.pnctx + mbpila.pnctx;
npct = kcwt.npctx + mbwt.npctx + mbpila.npctx;
nnct = kcwt.nnctx + mbwt.nnctx + mbpila.nnctx;

ppmeanphisx = (kcwt.ppmeanphisx.*kcwt.ppctx+mbwt.ppmeanphisx.*mbwt.ppctx+...
    mbpila.ppmeanphisx.*mbpila.ppctx)./ppct;
ppmeanphisy = (kcwt.ppmeanphisy.*kcwt.ppctx+mbwt.ppmeanphisy.*mbwt.ppctx+...
    mbpila.ppmeanphisy.*mbpila.ppctx)./ppct;

pnmeanphisx = (kcwt.pnmeanphisx.*kcwt.pnctx+mbwt.pnmeanphisx.*mbwt.pnctx+...
    mbpila.pnmeanphisx.*mbpila.pnctx)./pnct;
pnmeanphisy = (kcwt.pnmeanphisy.*kcwt.pnctx+mbwt.pnmeanphisy.*mbwt.pnctx+...
    mbpila.pnmeanphisy.*mbpila.pnctx)./pnct;

npmeanphisx = (kcwt.npmeanphisx.*kcwt.npctx+mbwt.npmeanphisx.*mbwt.npctx+...
    mbpila.npmeanphisx.*mbpila.npctx)./npct;
npmeanphisy = (kcwt.npmeanphisy.*kcwt.npctx+mbwt.npmeanphisy.*mbwt.npctx+...
    mbpila.npmeanphisy.*mbpila.npctx)./npct;

nnmeanphisx = (kcwt.nnmeanphisx.*kcwt.nnctx+mbwt.nnmeanphisx.*mbwt.nnctx+...
    mbpila.nnmeanphisx.*mbpila.nnctx)./nnct;
nnmeanphisy = (kcwt.nnmeanphisy.*kcwt.nnctx+mbwt.nnmeanphisy.*mbwt.nnctx+...
    mbpila.nnmeanphisy.*mbpila.nnctx)./nnct;

ppvarx = (kcwt.ppM2phisx + mbwt.ppM2phisx + mbpila.ppM2phisx)./ppct;
pnvarx = (kcwt.pnM2phisx + mbwt.pnM2phisx + mbpila.pnM2phisx)./pnct;
npvarx = (kcwt.npM2phisx + mbwt.npM2phisx + mbpila.npM2phisx)./npct;
nnvarx = (kcwt.nnM2phisx + mbwt.nnM2phisx + mbpila.nnM2phisx)./nnct;

ppvary = (kcwt.ppM2phisy + mbwt.ppM2phisy + mbpila.ppM2phisy)./ppct;
pnvary = (kcwt.pnM2phisy + mbwt.pnM2phisy + mbpila.pnM2phisy)./pnct;
npvary = (kcwt.npM2phisy + mbwt.npM2phisy + mbpila.npM2phisy)./npct;
nnvary = (kcwt.nnM2phisy + mbwt.nnM2phisy + mbpila.nnM2phisy)./nnct;

%%

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

imagesc(ax1,ppct/sum(ppct,'all')/(0.5*0.5))
imagesc(ax2,pnct/sum(pnct,'all')/(0.5*0.5))
imagesc(ax3,npct/sum(npct,'all')/(0.5*0.5))
imagesc(ax4,nnct/sum(nnct,'all')/(0.5*0.5))

I = imagesc(ax5,atan2(ppmeanphisy,ppmeanphisx));
Adata = ones(size(ppct));
Adata(abs(ppct)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax6,atan2(pnmeanphisy,pnmeanphisx)/3);
Adata = ones(size(pnct));
Adata(abs(pnct)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax7,atan2(npmeanphisy,npmeanphisx));
Adata = ones(size(npct));
Adata(abs(npct)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax8,atan2(nnmeanphisy,nnmeanphisx)/3);
Adata = ones(size(nnct));
Adata(abs(nnct)<1) = 0;
I.AlphaData = Adata;
axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8];
set(gcf,'Color','k')
for i = 1:4
    set(axs(i),'clim',[0 4E-3]);
    axes(axs(i))
    colorcet('L3')
end
for i = 5:8
    colormap(axs(i),colorcet('C9'));
end
for i = 1:8
    set(axs(i),'Visible','off');
    hold(axs(i),'on');
%    plot(axs(i),40.5,40.5,'w.','MarkerSize',30)
%    plot(axs(i),[40.5 47.5],[40.5 40.5],'w','LineWidth',3)
end

set(ax5,'clim',[-pi pi]);
set(ax6,'clim',[-pi/3 pi/3]);
set(ax7,'clim',[-pi pi]);
set(ax8,'clim',[-pi/3 pi/3]);

% for i = [3, 4, 7, 8]
%     plot(axs(i),[40.5 40.5+7*cos(2*pi/3)],[40.5 40.5+7*sin(2*pi/3)],'w','LineWidth',3)
%     plot(axs(i),[40.5 40.5+7*cos(4*pi/3)],[40.5 40.5+7*sin(4*pi/3)],'w','LineWidth',3)
% 
% end

%%
