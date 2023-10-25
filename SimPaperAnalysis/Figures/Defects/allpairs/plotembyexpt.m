kcwt = load("kcwt");
mbwt = load("mbwt");
mbpila = load("mbpila");
sim = load("sim");

%% 

pos = [0 0 500 500];
fig = figure('Position',pos);
ax1 = axes(fig,'Position',[0 0 0.25 0.25]);
ax2 = axes(fig,'Position',[0.25 0 0.25 0.25]);
ax3 = axes(fig,'Position',[0.5 0 0.25 0.25]);
ax4 = axes(fig,'Position',[0.75 0 0.25 0.25]);
ax5 = axes(fig,'Position',[0 0.25 0.25 0.25]);
ax6 = axes(fig,'Position',[0.25 0.25 0.25 0.25]);
ax7 = axes(fig,'Position',[0.5 0.25 0.25 0.25]);
ax8 = axes(fig,'Position',[0.75 0.25 0.25 0.25]);
ax9 = axes(fig,'Position',[0 0.5 0.25 0.25]);
ax10 = axes(fig,'Position',[0.25 0.5 0.25 0.25]);
ax11 = axes(fig,'Position',[0.5 0.5 0.25 0.25]);
ax12 = axes(fig,'Position',[0.75 0.5 0.25 0.25]);
ax13 = axes(fig,'Position',[0 0.75 0.25 0.25]);
ax14 = axes(fig,'Position',[0.25 0.75 0.25 0.25]);
ax15 = axes(fig,'Position',[0.5 0.75 0.25 0.25]);
ax16 = axes(fig,'Position',[0.75 0.75 0.25 0.25]);

I = imagesc(ax13,atan2(kcwt.nnmeanphisy,kcwt.nnmeanphisx)/3);
Adata = ones(size(kcwt.nnctx));
Adata(abs(kcwt.nnctx)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax9,kcwt.nnM2phisx./kcwt.nnctx);
I.AlphaData = Adata;
I = imagesc(ax5,kcwt.nnM2phisy./kcwt.nnctx);
I.AlphaData = Adata;
I = imagesc(ax1,kcwt.nnctx/sum(kcwt.nnctx,'all')/(0.5*0.5));
I.AlphaData = Adata;

I = imagesc(ax14,atan2(mbwt.nnmeanphisy,mbwt.nnmeanphisx)/3);
Adata = ones(size(mbwt.nnctx));
Adata(abs(mbwt.nnctx)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax10,mbwt.nnM2phisx./mbwt.nnctx);
I.AlphaData = Adata;
I = imagesc(ax6,mbwt.nnM2phisy./mbwt.nnctx);
I.AlphaData = Adata;
I = imagesc(ax2,mbwt.nnctx/sum(mbwt.nnctx,'all')/(0.5*0.5));
I.AlphaData = Adata;

I = imagesc(ax15,atan2(mbpila.nnmeanphisy,mbpila.nnmeanphisx)/3);
Adata = ones(size(mbpila.nnctx));
Adata(abs(mbpila.nnctx)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax11,mbpila.nnM2phisx./mbpila.nnctx);
I.AlphaData = Adata;
I = imagesc(ax7,mbpila.nnM2phisy./mbpila.nnctx);
I.AlphaData = Adata;
I = imagesc(ax3,mbpila.nnctx/sum(mbpila.nnctx,'all')/(0.5*0.5));
I.AlphaData = Adata;

I = imagesc(ax16,atan2(sim.nnmeanphisy,sim.nnmeanphisx)/3);
Adata = ones(size(sim.nnctx));
Adata(abs(sim.nnctx)<1) = 0;
I.AlphaData = Adata;
I = imagesc(ax12,sim.nnM2phisx./sim.nnctx);
I.AlphaData = Adata;
I = imagesc(ax8,sim.nnM2phisy./sim.nnctx);
I.AlphaData = Adata;
I = imagesc(ax4,sim.nnctx/sum(sim.nnctx,'all')/(0.5*0.5));
I.AlphaData = Adata;

axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,...
    ax15,ax16];
set(gcf,'Color','k')
for i = 1:numel(axs)
    set(axs(i),'Visible','off');
end

for i = 1:4
    set(axs(i),'clim',[0 4E-3]);
    axes(axs(i))
    colorcet('L3')
end
for i = 5:12
    set(axs(i),'clim',[0 0.8]);
    colormap(axs(i),ndefscmap);
end
for i = 13:16
    colormap(axs(i),colorcet('C9'));
end

set(ax13,'clim',[-pi/3 pi/3]);
set(ax14,'clim',[-pi/3 pi/3]);
set(ax15,'clim',[-pi/3 pi/3]);
set(ax16,'clim',[-pi/3 pi/3]);

% for i = [3, 4, 7, 8]
%     plot(axs(i),[40.5 40.5+7*cos(2*pi/3)],[40.5 40.5+7*sin(2*pi/3)],'w','LineWidth',3)
%     plot(axs(i),[40.5 40.5+7*cos(4*pi/3)],[40.5 40.5+7*sin(4*pi/3)],'w','LineWidth',3)
% 
% end

%%
