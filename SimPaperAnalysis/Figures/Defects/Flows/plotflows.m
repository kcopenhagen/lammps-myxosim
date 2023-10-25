kcwt = load("kcwt.mat");
mbwt = load("mbwtt.mat");
mbpila = load("mbflavot.mat");
sim = load("simt.mat");

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

sz = size(kcwt.pvx);
x = 1:sz;
y = x;
[x,y] = meshgrid(x,y);
dx = 5;
xinds = (sz(2) + 1)/2-dx*round(sz(2)/dx/2):dx:(sz(2) + 1)/2+dx*round(sz(2)/dx/2);
yinds = xinds;
[xx,yy] = meshgrid(xinds,yinds);
inds = sub2ind(sz,yy,xx);

sz = size(sim.pvx);
sx = 1:sz;
sy = sx;
[sx,sy] = meshgrid(sx,sy);
dx = 5;
xinds = (sz(2) + 1)/2-dx*floor(sz(2)/dx/2):dx:(sz(2) + 1)/2+dx*floor(sz(2)/dx/2);
yinds = xinds;
[xx,yy] = meshgrid(xinds,yinds);
sinds = sub2ind(sz,yy,xx);


imagesc(ax1,sqrt(kcwt.nvx.^2+kcwt.nvy.^2));
hold(ax1,'on');
quiver(ax1,x(inds),y(inds),kcwt.nvx(inds),kcwt.nvy(inds),'w');
imagesc(ax2,sqrt(mbwt.nvx.^2+mbwt.nvy.^2));
hold(ax2,'on');
quiver(ax2,x(inds),y(inds),mbwt.nvx(inds),mbwt.nvy(inds),'w');
imagesc(ax3,sqrt(mbpila.nvx.^2+mbpila.nvy.^2));
hold(ax3,'on');
quiver(ax3,x(inds),y(inds),mbpila.nvx(inds),mbpila.nvy(inds),'w');
imagesc(ax4,sqrt(sim.nvx.^2+sim.nvy.^2));
hold(ax4,'on');
quiver(ax4,sx(sinds),sy(sinds),sim.nvx(sinds),sim.nvy(sinds),'w');
imagesc(ax5,sqrt(kcwt.pvx.^2+kcwt.pvy.^2));
hold(ax5,'on');
quiver(ax5,x(inds),y(inds),kcwt.pvx(inds),kcwt.pvy(inds),'w');
imagesc(ax6,sqrt(mbwt.pvx.^2+mbwt.pvy.^2));
hold(ax6,'on');
quiver(ax6,x(inds),y(inds),mbwt.pvx(inds),mbwt.pvy(inds),'w');
imagesc(ax7,sqrt(mbpila.pvx.^2+mbpila.pvy.^2));
hold(ax7,'on');
quiver(ax7,x(inds),y(inds),mbpila.pvx(inds),mbpila.pvy(inds),'w');
imagesc(ax8,sqrt(sim.pvx.^2+sim.pvy.^2));
hold(ax8,'on');
quiver(ax8,sx(sinds),sy(sinds),sim.pvx(sinds),sim.pvy(sinds),'w');

axs = [ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8];

for i = 1:numel(axs)
    set(axs(i),'Visible','off');
end
for i = 1:4
    set(axs(i),'Colormap',magma);%,'clim',[0 0.3]);
end
for i = 5:8
    set(axs(i),'Colormap',magma);%,'clim',[0 0.6]);
end
